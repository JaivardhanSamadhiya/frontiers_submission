# PrecisionPhage

PrecisionPhage is a computational study of phage-host interaction prediction,
phage-cocktail set cover, and an illustrative resistance sensitivity model. The
active implementation is the v2 package under `src/precisionphage`; top-level
v1 scripts and archived v1 outputs are not part of the reported analysis.

## What the frozen analysis contains

- 8,849 experimentally assayed pairs from the VHIP input table: 2,770 positive
  and 6,079 negative, across 2,331 phages and 387 host identifiers.
- A sequence-covered subset of 1,947 NCBI_HR pairs: 1,488 positive and 459
  negative, with 1,418 phage and 323 host reference genomes.
- A locked, retraining-free test of the full sequence pipeline on all 1,053
  sequence-covered StaphStudy pairs (333 positive and 720 negative).
- Taxonomic, in-house sequence-cluster, and dual cold-start evaluations.
- XGBoost, GraphSAGE/edge-decoder, PHIST, and a RaFAH-inspired proxy.
- Host-cluster-grouped out-of-fold predictions for exploratory set cover.
- A deterministic ODE sensitivity analysis. It is not wet-lab validation and
  does not estimate clinical efficacy.

The three studies are present in the full four-feature baseline. The frozen
24-feature GBM was subsequently trained on NCBI_HR and applied without
retraining to StaphStudy: AUROC 0.637 (two-way sequence-cluster bootstrap 95%
CI 0.498-0.739), AUPRC 0.391, and ECE 0.499. This is source-held-out sequence
validation, but it is not a prospectively blinded or independently collected
cohort because StaphStudy was already represented in VHIP and in the earlier
four-feature audit. GNN, comparator, cocktail, and temporal analyses remain
single-source. See `EXTERNAL_VALIDATION_PROTOCOL.md`,
`data/results_v2/external_staph_validation_summary.csv`, and `DATA.md`.

## Reproducible setup

Clone with the PHIST submodule and install the locked environment:

```bash
git clone --recurse-submodules <repository-url>
cd full_pipeline
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements-lock.txt
python -m pip install -e . --no-deps
```

## Data acquisition

Large genome files are intentionally not stored in Git. The compact frozen
interaction snapshot at `data/raw/VirusHostInter.csv` comes from Supporting
Table S2 of Bastien et al. (2024):

- Article and data availability: https://doi.org/10.1371/journal.pcbi.1011649
- S2 interaction/model-input table: https://doi.org/10.1371/journal.pcbi.1011649.s010
- S5 viral and host accession list: https://doi.org/10.1371/journal.pcbi.1011649.s013
- VHIP analysis repository: https://github.com/DuhaimeLab/VHIP_analyses_Bastien_et_al_2023

The 1,418 phage and 323 host FASTAs used for the frozen sequence analysis were
retrieved from NCBI records represented in the VHIP accessions and staged under
`external/phist_run/phages` and `external/phist_run/hosts`. Those FASTAs are
excluded from this repository because of their size. The compact staging maps
remain in `external/phist_run/phage_map.json` and `host_map.json`.

The existing fetch scripts obtain optional Virus-Host DB, PhagesDB, Nahant, and
host caches. `experiments/fetch_staph_validation_sources.py` stages the pinned
upstream VHIP repositories that contain the StaphStudy example genomes and
identifier crosswalk; experiment 26 then creates the 39 analysis aliases. The
script does **not** recreate the separate NCBI_HR training FASTA bundle. A fresh
NCBI download should be treated as a new dataset version and checksummed.

For a code-only clone without the genome bundle, run:

```bash
python experiments/00_validate_release.py --allow-missing-genomes
python -m unittest discover -s tests -v
```

After staging the exact FASTAs locally, omit `--allow-missing-genomes` to run
the strict coverage audit.

Then run numbered experiments in order. Steps 01-14 train models or compute
analyses; step 16 redraws the primary figures from saved result artifacts.

```bash
python experiments/01_build_dataset.py
python experiments/02_baseline.py
# continue through experiments/14_external_report.py
python experiments/16_regenerate_figures.py
python experiments/fetch_staph_validation_sources.py
python experiments/25_external_sequence_overlap_audit.py
python experiments/26_frozen_cross_source_validation.py
```

The manuscript source and submission package are intentionally not published in
this code repository. The frozen result files live under `data/results_v2`. Do
not combine them with `data/archive_results_v1`.

After rerunning the audited analyses, maintainers can create the revised DOCX
and submission-format figures with:

```bash
python experiments/17_revise_manuscript.py
python experiments/18_package_submission.py
```

These manuscript utilities require the author's local, gitignored source DOCX;
they are included for auditability but are not part of the code-only CI run.

## Important interpretation limits

- Host identifiers in the sequence-covered subset are species/taxon labels, not
  a strain-resolved therapeutic panel.
- Only experimentally tested pairs are labeled. Missing matrix cells are
  unobserved, not confirmed negatives.
- The custom bottom-k MinHash clustering uses a Mash-form distance; it is not a
  run of Mash and its threshold is not a direct ANI measurement.
- The fixed XGBoost configuration is the headline baseline. Nested tuning was
  evaluated separately and was not used for the saved headline estimates.
- Row-wise DeLong, McNemar, permutation, and bootstrap results are exploratory
  because observations share phage and host entities.
- The locked NCBI_HR-to-StaphStudy test supports modest but uncertain
  cross-source ranking, not portable calibration or prospective validation.
- The RaFAH-style result is an in-house protein-feature proxy, not the published
  pretrained RaFAH model.
- The temporal model uses unfitted parameters to illustrate complete
  cross-resistance, partial dependence, and independent resistance. Its grid
  outputs are conditional equation behavior, not biological prediction or
  efficacy validation.

## Layout

| Path | Purpose |
|---|---|
| `configs/default.yaml` | Frozen v2 paths and parameters |
| `src/precisionphage/` | Active implementation |
| `experiments/` | Numbered analysis drivers |
| `tests/` | Release and leakage-invariant tests |
| `external/phist_run/` | Compact PHIST staging maps; large FASTAs are local-only |
| `data/interim_v2/` | Canonical cleaned data and feature caches |
| `data/results_v2/` | Saved tables, compact predictions, and figures |

## Hardware

The reported pipeline ran on CPU. A GPU can accelerate GraphSAGE, but it is not
required. An AMD ROCm environment must be validated against the selected
PyTorch/PyG wheels before paid use.
