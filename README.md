# PrecisionPhage

Scientific code and compact data for leakage-controlled phage–host prediction
and retrospective scoring of published, fixed-bank *Escherichia* cocktail
assays. This is a research audit, **not a clinical recommender**.

The dated scientific snapshot is `v2026.10.15.1`. It extends the existing public
v2.0.1 history without replacing historical frozen predictions. Source code,
saved predictions, analysis protocols, accession provenance and verification
helpers are included. Manuscripts, private peer review, submission materials,
paper-writing code and large genome files are not part of this release.

## What the evidence supports

| Analysis | Evaluated data | Main finding and boundary |
| --- | --- | --- |
| Full interaction-table audit | 8,849 experimentally assayed pairs from three source studies | Cross-source performance depends on source and representation; the four-feature neural baseline is retained alongside the tree baseline. |
| Within-source sequence evaluation | 1,947 covered NCBI_HR pairs; 1,418 phages and 323 host reference genomes | Saved tree AUROC is 0.959 in eligible species holdouts and 0.785 in five dual-exclusion blocks. Cohorts differ; later matched-cohort and canonical-distance sensitivities are explicitly retrospective. |
| Frozen NCBI_HR → StaphStudy | 1,053 pairs, positive prevalence 0.316 | Full-input AUROC 0.637, two-way sequence-cluster interval 0.498–0.739; AP 0.391 and ECE 0.499. Broad transfer and portable calibrated probabilities are not established. |
| Frozen sequence-only → BASEL–ECOR | 1,240 pairs, positive prevalence 0.244 | AUROC 0.470 and AP 0.222. This further public-data transfer result is negative; the full model's four source-supplied columns were unavailable. |
| Historical fixed-bank cocktail scoring | 100 source-selected host–triplet records, 82 active outcomes, 18 distinct compositions | AUROC 0.807 and AP 0.951 against prevalence 0.820. Source authors chose the triplets; this does not validate de novo cocktail selection. |

The cocktail point estimates are unchanged. The original host-only AUROC
interval is 0.693–0.906. A **new retrospective crossed host-cluster × composition
sensitivity** gives 0.543–0.965, AP interval 0.852–0.998, and AP-minus-prevalence
0.131 with interval 0.024–0.257. These wider intervals address repeated menus
under an additional resampling assumption; they do not establish transport to
new banks, compositions or selection procedures. Distinct menus can also share
phages, a dependence not fully covered by whole-menu grouping.

No resolved advantage over simple constituent aggregation or the stronger
matched controls is claimed. The measured two-option comparison has only 17
hosts with different binary outcomes. Its historical policy gain over the one
fixed baseline is 0.050; the new crossed interval is **−0.039–0.194**, leaving
selection benefit unresolved.

The original-source comparison is complete for **64 hosts only**. Of 300
source-selected constituent records, 233 have released numeric decision scores;
67 rescue-stage values are absent. Inspection of the source release and all 96
prediction-archive members did not recover authentic external scores or a
faithfully runnable all-100 predictor. Missing scores were not invented. See
`data/results_v2/frontiers_revision_2026_10_15/official_gaborieau/` for the public
feasibility record.

## Quick start: audit supplied predictions

Use Python 3.11. The historical `requirements-lock.txt` is distinct from the
actual revision-audit environment in
`reproducibility/requirements-revision-audit.txt`. The audit requirements include
sourmash 4.9.3 for canonical-sketch tests and openpyxl 3.1.5 for public
accession-table preparation. The `audit` package extra also records sourmash.

```powershell
git clone https://github.com/JaivardhanSamadhiya/frontiers_submission.git
cd frontiers_submission
git checkout v2026.10.15.1
python -m venv .venv-audit
.venv-audit\Scripts\Activate.ps1
python -m pip install -r reproducibility/requirements-revision-audit.txt
python -m unittest discover -s tests -v
$env:PYTHONPATH = "src"
python -m pytest tests/test_screening_diagnostics.py -q
python experiments/00_validate_release.py --allow-missing-genomes
python reproducibility/verify_public_snapshot.py
python reproducibility/verify_saved_predictions.py
```

Linux/macOS users activate `.venv-audit/bin/activate` and run the separate
screening function tests with
`PYTHONPATH=src python -m pytest tests/test_screening_diagnostics.py -q`.
The suites cover 38 unittest cases plus four screening function tests.
No GPU or paid service is needed for the saved-score checks. The tests do not replace scientific
review or certify biological usefulness. The legacy release validator checks
historical artifacts; the public snapshot inventory checks this dated payload.

## Reproduce revision analyses

The following commands use saved predictions and compact tables. They do not
refit a predictive model:

```text
python scripts/frontiers_revision/submission_audit.py
python scripts/frontiers_revision/cocktail_reporting.py --self-test-only
python scripts/frontiers_revision/cocktail_reporting.py
python scripts/frontiers_revision/crossed_cocktail_uncertainty.py
python scripts/frontiers_revision/crossed_two_options.py
```

Reruns regenerate outputs in their named analysis directories. Run them in a
separate clone if preserving a byte-identical release checkout is important;
timestamps and serialized metadata may change. Compare numerical predictions
and input checksums, not newly generated metadata bytes alone.

The exact phage-alias check additionally needs two small public sequence
records. Fetch only the pinned, checksum-checked sources before rerunning it:

```text
python reproducibility/stage_public_sources.py
python scripts/frontiers_revision/cocktail_source_comparison.py
```

The source-predictor feasibility audit is retained as a saved non-success
record. Rerunning `official_gaborieau_audit.py` requires a full pristine upstream
Git checkout at its pinned revision, including the released prediction RAR,
not merely the compact source CSVs supplied here. This optional forensic rerun
also needs seaborn for the upstream script and a RAR-capable libarchive/bsdtar
executable available as `tar`; GNU tar alone cannot read that archive. It performs no model fitting
and does not fill missing rescue scores. Public summary records omit contact
drafts, personal commit metadata and private machine-path details.

## Files and stage separation

- `src/precisionphage/`: active v2 model, feature, split and evaluation code.
- `configs/default.yaml`: frozen analysis definitions; not a claim of an
  optimally selected model. The original GBM selection rationale is incompletely
  recoverable from surviving development records.
- `data/results_v2/external_cocktail_validation/`: original staged predictions
  and evaluation artifacts; historical values remain unchanged.
- `data/results_v2/frontiers_revision_2026_10_14/`: corrected cohorts, coverage,
  controls, source-score comparisons, canonical-cutoff/feature ablations,
  target-label recalibration and genome provenance.
- `data/results_v2/frontiers_revision_2026_10_15/`: crossed-dependence checks and
  official-source comparator recovery attempt.
- `reproducibility/`: source notices, pinned input retrieval, actual audit
  environment, saved-prediction checks and the public-payload checksum inventory.
- `DATA.md` and `DESIGN.md`: source provenance and analysis definitions, with
  dated corrections identifying superseded historical interpretations.

Historical row-independent p-values, RaFAH-inspired proxies, mathematical set
cover and unfitted resistance outputs may remain in the retained scientific
history. They are **not evidence for the current contribution**, released-tool
superiority, calibrated clinical probabilities or resistance suppression.
Top-level v1 scripts and archived v1 data are not the active analysis.

## Genomes and raw-data limits

Large FASTAs are intentionally absent. The frozen NCBI_HR bundle has versioned
header/accession/order/length/sequence-hash records for all 1,418 phage FASTAs and
323 host FASTAs (666 ordered host records, 656 distinct accessions). Its original
file hashes matched the preserved inventory. A fresh download does not
automatically recreate identical FASTA bytes. Use the opt-in recovery helper
described in `DATA.md`; it refuses mismatches and does not overwrite originals.

Optional raw MinHash caches are omitted for size; complete cached distance
matrices and their source parameters remain available. Genome-dependent raw
feature regeneration and model refitting require the verified original inputs
or a separately versioned reconstructed sequence release. The Nahant feasibility
audit is **cache-limited**: zero pairs had both genomes in the checked caches,
not proof that those genomes are unavailable everywhere.

The full ECOR source-article XML is also not redistributed. Its primary URL and
historical checksum remain recorded; `prepare_basel_ecor_genomes.py` retrieves
it if needed. Derived accession mappings are supplied without the article's
correspondence details.

Historical source/cocktail protocols describe local staged locks. They are not
independent prospective preregistrations, and public outcomes had been examined
in other development analyses. Later controls and repairs are retrospective;
rerunning them cannot recreate original outcome blindness.

## Sources, licensing and citation

Code is MIT licensed. Third-party data and tools retain their own terms; see
`LICENSE` and `reproducibility/SOURCE_NOTICES.md`. Main data sources are
[VHIP](https://doi.org/10.1371/journal.pcbi.1011649),
[Gaborieau et al.](https://doi.org/10.1038/s41564-024-01832-5), and
[GenoPHI](https://doi.org/10.1038/s41564-026-02482-5). The PHIST submodule remains
pinned but is not needed for saved-prediction checks. Its source may be obtained
with `git submodule update --init --recursive` when rebuilding that comparator.

Use `CITATION.cff` and cite the exact tag or commit used. `.zenodo.json` supplies
archive metadata but **does not itself create a deposit or DOI**. No unissued
archival identifier is claimed. Scientific snapshot version and historical
dataset/model version are distinct.

Patch `v2026.10.15.1` corrects archive portability of the checksum verifier by
recording exact staged Git-blob hashes alongside working-tree hashes. It does
not change scientific inputs, predictions or results; `v2026.10.15` is retained.
