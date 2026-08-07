# PrecisionPhage

Leakage-controlled phage–host interaction prediction, phage-cocktail optimization, and resistance-aware therapy simulation.

All headline numbers, tables, and figures are produced under `data/results_v2/`. The older `data/archive_results_v1/` directory contains development experiments from the v1 flat pipeline and should **not** be used.

## Quick start (v2 pipeline)

```bash
pip install -e .

# Run the numbered experiment scripts in order (see experiments/)
python experiments/01_build_dataset.py
python experiments/04_leakage_splits.py
# ... through experiments/14_external_report.py
```

## Repository layout

| Path | Purpose |
|------|---------|
| `experiments/` | Numbered v2 pipeline scripts |
| `src/precisionphage/` | Package implementation |
| `configs/default.yaml` | All v2 hyperparameters and paths |
| `data/results_v2/` | **Canonical results** (CSVs, JSON, figures) |
| `data/interim_v2/` | Cleaned datasets and feature caches |
| `data/archive_results_v1/` | Deprecated v1 outputs (do not cite) |
| `data/raw/` | Manual inputs (VHI CSV, host FASTAs, etc.) |

## Manual inputs required

- `data/raw/VirusHostInter.csv`
- `data/raw/phage_ncbi_refseq_def_info.txt`
- `data/fastas/hosts/` (one `.fasta` per host species)

See [DATA.md](DATA.md) for the full data policy.

## Reproducing figures

Result figures live in `data/results_v2/`:

| Panel | File | Generator |
|-------|------|-----------|
| Main (3-panel) | `figure_main.png` | `experiments/10_report.py` |
| Feature importance | `fig_feature_importance.png` | `experiments/11_figures.py` or `16_regenerate_figures.py` |
| External baselines | `fig_phist_compare.png` | `experiments/16_regenerate_figures.py` |
| Cocktail coverage | `cocktail_coverage.png` | `experiments/08_cocktail.py` or `16_regenerate_figures.py` |
| Therapy simulation | `temporal_dynamics.png` | `experiments/09_temporal.py` then `16_regenerate_figures.py` |

Regenerate all five from saved result tables (no model retraining):

```bash
python experiments/16_regenerate_figures.py
```

Figure 5 additionally requires `temporal_trajectory.npz`, written by `experiments/09_temporal.py` (needs host FASTAs).

## Host genomes (required for full pipeline)

Host reference FASTAs are **not** stored in git (~1.4 GB). Fetch them once:

```bash
# Edit configs/default.yaml -> fetch.ncbi_email with your NCBI email
python experiments/fetch_host_genomes.py
```

This writes `data/fastas/hosts/*.fasta` and `data/raw/host_genome_manifest.json`.

## Environment

See `environment.yml` or install from `pyproject.toml`:

```bash
pip install -e ".[extras]"
```

Python ≥ 3.10; PyTorch and PyTorch Geometric required for GNN experiments.
