# PrecisionPhage

Leakage-controlled phage–host interaction prediction, phage-cocktail optimization, and resistance-aware therapy simulation.

**Manuscript:** All numbers, tables, and figures in the Frontiers manuscript come from `data/results_v2/`. The older `data/archive_results_v1/` directory contains development experiments from the v1 flat pipeline and should **not** be used to reproduce the manuscript.

## Quick start (v2 pipeline)

```bash
pip install -e .
pip install python-docx   # for manuscript generation only

# Run the numbered experiment scripts in order (see experiments/)
python experiments/01_build_dataset.py
python experiments/04_leakage_splits.py
# ... through experiments/15_write_paper.py
```

## Repository layout

| Path | Purpose |
|------|---------|
| `experiments/` | Numbered v2 pipeline scripts |
| `src/precisionphage/` | Package implementation |
| `configs/default.yaml` | All v2 hyperparameters and paths |
| `data/results_v2/` | **Canonical manuscript results** (CSVs, JSON, figures) |
| `data/interim_v2/` | Cleaned datasets and feature caches |
| `data/archive_results_v1/` | Deprecated v1 outputs (do not cite) |
| `data/raw/` | Manual inputs (VHI CSV, host FASTAs, etc.) |

## Manual inputs required

- `data/raw/VirusHostInter.csv`
- `data/raw/phage_ncbi_refseq_def_info.txt`
- `data/fastas/hosts/` (one `.fasta` per host species)

See [DATA.md](DATA.md) for the full data policy.

## Reproducing figures

Manuscript figures live in `data/results_v2/`:

| Figure | File | Generator |
|--------|------|-----------|
| 1 | `figure_main.png` | `experiments/10_report.py` |
| 2 | `fig_feature_importance.png` | `experiments/11_figures.py` or `16_regenerate_manuscript_figures.py` |
| 3 | `fig_phist_compare.png` | `experiments/16_regenerate_manuscript_figures.py` |
| 4 | `cocktail_coverage.png` | `experiments/08_cocktail.py` or `16_regenerate_manuscript_figures.py` |
| 5 | `temporal_dynamics.png` | `experiments/09_temporal.py` then `16_regenerate_manuscript_figures.py` |

Regenerate all five from saved result tables (no model retraining):

```bash
python experiments/16_regenerate_manuscript_figures.py
python experiments/15_write_paper.py   # re-embed figures into the Word manuscript
```

Figure 5 additionally requires `temporal_trajectory.npz`, written by `experiments/09_temporal.py` (needs host FASTAs).

## Host genomes (required for full pipeline)

Host reference FASTAs are **not** stored in git (~1.4 GB). Fetch them once:

```bash
# Edit configs/default.yaml -> fetch.ncbi_email with your NCBI email
python experiments/fetch_host_genomes.py
```

This writes `data/fastas/hosts/*.fasta` and `data/raw/host_genome_manifest.json`.

## GitHub push checklist

1. Set `fetch.ncbi_email` locally; do **not** commit secrets.
2. Run `python experiments/16_regenerate_manuscript_figures.py` and confirm five PNGs in `data/results_v2/`.
3. Commit: `src/`, `experiments/`, `configs/`, `data/results_v2/` (tables + figures), `README.md`, `DATA.md`, `environment.yml`, `pyproject.toml`.
4. Exclude (via `.gitignore`): `data/phages/`, `data/fastas/hosts/`, `.venv/`.
5. Push to `https://github.com/JaivardhanSamadhiya/frontiers_submission` and verify figures render on GitHub.

## Environment

See `environment.yml` or install from `pyproject.toml`:

```bash
pip install -e ".[extras]"
```

Python ≥ 3.10; PyTorch and PyTorch Geometric required for GNN experiments.
