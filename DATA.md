# Data policy (PrecisionPhage)

## One command does everything

```bash
cd SciFair2026
pip install -r scripts/requirements.txt
python scripts/run_pipeline.py
```

Pipeline steps: **verify inputs → fetch online data → assemble dataset → features → models → stats → ablations → plots → summary docs**.

Use `--skip-download` only when offline and `virushostdb_raw.tsv` + `phagesdb_records_cache.json` are already cached.

---

## Manual inputs (you provide — never downloaded)

| File / folder | Purpose |
|---------------|---------|
| `data/raw/VirusHostInter.csv` | Millard lab VHI interaction table |
| `data/raw/phage_ncbi_refseq_def_info.txt` | NCBI RefSeq metadata (accession list for phage downloads) |
| `data/fastas/hosts/*.fasta` | Host reference genomes (~415 species, ~1.4 GB). Fetch with `python experiments/fetch_host_genomes.py` (set `fetch.ncbi_email` in `configs/default.yaml`) or copy from your local archive. |
| `data/raw/saureus_*_interactions.csv` | Optional; auto-derived from VHI if absent |

These are checked at startup by `fetch_data.verify_manual_inputs()`.

---

## Online sources (dedicated functions in `scripts/fetch_data.py`)

| Function | Output | Source |
|----------|--------|--------|
| `fetch_virushostdb()` | `data/raw/virushostdb_raw.tsv` | KEGG FTP |
| `fetch_phagesdb_cache()` | `data/raw/phagesdb_records_cache.json` | PhagesDB REST API |
| `fetch_phage_fastas()` | `data/phages/NC_*.fasta` | NCBI Entrez |
| `fetch_inphared_tsv()` | `data/raw/inphared/*.tsv` | Optional URL in `config.yaml` |
| `ensure_validation_csvs()` | `data/raw/saureus_*.csv` | Derived from VHI (no network) |

`data_collection.py` **only reads local files** — it never downloads.

---

## Generated on each run (gitignored)

- `data/results/` — metrics, frozen dataset, logs, markdown summaries
- `data/features/` — feature CSVs
- `data/plots/` — all publication figures + `captions.md`

---

## Optional INPHARED

Place a `*_data_excluding_refseq.tsv` in `data/raw/inphared/`, or set `fetch.inphared_tsv_url` in `config.yaml` to a direct download link from the [Millard lab INPHARED page](https://millardlab.org/bacteriophage-genomics/inphared/).
