"""
summary.py - Methods skeleton and reproducibility card
=======================================================
generate_summary_documents(...) writes:
    results/methods_skeleton.md
    results/reproducibility_card.md
    results/final_summary.csv

All numeric values are filled in from the live results dict, so a reviewer
can copy the methods text directly into a paper.
"""
from __future__ import annotations

import logging
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def _summary_table(all_results: dict) -> pd.DataFrame:
    rows = []
    for name, res in all_results.items():
        rows.append({
            "model": name,
            "loso_mean_auc": round(float(res.get("loso_mean",
                                                   res.get("mean_auc", 0))), 4),
            "loso_std_auc": round(float(res.get("loso_std", 0)), 4),
            "loso_pooled_auc": round(float(res.get("loso_pooled", 0)), 4),
            "logo_mean_auc": round(float(res.get("logo_mean", 0)), 4)
                if "logo_mean" in res else None,
            "mc_unseen_auc": round(float(res.get("mc_auc", 0)), 4)
                if "mc_auc" in res else None,
        })
    return pd.DataFrame(rows)


def generate_summary_documents(dataset: pd.DataFrame, all_results: dict,
                               stat_results: dict, ablation_results: dict,
                               name_ablation: dict, exp_results,
                               valid_species: list[str],
                               feature_meta: dict, run_meta: dict,
                               sha256: str, cfg: dict) -> None:
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)

    # Final summary CSV
    summary_df = _summary_table(all_results)
    summary_df.to_csv(results_dir / "final_summary.csv", index=False)

    # Dataset stats
    n_total = int(len(dataset))
    n_pos = int((dataset["label"] == 1).sum())
    n_neg = int((dataset["label"] == 0).sum())
    n_phages = int(dataset["phage"].nunique())
    n_hosts = int(dataset["host"].nunique())
    n_genera = int(dataset["host_genus"].nunique()) if "host_genus" in dataset.columns else 0
    if "source" in dataset.columns:
        n_within = int((dataset["source"] == "constructed_within_genus").sum())
        n_cross = int((dataset["source"] == "constructed_cross_genus").sum())
        n_neg_constructed = int(
            (dataset.get("is_constructed_negative",
                         pd.Series([False] * len(dataset))) == True).sum())
    else:
        n_within = n_cross = n_neg_constructed = 0

    # Experimental validation stats
    if exp_results is not None:
        n_exp_val = int(exp_results.get("n_pairs", 0))
        per_model = exp_results.get("per_model", {})
        if per_model:
            first = next(iter(per_model.values()))
            n_exp_pos = int(first.get("n_pos", 0))
            n_exp_neg = int(first.get("n_neg", 0))
        else:
            n_exp_pos = n_exp_neg = 0
        exp_val_source = "saureus_experimental + saureus_phage curated"
    else:
        n_exp_val = n_exp_pos = n_exp_neg = 0
        exp_val_source = "(none provided)"

    # Feature stats
    n_features = int(feature_meta.get("phage_dim", 0) + feature_meta.get("host_dim", 0)
                     + len(feature_meta.get("pair_features", [])))
    n_pair = len(feature_meta.get("pair_features", []))
    n_phages_with_seq = int(feature_meta.get("n_phages_with_seq", 0))
    n_hosts_with_seq = int(feature_meta.get("n_hosts_with_seq", 0))
    n_missing_seq = n_phages - n_phages_with_seq
    n_missing_host_seq = n_hosts - n_hosts_with_seq

    # Name-leakage results
    a_w = float(name_ablation.get("with_names", 0))
    a_s = float(name_ablation.get("sanitized_names", 0))
    a_n = float(name_ablation.get("without_names", 0))
    leakage_conclusion = name_ablation.get(
        "conclusion",
        "Inconclusive (see results/name_leakage_ablation.csv).")

    # Pairwise tests count
    n_comparisons = sum(1 for r in stat_results.get("pairwise", []))

    # Best model
    if len(summary_df):
        best_row = summary_df.sort_values("loso_mean_auc", ascending=False).iloc[0]
        best_model = str(best_row["model"])
        best_auc = float(best_row["loso_mean_auc"])
    else:
        best_model = "n/a"
        best_auc = 0.0

    # Methods skeleton
    methods = f"""## Methods

### Dataset
We assembled a phage-host interaction dataset comprising {n_pos} experimentally
verified phage-host pairs spanning {n_phages} bacteriophages and {n_hosts}
bacterial host species across {n_genera} genera. Primary data were sourced
from VirusHostInteract [CITE]; supplementary pairs were obtained from INPHARED
[CITE], Virus-Host DB [CITE], and NCBI [CITE]. The complete assembled dataset
({n_total} total pairs including {n_neg_constructed} computationally generated
negative pairs) is archived at [ZENODO DOI - INSERT BEFORE SUBMISSION]
(version {cfg["reproducibility"]["dataset_version"]}, SHA-256: {sha256}).

### Negative Pair Construction
Negative (non-infecting) phage-host pairs were computationally generated
using genus-level host-range constraints: for each phage, candidate negative
hosts were drawn from (i) the same bacterial genus as known positive hosts
(within-genus negatives, n={n_within}) and (ii) genera not associated with
any known positive interaction for that phage (cross-genus negatives,
n={n_cross}). All constructed negatives are flagged with
is_constructed_negative=True in the archived dataset. We note that absence
of a recorded interaction does not guarantee biological non-infectivity;
this is a fundamental limitation of all computational phage-host
prediction studies [CITE].

### Independent Experimental Validation
To assess performance on experimentally validated interactions not used
during model development, we evaluated all models on {n_exp_val} phage-host
pairs from {exp_val_source} ({n_exp_pos} positive, {n_exp_neg} negative).
These pairs were strictly held out from all training and hyperparameter
selection steps.

### Feature Engineering
{n_features} features were computed per phage-host pair, comprising
{feature_meta.get("phage_dim", 0)} phage-level features
(di/tri/tetra-nucleotide composition reduced via PCA, codon usage bias,
GC content, genome length, and SVD character n-gram name embeddings),
{feature_meta.get("host_dim", 0)} host-level features, and {n_pair}
pair-level interaction features (tetra-correlation, codon-usage distance,
GC-match, length ratio). For pairs where phage or host sequences were
unavailable ({n_missing_seq} phages, {n_missing_host_seq} hosts),
sequence-derived features were set to zero and flagged with
feature_was_computable=False. A complete feature manifest is provided in
results/feature_manifest.json.

### Name Embedding Leakage Audit
Phage name character n-gram embeddings (SVD, d={cfg["features"]["svd_dim"]})
were evaluated for potential label leakage arising from taxonomic naming
conventions (e.g. phages named after their host). Three conditions were
compared: original names (AUC={a_w:.4f}), names with host-genus tokens
replaced by a neutral token (AUC={a_s:.4f}), and no name features
(AUC={a_n:.4f}). {leakage_conclusion}

### Model Architecture
We trained classical baselines (RandomForest, HistGradientBoosting, and
XGBoost when available) on the concatenated per-pair feature vector, and
two graph neural network variants on a bipartite phage-host graph:
GAT (Graph Attention Network, {cfg["models"]["gat_heads"]} heads,
{cfg["models"]["gat_layers"]} layers, hidden dim {cfg["models"]["hidden_dim"]})
and GraphSAGE (mean aggregator, {cfg["models"]["gat_layers"]} layers).
Each GNN includes a parallel feed-forward bypass branch combined via a
learned sigmoid gate alpha. Per-fold alpha values are logged to
results/gnn_*_alpha_per_fold.csv. All models were trained with Adam
(lr={cfg["models"]["learning_rate"]},
weight_decay={cfg["models"]["weight_decay"]}, dropout={cfg["models"]["dropout"]})
with early stopping (patience={cfg["models"]["patience"]} epochs of no
validation-AUC improvement, max {cfg["models"]["epochs"]} epochs).

### Evaluation Protocol
Primary evaluation used Leave-One-Species-Out cross-validation (LOSO-CV)
over {len(valid_species)} host species (excluding species with fewer than
{cfg["data"]["min_positive_per_species"]} positive or
{cfg["data"]["min_negative_per_species"]} negative pairs; see Supplementary
Table species_inclusion_analysis.csv). We note that in LOSO-CV, phages
infecting multiple host species may appear in both training (via a retained
species) and the held-out test species; this evaluates host-species
generalization rather than phage generalization, which we disclose as a
limitation. Generalization across bacterial genera was assessed using
Leave-One-Genus-Out CV (LOGO-CV). All feature scaling (StandardScaler) was
fit on training folds only.

### Statistical Analysis
Model comparisons used the Wilcoxon signed-rank test on per-species LOSO-CV
AUC values (paired by species). P-values were corrected for multiple
comparisons using the Benjamini-Hochberg FDR procedure (alpha=0.05;
{n_comparisons} pairwise tests). AUC comparisons additionally used the
DeLong et al. (1988) method on pooled predictions and McNemar's test on
binary predictions. Effect sizes are reported as Cohen's d on paired
AUC differences. 95% confidence intervals used cluster bootstrap
(B={cfg["evaluation"]["bootstrap_n"]} resamples, clustered by bacterial
genus) to account for taxonomic non-independence of per-species
observations.

### Reproducibility
All analyses used random seed {cfg["reproducibility"]["seed"]}. The complete
pipeline is available at [GITHUB URL - INSERT BEFORE SUBMISSION]. The frozen
dataset (SHA-256: {sha256}) and this pipeline together enable complete
reproduction of all reported results via:

    python scripts/run_pipeline.py \\
        --frozen-dataset data/results/frozen_dataset_v{cfg["reproducibility"]["dataset_version"]}.csv \\
        --skip-download
"""
    (results_dir / "methods_skeleton.md").write_text(methods, encoding="utf-8")
    log.info(f"  wrote methods_skeleton.md")

    # Reproducibility card
    pkg_versions = run_meta.get("package_versions", {})
    card_rows = [
        ("Run timestamp", run_meta.get("timestamp", "")),
        ("Random seed", str(cfg["reproducibility"]["seed"])),
        ("Dataset version", str(cfg["reproducibility"]["dataset_version"])),
        ("Dataset SHA-256", sha256 or ""),
        ("Python", pkg_versions.get("python", sys.version.split()[0])),
        ("Platform", f"{platform.system()} {platform.release()} ({platform.machine()})"),
        ("numpy", pkg_versions.get("numpy", "")),
        ("pandas", pkg_versions.get("pandas", "")),
        ("scikit-learn", pkg_versions.get("sklearn", "")),
        ("scipy", pkg_versions.get("scipy", "")),
        ("matplotlib", pkg_versions.get("matplotlib", "")),
        ("seaborn", pkg_versions.get("seaborn", "")),
        ("xgboost", pkg_versions.get("xgboost", "")),
        ("torch", pkg_versions.get("torch", "")),
        ("torch_geometric", pkg_versions.get("torch_geometric", "")),
        ("Pairs total", str(n_total)),
        ("Pairs positive", str(n_pos)),
        ("Pairs negative", str(n_neg)),
        ("Negative within-genus", str(n_within)),
        ("Negative cross-genus", str(n_cross)),
        ("Unique phages", str(n_phages)),
        ("Unique hosts", str(n_hosts)),
        ("Genera", str(n_genera)),
        ("Phages with sequence", str(n_phages_with_seq)),
        ("Hosts with sequence", str(n_hosts_with_seq)),
        ("Features per pair", str(n_features)),
        ("Pair-level features", ", ".join(feature_meta.get("pair_features", []))),
        ("Valid LOSO species", str(len(valid_species))),
        ("Best model (LOSO)", best_model),
        ("Best LOSO mean AUC", f"{best_auc:.4f}"),
        ("Name-leakage with", f"{a_w:.4f}"),
        ("Name-leakage sanit", f"{a_s:.4f}"),
        ("Name-leakage without", f"{a_n:.4f}"),
    ]
    lines = ["# Reproducibility Card\n",
             f"_Generated: {datetime.now().isoformat()}_\n",
             "| Parameter | Value |",
             "|-----------|-------|"]
    for k, v in card_rows:
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("To reproduce all reported numbers run:")
    lines.append("")
    lines.append("```bash")
    lines.append("python scripts/run_pipeline.py \\")
    lines.append(f"    --frozen-dataset data/results/frozen_dataset_v{cfg['reproducibility']['dataset_version']}.csv \\")
    lines.append("    --skip-download")
    lines.append("```")
    (results_dir / "reproducibility_card.md").write_text("\n".join(lines), encoding="utf-8")
    log.info(f"  wrote reproducibility_card.md")
