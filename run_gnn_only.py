#!/usr/bin/env python3
"""
run_gnn_only.py - GNN evaluation on the frozen genomically-grounded subset
==========================================================================
Reuses the frozen dataset produced by the main pipeline (same subset, same
seed), recomputes features deterministically, then runs the GAT and SAGE
pipelines so GNN results can be compared against the classical backbone
WITHOUT re-running the slow classical + ablation stages.

Usage:
  env -u PYTHONHOME -u PYTHONPATH .venv/bin/python run_gnn_only.py \
      --config scripts/config.yaml \
      --frozen-dataset data/results/frozen_dataset_v1.0.0.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "scripts" / "config.yaml"))
    ap.add_argument("--frozen-dataset",
                    default=str(ROOT / "data" / "results" / "frozen_dataset_v1.0.0.csv"))
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    log = logging.getLogger("gnn_only")

    import run_pipeline as rp
    cfg = rp.resolve_paths(yaml.safe_load(open(args.config)))

    import gnn as gnn_mod
    log.info("HAS_TORCH=%s", gnn_mod.HAS_TORCH)

    dataset = pd.read_csv(args.frozen_dataset)
    log.info("Loaded frozen dataset: %s rows (pos=%d neg=%d)",
             len(dataset), int((dataset.label == 1).sum()),
             int((dataset.label == 0).sum()))

    from data_enrichment import compute_all_features
    from model import build_feature_matrices
    from evaluation import get_valid_species, get_valid_genera

    seed = int(cfg["reproducibility"]["seed"])
    dataset, phage_feat_df, host_feat_df, feature_meta = compute_all_features(
        dataset, cfg, seed=seed)
    matrices = build_feature_matrices(dataset, phage_feat_df, host_feat_df,
                                      feature_meta, cfg, seed=seed)
    valid_species = get_valid_species(dataset, cfg)
    valid_genera = get_valid_genera(dataset, cfg)
    log.info("valid LOSO species=%d  valid LOGO genera=%d",
             len(valid_species), len(valid_genera))

    from gnn import run_gnn_pipeline
    summary = {}
    for arch in ["GAT", "SAGE"]:
        log.info("=== Running %s ===", arch)
        res = run_gnn_pipeline(arch, dataset, matrices, valid_species,
                               valid_genera, cfg, seed=seed)
        alphas = [a for a in res.get("alpha_values", []) if a == a]  # drop nan
        summary[arch] = {
            "loso_mean": res["loso_mean"],
            "loso_std": res["loso_std"],
            "loso_pooled": res["loso_pooled"],
            "logo_mean": res.get("logo_mean", 0.0),
            "mc_unseen_auc": res.get("mc_auc", 0.0),
            "used_fallback": res.get("used_fallback"),
            "n_fallback_folds": res.get("n_fallback_folds"),
            "arch_label": res.get("arch_label"),
            "mean_alpha": (sum(alphas) / len(alphas)) if alphas else None,
            "n_alpha": len(alphas),
        }
        log.info("%s: LOSO mean=%.4f std=%.4f pooled=%.4f | LOGO=%.4f | "
                 "unseen=%.4f | fallback=%s | mean_alpha=%s",
                 arch, summary[arch]["loso_mean"], summary[arch]["loso_std"],
                 summary[arch]["loso_pooled"], summary[arch]["logo_mean"],
                 summary[arch]["mc_unseen_auc"], summary[arch]["used_fallback"],
                 summary[arch]["mean_alpha"])

    out = Path(cfg["paths"]["results_dir"]) / "gnn_only_summary.json"
    out.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    log.info("Wrote %s", out)


if __name__ == "__main__":
    main()
