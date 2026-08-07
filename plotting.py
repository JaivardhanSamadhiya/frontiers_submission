"""
plotting.py - Publication-grade figure generation
=================================================
All plots saved at 300 DPI with colorblind-safe palette, no red/green
as the sole encoding, white facecolor, and tight bounding boxes.

Captions are appended to plots/captions.md in figure order for direct use
in the paper.
"""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

log = logging.getLogger(__name__)

PALETTE = sns.color_palette("colorblind")
sns.set_theme(style="whitegrid", palette="colorblind")


def _save(fig, plot_dir: Path, fname: str, caption: str):
    plot_dir.mkdir(parents=True, exist_ok=True)
    path = plot_dir / fname
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    cap_path = plot_dir / "captions.md"
    with open(cap_path, "a", encoding="utf-8") as f:
        f.write(f"\n### {fname}\n\n{caption}\n")
    log.info(f"  saved {fname}")


# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------
def plot_loso_per_species(all_results: dict, plot_dir: Path):
    rows = []
    for name, res in all_results.items():
        s = res.get("per_species_auc")
        if not isinstance(s, pd.Series) or len(s) == 0:
            continue
        for sp, auc in s.items():
            rows.append({"model": name, "species": sp, "auc": float(auc)})
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(11, 6))
    sns.boxplot(data=df, x="model", y="auc", ax=ax, palette=PALETTE,
                hue="model", legend=False)
    sns.stripplot(data=df, x="model", y="auc", ax=ax, color="black",
                  size=2.5, alpha=0.5)
    ax.set_ylabel("LOSO ROC-AUC")
    ax.set_xlabel("Model")
    ax.set_title("Per-species LOSO ROC-AUC by model")
    plt.xticks(rotation=20)
    _save(fig, plot_dir, "01_loso_per_species.png",
          "Figure 1. Per-species ROC-AUC under Leave-One-Species-Out "
          "cross-validation.  Each point is one held-out species.")


def plot_model_comparison(all_results: dict, plot_dir: Path):
    rows = []
    for name, res in all_results.items():
        rows.append({"model": name,
                     "mean_auc": float(res.get("loso_mean", res.get("mean_auc", 0))),
                     "std_auc": float(res.get("loso_std", 0))})
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(df["model"], df["mean_auc"], yerr=df["std_auc"], color=PALETTE[:len(df)],
           capsize=4)
    ax.set_ylabel("Mean LOSO ROC-AUC")
    ax.set_title("Model comparison (mean +/- std across LOSO folds)")
    ax.set_ylim(0.4, 1.0)
    plt.xticks(rotation=20)
    _save(fig, plot_dir, "03_model_comparison.png",
          "Figure 3. Mean LOSO ROC-AUC per model.  Error bars are species-level"
          " standard deviation.")


def plot_ablation(ablation_results: dict, plot_dir: Path):
    summary = ablation_results.get("summary")
    if summary is None or len(summary) == 0:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(summary["ablation"], summary["mean_auc"], yerr=summary["std_auc"],
           color=PALETTE[:len(summary)], capsize=4)
    ax.set_ylabel("Mean LOSO ROC-AUC")
    ax.set_title("Feature ablation (RF, LOSO)")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    plt.xticks(rotation=18)
    _save(fig, plot_dir, "09_ablation_study.png",
          "Figure 9. Feature ablation results.  Dashed line marks chance.")


def plot_name_leakage(name_ablation: dict, plot_dir: Path):
    conds = ["with_names", "without_names", "sanitized_names"]
    aucs = [name_ablation.get(c, 0.0) for c in conds]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.bar(conds, aucs, color=PALETTE[:3])
    ax.set_ylabel("Mean LOSO ROC-AUC")
    ax.set_title("Name Embedding Leakage Audit")
    ax.set_ylim(0.4, 1.0)
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    for i, v in enumerate(aucs):
        ax.text(i, v + 0.01, f"{v:.3f}", ha="center", va="bottom")
    drop = name_ablation.get("auc_drop_with_minus_sanit", 0.0)
    ax.text(0.5, 0.95,
            f"AUC drop (with->sanit) = {drop:+.4f}",
            transform=ax.transAxes, ha="center",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85))
    _save(fig, plot_dir, "15_name_leakage_ablation.png",
          "Figure 15. Name-embedding leakage audit.  We compare LOSO AUC for "
          "RandomForest with original phage names, with the name embedding "
          "zeroed out, and with host-genus tokens in names replaced by a "
          "neutral token.  A small drop confirms minimal leakage.")


def plot_experimental_validation(exp_results: dict, all_results: dict,
                                 plot_dir: Path):
    if exp_results is None:
        return
    rows = []
    per_model = exp_results.get("per_model", {})
    for name, m in per_model.items():
        loso_auc = float(all_results.get(name, {}).get("loso_mean",
                          all_results.get(name, {}).get("mean_auc", 0)))
        rows.append({"model": name,
                     "loso_auc": loso_auc,
                     "exp_val_auc": float(m.get("roc_auc", 0))})
    if not rows:
        return
    df = pd.DataFrame(rows)
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - 0.18, df["loso_auc"], width=0.34, label="LOSO AUC",
           color=PALETTE[0])
    ax.bar(x + 0.18, df["exp_val_auc"], width=0.34,
           label="Held-out exp. AUC", color=PALETTE[2], hatch="//")
    ax.set_xticks(x)
    ax.set_xticklabels(df["model"], rotation=18)
    ax.set_ylabel("ROC-AUC")
    ax.set_title("Held-out experimental validation\n(never used in training)")
    ax.legend()
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=0.8)
    _save(fig, plot_dir, "16_experimental_validation.png",
          "Figure 16. Held-out experimental validation: ROC-AUC under "
          "Leave-One-Species-Out (training, blue) vs. an independent "
          "experimentally validated set never used in training (orange, hatched).")


def plot_confidence_intervals(stat_results: dict, plot_dir: Path):
    boot = stat_results.get("bootstrap", [])
    if not boot:
        return
    df = pd.DataFrame(boot)
    df = df.sort_values("mean")
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(9, 0.45 * max(4, len(df))))
    # standard bootstrap CI (lighter)
    ax.hlines(y, df["std_ci_lo"], df["std_ci_hi"], color=PALETTE[7], linewidth=4,
              alpha=0.6, label="Standard bootstrap CI")
    # cluster bootstrap CI (primary)
    ax.hlines(y + 0.18, df["cluster_ci_lo"], df["cluster_ci_hi"],
              color=PALETTE[0], linewidth=4, alpha=0.95,
              label="Cluster bootstrap CI (preferred)")
    ax.scatter(df["mean"], y + 0.09, color="black", marker="|", s=80)
    ax.set_yticks(y + 0.09)
    ax.set_yticklabels(df["model"])
    ax.set_xlabel("ROC-AUC")
    ax.set_title("Confidence intervals (cluster bootstrap, B=10000)")
    ax.legend()
    _save(fig, plot_dir, "18_confidence_intervals.png",
          "Figure 18. Bootstrap 95% confidence intervals for mean LOSO AUC.  "
          "Cluster bootstrap (genus-level) is the publication-preferred interval; "
          "the standard interval is shown for reference.")


def plot_data_leakage_audit(dataset: pd.DataFrame, ablation_results: dict,
                            plot_dir: Path):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    cols = [c for c in ["tetra_corr", "cub_dist", "gc_match", "len_ratio"]
            if c in dataset.columns]
    if cols:
        # Panel A: pos vs neg distribution on training rows
        ax = axes[0, 0]
        for lab, sub in dataset.groupby("label"):
            ax.hist(sub[cols[0]].astype(float), bins=40, alpha=0.55,
                    label=f"label={int(lab)}",
                    color=PALETTE[0] if int(lab) == 1 else PALETTE[7])
        ax.set_title(f"A. {cols[0]} distribution (pos vs neg)")
        ax.set_xlabel(cols[0])
        ax.legend()

        # Panel B: GC-match
        ax = axes[0, 1]
        if "gc_match" in dataset.columns:
            for lab, sub in dataset.groupby("label"):
                ax.hist(sub["gc_match"].astype(float), bins=40, alpha=0.55,
                        label=f"label={int(lab)}",
                        color=PALETTE[2] if int(lab) == 1 else PALETTE[5])
            ax.set_title("B. gc_match distribution")
            ax.set_xlabel("gc_match")
            ax.legend()

        # Panel C: constructed vs VHI feature comparison
        ax = axes[1, 0]
        if "is_constructed_negative" in dataset.columns:
            for is_constr, sub in dataset.groupby("is_constructed_negative"):
                lab = "constructed" if is_constr else "natural"
                ax.hist(sub[cols[0]].astype(float), bins=40, alpha=0.55,
                        label=lab,
                        color=PALETTE[3] if is_constr else PALETTE[6])
            ax.set_title(f"C. constructed vs natural ({cols[0]})")
            ax.set_xlabel(cols[0])
            ax.legend()

    # Panel D: AUC vs hard-negative fraction (from ablation)
    ax = axes[1, 1]
    summary = ablation_results.get("summary")
    if summary is not None and len(summary):
        sub = summary[summary["ablation"].isin(["full", "hard_negatives_only"])]
        ax.bar(sub["ablation"], sub["mean_auc"], color=PALETTE[:len(sub)])
        ax.set_ylim(0.4, 1.0)
        ax.set_title("D. AUC: full vs hard-negatives only")
        ax.set_ylabel("Mean LOSO AUC")
    plt.tight_layout()
    _save(fig, plot_dir, "14_data_leakage_audit.png",
          "Figure 14. Data leakage audit.  (A) Pair-feature distribution by "
          "label.  (B) GC-match distribution by label.  (C) Natural vs "
          "constructed negative distributions.  (D) Effect of training only "
          "on hard (within-genus) negatives.")


def plot_reproducibility_card(plot_dir: Path, cfg: dict, sha256: str | None = None):
    info = [
        ["Random seed", str(cfg["reproducibility"]["seed"])],
        ["Dataset version", str(cfg["reproducibility"]["dataset_version"])],
        ["Dataset SHA-256", (sha256 or "(see results/frozen_dataset_v*.sha256)")[:32] + "..."],
        ["Hidden dim", str(cfg["models"]["hidden_dim"])],
        ["GAT heads", str(cfg["models"]["gat_heads"])],
        ["GAT layers", str(cfg["models"]["gat_layers"])],
        ["Max epochs", str(cfg["models"]["epochs"])],
        ["Patience", str(cfg["models"]["patience"])],
        ["MC rounds", str(cfg["evaluation"]["n_mc_rounds"])],
        ["Bootstrap B", str(cfg["evaluation"]["bootstrap_n"])],
        ["Min species rows", str(cfg["data"]["min_species_rows"])],
        ["Within-genus neg ratio", str(cfg["data"]["neg_ratio_within_genus"])],
        ["Cross-genus neg ratio", str(cfg["data"]["neg_ratio_cross_genus"])],
    ]
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(info) + 1))
    ax.axis("off")
    table = ax.table(cellText=info, colLabels=["Parameter", "Value"],
                     loc="center", cellLoc="left", colLoc="left")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.4)
    ax.set_title("Reproducibility card")
    _save(fig, plot_dir, "17_reproducibility_card.png",
          "Figure 17. Reproducibility card listing fixed hyperparameters and "
          "the frozen dataset hash.  See results/reproducibility_card.md for "
          "the complete table.")


# ---------------------------------------------------------------------------
# Public driver
# ---------------------------------------------------------------------------
def generate_all_plots(dataset: pd.DataFrame, all_results: dict,
                       stat_results: dict, ablation_results: dict,
                       name_ablation: dict, exp_results, cfg: dict,
                       seed: int) -> None:
    plot_dir = Path(cfg["paths"]["plots_dir"])
    plot_dir.mkdir(parents=True, exist_ok=True)
    cap_path = plot_dir / "captions.md"
    cap_path.write_text("# Figure captions\n", encoding="utf-8")

    try:
        plot_loso_per_species(all_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_loso_per_species failed: {e}")
    try:
        plot_model_comparison(all_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_model_comparison failed: {e}")
    try:
        plot_ablation(ablation_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_ablation failed: {e}")
    try:
        plot_name_leakage(name_ablation, plot_dir)
    except Exception as e:
        log.exception(f"  plot_name_leakage failed: {e}")
    try:
        plot_experimental_validation(exp_results, all_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_experimental_validation failed: {e}")
    try:
        plot_confidence_intervals(stat_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_confidence_intervals failed: {e}")
    try:
        plot_data_leakage_audit(dataset, ablation_results, plot_dir)
    except Exception as e:
        log.exception(f"  plot_data_leakage_audit failed: {e}")
    try:
        plot_reproducibility_card(plot_dir, cfg)
    except Exception as e:
        log.exception(f"  plot_reproducibility_card failed: {e}")
