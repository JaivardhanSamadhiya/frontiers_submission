"""Plot observed matched-cohort and nested-selection diagnostics on CPU."""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "data/results_v2/post_rejection_audit"
OUTPUT = ROOT / "output/scientific_revision_2026_10_05"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans", "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.15), constrained_layout=True)
    frame = pd.read_csv(RESULT / "matched_cohort_summary.csv")
    keys = ["shared_entities_allowed", "host_cluster_excluded", "phage_cluster_excluded", "both_clusters_excluded"]
    for i, key in enumerate(keys):
        values = frame.loc[frame.regime == key, "auroc"].to_numpy()
        axes[0].scatter(i + np.linspace(-.07, .07, len(values)), values, s=35, color="#214E74")
    axes[0].set_xticks(range(4), ["Shared\nallowed", "Host\nexcluded", "Phage\nexcluded", "Both\nexcluded"])
    axes[0].set_ylim(.7, 1)
    axes[0].set_ylabel("Pooled AUROC on the same 398 test rows")
    axes[0].set_title("A  Matched training class counts", loc="left", fontweight="bold", fontsize=10)
    rows = json.loads((RESULT / "corrected_nested_selector_summary.json").read_text())["summaries"]
    for i, row in enumerate(rows):
        delta = 100 * (row["design_lower"] - row["baseline_lower"])
        low, high = np.array(row["paired_bootstrap"]["missingness_robust_difference_interval_with_sampling"]) * 100
        axes[1].errorbar(i, delta, yerr=[[delta-low], [high-delta]], fmt="o", color="#214E74", capsize=4)
    axes[1].axhline(0, color="gray", linewidth=1, linestyle="--")
    axes[1].set_xticks(range(3), ["0.0001", "0.001", "0.01"])
    axes[1].set_xlabel("Published core-genome grouping cutoff")
    axes[1].set_ylabel("Selected minus generalist coverage\npercentage points")
    axes[1].set_title("B  Nested raw-score constituent selection", loc="left", fontweight="bold", fontsize=10)
    for ax in axes:
        ax.grid(axis="y", alpha=.18)
    fig.savefig(OUTPUT / "reassessment_figure.png", dpi=260)
    fig.savefig(OUTPUT / "reassessment_figure.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
