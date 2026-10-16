"""Retrospective reporting controls for historical ColoColi cocktail scores.

No model is fitted and no source file is changed. The outcome is the recorded
positive assay of the one historical selected triple per host. Methods added
here are retrospective controls. Their scores are not calibrated probabilities.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.metrics import (
    average_precision_score, precision_recall_curve, roc_auc_score, roc_curve,
)

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "data/results_v2/external_cocktail_validation"
UPSTREAM = ROOT / "external/gaborieau_coli_phage_interactions_2023"
OUT = ROOT / "data/results_v2/frontiers_revision_2026_10_14/cocktail_reporting"
SEED = 42030  # Historical experiment 30 bootstrap seed.
REPLICATES = 2000
METHODS = {
    "historical_product": "Historical product score",
    "maximum_constituent": "Maximum constituent",
    "mean_constituent": "Mean constituent",
    "logistic_component": "Logistic component",
    "neighbour_component": "Neighbour component",
    "training_phage_mean": "Training phage susceptibility",
    "host_phylogroup": "Host phylogroup only",
}
INPUTS = {
    "cocktails": SOURCE / "heldout_cocktail_predictions_and_outcomes.csv",
    "constituents": SOURCE / "heldout_selected_pair_predictions_and_outcomes.csv",
    "all_frozen_pairs": SOURCE / "frozen_pair_predictions.csv",
    "historical_summary": SOURCE / "external_cocktail_validation_summary.json",
    "frozen_manifest": SOURCE / "frozen_prediction_manifest.json",
    "training_matrix": UPSTREAM / "data/interactions/interaction_matrix.csv",
    "training_metadata": UPSTREAM / "data/genomics/bacteria/picard_collection.csv",
    "external_metadata": UPSTREAM / "dev/cocktails/data/test_collection.csv",
    "distances": UPSTREAM / "dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv",
    "host_cluster_definition": ROOT / "experiments/30_evaluate_external_cocktails.py",
    "composition_shuffle": ROOT / "data/results_v2/post_rejection_audit/composition_shuffle_diagnostic.csv",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def dump_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def historical_clusters(names: list[str], distance: pd.DataFrame) -> np.ndarray:
    spec = importlib.util.spec_from_file_location("historical_cocktail_evaluation", INPUTS["host_cluster_definition"])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.host_clusters(names, distance)


def weighted_metrics_many(y: np.ndarray, score: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Weighted AUROC/AP/prevalence/AP-minus-prevalence, preserving score ties.

    Rows of weights represent bootstrap draws; columns represent original hosts.
    Grouping tied scores reproduces sklearn's AP and half-credit tied AUROC.
    """
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    weights = np.atleast_2d(np.asarray(weights, dtype=float))
    if len(y) != len(score) or weights.shape[1] != len(y):
        raise ValueError("Label/score/weight dimensions differ")
    if not np.isin(y, [0, 1]).all() or not np.isfinite(score).all():
        raise ValueError("Invalid labels or scores")
    if not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError("Invalid bootstrap weights")
    order = np.argsort(-score, kind="stable")
    sorted_score = score[order]
    starts = np.r_[0, np.flatnonzero(sorted_score[1:] != sorted_score[:-1]) + 1]
    positive = np.add.reduceat(weights[:, order] * y[order], starts, axis=1)
    negative = np.add.reduceat(weights[:, order] * (1 - y[order]), starts, axis=1)
    total_positive = positive.sum(axis=1)
    total_negative = negative.sum(axis=1)
    if (total_positive <= 0).any() or (total_negative <= 0).any():
        raise ValueError("Both outcome classes are required")
    cumulative_positive = np.cumsum(positive, axis=1)
    cumulative_negative = np.cumsum(negative, axis=1)
    auroc = np.sum(positive * (
        total_negative[:, None] - cumulative_negative + 0.5 * negative
    ), axis=1) / (total_positive * total_negative)
    precision = np.divide(cumulative_positive, cumulative_positive + cumulative_negative,
                          out=np.zeros_like(cumulative_positive),
                          where=(cumulative_positive + cumulative_negative) > 0)
    ap = np.sum(positive * precision, axis=1) / total_positive
    prevalence = total_positive / (total_positive + total_negative)
    return np.column_stack([auroc, ap, prevalence, ap - prevalence])


def bootstrap_weights(y: np.ndarray, clusters: np.ndarray) -> tuple[np.ndarray, int]:
    unique, inverse = np.unique(clusters, return_inverse=True)
    rng = np.random.default_rng(SEED)
    draws = []
    attempts = 0
    while len(draws) < REPLICATES and attempts < 50 * REPLICATES:
        attempts += 1
        sample = rng.choice(unique, len(unique), replace=True)
        counts = np.array([np.sum(sample == group) for group in unique])
        weight = counts[inverse]
        if np.unique(y[weight > 0]).size == 2:
            draws.append(weight)
    if len(draws) != REPLICATES:
        raise RuntimeError("Insufficient valid two-class bootstrap draws")
    return np.asarray(draws, dtype=float), attempts


def score_frame() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    cocktails = pd.read_csv(INPUTS["cocktails"])
    pairs = pd.read_csv(INPUTS["constituents"])
    all_pairs = pd.read_csv(INPUTS["all_frozen_pairs"])
    manifest = json.loads(INPUTS["frozen_manifest"].read_text(encoding="utf-8"))
    summary = json.loads(INPUTS["historical_summary"].read_text(encoding="utf-8"))
    if sha256(INPUTS["all_frozen_pairs"]) != manifest["checksums"]["frozen_pair_predictions"]:
        raise RuntimeError("Frozen prediction checksum mismatch")
    if manifest["checksums"]["frozen_pair_predictions"] != summary["frozen_prediction_sha256"]:
        raise RuntimeError("Historical manifests disagree")
    assert len(cocktails) == 100 and cocktails.bacteria.is_unique
    assert len(pairs) == 300 and pairs.groupby("bacteria").size().eq(3).all()
    assert not pairs.duplicated(["bacteria", "phage"]).any()
    assert len(all_pairs) == 9600 and all_pairs.groupby("bacteria").size().eq(96).all()
    assert int(cocktails.productive.sum()) == 82 and int(pairs.productive.sum()) == 170
    assert set(cocktails.bacteria) == set(pairs.bacteria) == set(all_pairs.bacteria)
    frame = cocktails[["bacteria", "cocktail_mlc", "productive", "predicted_constituents"]].copy()
    frame["historical_product"] = cocktails.score_independence
    frame["maximum_constituent"] = cocktails.score_max
    pair_groups = dict(tuple(pairs.groupby("bacteria", sort=False)))
    for host, row in frame.set_index("bacteria").iterrows():
        group = pair_groups[host]
        assert set(group.phage) == set(row.predicted_constituents.split("|"))
        for field in ["score_ensemble", "score_logistic", "score_knn", "training_prevalence"]:
            assert np.isfinite(group[field]).all() and group[field].between(0, 1).all()
    aggregate = pairs.groupby("bacteria", sort=False).agg(
        mean_constituent=("score_ensemble", lambda values: float(np.mean(values.to_numpy(float)))),
        training_phage_mean=("training_prevalence", lambda values: float(np.mean(values.to_numpy(float)))),
    )
    for input_key, output_key in [("score_logistic", "logistic_component"),
                                   ("score_knn", "neighbour_component")]:
        aggregate[output_key] = pairs.groupby("bacteria", sort=False)[input_key].apply(
            lambda values: 1.0 - np.prod(1.0 - values.to_numpy(float)))
    recomputed_product = pairs.groupby("bacteria", sort=False).score_ensemble.apply(
        lambda values: 1.0 - np.prod(1.0 - values.to_numpy(float)))
    recomputed_max = pairs.groupby("bacteria", sort=False).score_ensemble.max()
    np.testing.assert_allclose(frame.historical_product, recomputed_product.loc[frame.bacteria], atol=1e-14, rtol=0)
    np.testing.assert_allclose(frame.maximum_constituent, recomputed_max.loc[frame.bacteria], atol=1e-14, rtol=0)
    frame = frame.join(aggregate, on="bacteria", validate="one_to_one")

    matrix = pd.read_csv(INPUTS["training_matrix"], sep=";", index_col=0)
    training_meta = pd.read_csv(INPUTS["training_metadata"], sep=";").set_index("bacteria")
    external_meta = pd.read_csv(INPUTS["external_metadata"], sep=";").set_index("bacteria")
    assert matrix.shape == (402, 96) and matrix.index.is_unique
    assert training_meta.index.is_unique and external_meta.index.is_unique
    training_meta = training_meta.loc[matrix.index]
    external_meta = external_meta.loc[frame.bacteria]
    observed = matrix.notna()
    positive = matrix.gt(0) & observed
    global_mean = float(positive.to_numpy().sum() / observed.to_numpy().sum())
    assert int(observed.to_numpy().sum()) == 38435
    np.testing.assert_allclose(global_mean, manifest["training_positive_fraction"], atol=1e-15, rtol=0)
    training_prevalence = positive.sum(axis=0) / observed.sum(axis=0)
    np.testing.assert_allclose(pairs.training_prevalence, training_prevalence.loc[pairs.phage], atol=1e-15, rtol=0)
    phylo_counts = []
    phylo_score = {}
    for phylo, part in training_meta.groupby("Clermont_Phylo", dropna=True, sort=True):
        names = part.index
        n_observed = int(observed.loc[names].to_numpy().sum())
        n_positive = int(positive.loc[names].to_numpy().sum())
        score = n_positive / n_observed if n_observed else global_mean
        phylo_score[str(phylo)] = score
        phylo_counts.append({"phylogroup": str(phylo), "training_hosts": len(names),
                             "observed_individual_assays": n_observed,
                             "positive_individual_assays": n_positive,
                             "score": score})
    frame["phylogroup"] = external_meta.Clermont_Phylo.astype("string").fillna("Unknown").to_numpy()
    frame["phylogroup_fallback"] = ~frame.phylogroup.isin(phylo_score)
    frame["host_phylogroup"] = frame.phylogroup.map(phylo_score).fillna(global_mean)
    distance = pd.read_csv(INPUTS["distances"], sep="\t", index_col=0)
    frame["host_cluster"] = historical_clusters(frame.bacteria.tolist(), distance)
    frame["missing_distance_row"] = ~frame.bacteria.isin(distance.index)
    assert frame.host_cluster.nunique() == 81
    assert frame.missing_distance_row.sum() == 3
    assert np.isfinite(frame[list(METHODS)].to_numpy()).all()
    return frame, pd.DataFrame(phylo_counts), {
        "training_hosts": 402, "bank_phages": 96,
        "observed_training_individual_assays": int(observed.to_numpy().sum()),
        "positive_training_individual_assays": int(positive.to_numpy().sum()),
        "global_training_individual_activity": global_mean,
        "phylogroup_fallback_hosts": frame.loc[frame.phylogroup_fallback, "bacteria"].tolist(),
        "missing_distance_hosts_given_separate_clusters": frame.loc[frame.missing_distance_row, "bacteria"].tolist(),
        "phylogroup_baseline_rule": "Mean binary individual-assay activity over observed Picard host-phage pairs within the host phylogroup; missing/unseen phylogroups use global observed-training mean. No cocktail labels, tuning, or model fitting.",
    }


def metrics_for_cohort(frame: pd.DataFrame, cohort: str) -> tuple[list[dict], list[dict], np.ndarray, dict]:
    y = frame.productive.to_numpy(int)
    weights, attempts = bootstrap_weights(y, frame.host_cluster.to_numpy())
    metric_names = ["auroc", "average_precision", "positive_prevalence", "ap_minus_prevalence"]
    rows, boot, point = [], {}, {}
    for method in METHODS:
        score = frame[method].to_numpy(float)
        point[method] = weighted_metrics_many(y, score, np.ones((1, len(y))))[0]
        np.testing.assert_allclose(point[method][:2], [roc_auc_score(y, score), average_precision_score(y, score)], atol=1e-13, rtol=0)
        boot[method] = weighted_metrics_many(y, score, weights)
        lo, hi = np.quantile(boot[method], [.025, .975], axis=0)
        row = {"cohort": cohort, "method": method, "label": METHODS[method],
               "n_hosts": len(frame), "positive": int(y.sum()),
               "negative": int(len(y) - y.sum()), "host_clusters": int(frame.host_cluster.nunique())}
        for j, name in enumerate(metric_names):
            row[name] = float(point[method][j])
            row[name + "_ci_low"] = float(lo[j])
            row[name + "_ci_high"] = float(hi[j])
        rows.append(row)
    delta_rows = []
    for first, second in combinations(METHODS, 2):
        difference = boot[first] - boot[second]
        lo, hi = np.quantile(difference, [.025, .975], axis=0)
        row = {"cohort": cohort, "method_a": first, "method_b": second,
               "difference_definition": "method_a minus method_b"}
        for j in [0, 1, 3]:
            name = metric_names[j]
            row[name + "_difference"] = float(point[first][j] - point[second][j])
            row[name + "_difference_ci_low"] = float(lo[j])
            row[name + "_difference_ci_high"] = float(hi[j])
        delta_rows.append(row)
    expected_random_ap = ((int(y.sum()) - 1) / (len(y) - 1)
                          + (len(y) - int(y.sum())) * sum(1.0 / k for k in range(1, len(y) + 1))
                          / (len(y) * (len(y) - 1)))
    return rows, delta_rows, weights, {
        "cohort": cohort, "n_hosts": len(frame), "positive": int(y.sum()),
        "host_clusters": int(frame.host_cluster.nunique()),
        "valid_replicates": len(weights), "attempted_replicates": attempts,
        "seed": SEED, "confidence_interval": "2.5th and 97.5th percentiles",
        "constant_score_ap": float(y.mean()),
        "expected_uniform_random_order_ap": expected_random_ap,
        "expected_random_ap_formula": "(P-1)/(N-1) + (N-P)*H_N/(N*(N-1)) for fixed P positives among N records",
        "resampling_scope": "Paired resampling of the historical external host clusters at distance 1e-4; all methods use identical multiplicity weights; row-weighted pooled metrics; conditional on the saved model, historical cohort normalization, chosen grouping and fixed phage bank.",
        "ap_excess_covariance": "Prevalence is recomputed using the same weights within every replicate before subtracting it from AP.",
    }


def triage_for_cohort(frame: pd.DataFrame, cohort: str) -> tuple[list[dict], list[dict], list[dict]]:
    topk_rows, recovery_rows, curve_rows = [], [], []
    n = len(frame)
    total_positive = int(frame.productive.sum())
    for method in METHODS:
        ranked = frame.sort_values([method, "bacteria"], ascending=[False, True], kind="stable")
        cumulative = ranked.productive.to_numpy(int).cumsum()
        for k in [10, 25, 50]:
            hits = int(cumulative[k - 1])
            topk_rows.append({
                "cohort": cohort, "method": method, "k": k, "n_hosts": n,
                "positive_total": total_positive, "observed_positive_recovered": hits,
                "precision_at_k": hits / k, "recall_at_k": hits / total_positive,
                "random_expected_positives_at_k": k * total_positive / n,
                "random_expected_precision_at_k": total_positive / n,
                "random_expected_recall_at_k": k / n,
            })
        for target in [.5, .8, .9]:
            required = int(math.ceil(target * total_positive))
            rank = int(np.searchsorted(cumulative, required) + 1)
            random_rank = required * (n + 1) / (total_positive + 1)
            recovery_rows.append({
                "cohort": cohort, "method": method, "target_recall": target,
                "positive_records_required": required, "actual_recall_at_requirement": required / total_positive,
                "observed_rank_to_recover": rank,
                "expected_random_rank_to_recover": random_rank,
                "random_expected_rank_minus_observed": random_rank - rank,
            })
        for index, (host, observed) in enumerate(zip(ranked.bacteria, ranked.productive), start=1):
            curve_rows.append({
                "cohort": cohort, "method": method, "rank": index,
                "bacteria": host, "observed_positive": int(observed),
                "cumulative_positive": int(cumulative[index - 1]),
                "precision": cumulative[index - 1] / index,
                "recall": cumulative[index - 1] / total_positive,
                "random_expected_recall": index / n,
            })
    return topk_rows, recovery_rows, curve_rows


def self_checks() -> dict:
    rng = np.random.default_rng(20261014)
    cases = 0
    for n in [4, 10, 31, 100]:
        for trial in range(12):
            y = rng.integers(0, 2, n)
            y[0], y[-1] = 0, 1
            score = rng.integers(0, 4, n).astype(float) if trial % 2 else rng.random(n)
            weight = rng.integers(0, 4, (5, n)).astype(float)
            weight[:, 0] = weight[:, -1] = 1
            actual = weighted_metrics_many(y, score, weight)
            for i, w in enumerate(weight):
                expected = [roc_auc_score(y, score, sample_weight=w),
                            average_precision_score(y, score, sample_weight=w),
                            np.average(y, weights=w)]
                np.testing.assert_allclose(actual[i, :3], expected, atol=1e-13, rtol=0)
                np.testing.assert_allclose(actual[i, 3], expected[1] - expected[2], atol=1e-13, rtol=0)
                cases += 1
    constant = weighted_metrics_many(np.array([0, 1, 1, 1]), np.ones(4), np.ones((1, 4)))[0]
    np.testing.assert_allclose(constant, [.5, .75, .75, 0], atol=1e-15, rtol=0)
    from itertools import combinations as choose
    random_recovery_cases = 0
    for n in range(3, 8):
        for p in range(1, n):
            locations = np.asarray(list(choose(range(1, n + 1), p)))
            for r in range(1, p + 1):
                np.testing.assert_allclose(locations[:, r - 1].mean(), r * (n + 1) / (p + 1), atol=1e-14, rtol=0)
                random_recovery_cases += 1
    return {"sklearn_weighted_metric_cases": cases,
            "constant_score_ap_and_zero_lift": True,
            "exhaustive_expected_random_recovery_cases": random_recovery_cases,
            "passed": True}


def figures(frame: pd.DataFrame, summary: pd.DataFrame, curves: pd.DataFrame) -> list[str]:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12,
                         "axes.labelsize": 13, "axes.titlesize": 14,
                         "legend.fontsize": 10.5, "xtick.labelsize": 11,
                         "ytick.labelsize": 11, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    chosen = ["historical_product", "maximum_constituent", "mean_constituent",
              "training_phage_mean", "host_phylogroup"]
    colors = {"historical_product": "#145DA0", "maximum_constituent": "#D88022",
              "mean_constituent": "#20886E", "training_phage_mean": "#8C519B",
              "host_phylogroup": "#565656"}
    lines = {"historical_product": "-", "maximum_constituent": "--",
             "mean_constituent": "-.", "training_phage_mean": ":", "host_phylogroup": "--"}
    y = frame.productive.to_numpy(int)
    main_summary = summary[summary.cohort.eq("historical_100")].set_index("method")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.7))
    for method in chosen:
        score = frame[method].to_numpy(float)
        fpr, tpr, _ = roc_curve(y, score)
        precision, recall, _ = precision_recall_curve(y, score)
        axes[0].plot(fpr, tpr, color=colors[method], ls=lines[method], lw=2,
                     label=f"{METHODS[method]} ({main_summary.loc[method, 'auroc']:.3f})")
        axes[1].step(recall, precision, where="post", color=colors[method], ls=lines[method], lw=2,
                     label=f"{METHODS[method]} ({main_summary.loc[method, 'average_precision']:.3f})")
    axes[0].plot([0, 1], [0, 1], color="0.65", ls=":", lw=1.5)
    axes[1].axhline(y.mean(), color="0.45", ls=":", lw=1.5,
                    label=f"Positive prevalence ({y.mean():.3f})")
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="A  ROC (legend: AUROC)", xlim=(0, 1), ylim=(0, 1.02))
    axes[1].set(xlabel="Recall", ylabel="Precision", title="B  Precision–recall (legend: AP)", xlim=(0, 1), ylim=(0, 1.02))
    for ax in axes:
        ax.grid(alpha=.15)
        ax.legend(loc="lower right", frameon=False)
    fig.suptitle("Historical selected triples: 100 hosts, 82 positive assays", fontsize=15)
    fig.subplots_adjust(left=.07, right=.99, bottom=.14, top=.85, wspace=.26)
    paths = []
    for extension in ["png", "pdf"]:
        path = OUT / f"cocktail_roc_pr.{extension}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path.name)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.7))
    for method in chosen:
        part = curves[curves.cohort.eq("historical_100") & curves.method.eq(method)]
        axes[0].step(np.r_[0, part["rank"]], np.r_[0, part.recall], where="post",
                     color=colors[method], ls=lines[method], lw=2, label=METHODS[method])
        axes[1].plot(part["rank"], part.precision, color=colors[method], ls=lines[method], lw=2)
    axes[0].plot([0, len(frame)], [0, 1], color="0.4", ls=":", lw=1.5, label="Random-order expectation")
    axes[1].axhline(y.mean(), color="0.4", ls=":", lw=1.5)
    axes[0].set(xlabel="Host–cocktail records screened", ylabel="Fraction of positive records recovered",
                title="A  Retrospective recovery", xlim=(0, 100), ylim=(0, 1.02))
    axes[1].set(xlabel="Host–cocktail records screened", ylabel="Precision among screened records",
                title="B  Precision by screening budget", xlim=(0, 100), ylim=(0, 1.02))
    axes[0].legend(loc="lower right", frameon=False)
    for ax in axes:
        ax.grid(alpha=.15)
    fig.suptitle("Cohort triage of previously tested triples", fontsize=15)
    fig.subplots_adjust(left=.07, right=.99, bottom=.14, top=.85, wspace=.26)
    for extension in ["png", "pdf"]:
        path = OUT / f"cocktail_triage.{extension}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path.name)
    plt.close(fig)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    checks = self_checks()
    if args.self_test_only:
        print(json.dumps(checks, indent=2))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    hashes_before = {name: sha256(path) for name, path in INPUTS.items()}
    frame, phylo_reference, baseline = score_frame()
    all_metric_rows, all_difference_rows, all_topk, all_recovery, all_curves, cohorts = [], [], [], [], [], []
    for cohort, part in [("historical_100", frame),
                         ("exclude_AN03_99", frame.loc[frame.bacteria.ne("AN03")].copy())]:
        assert len(part) == (100 if cohort == "historical_100" else 99)
        metric_rows, difference_rows, weights, meta = metrics_for_cohort(part, cohort)
        all_metric_rows.extend(metric_rows)
        all_difference_rows.extend(difference_rows)
        topk, recovery, curves = triage_for_cohort(part, cohort)
        all_topk.extend(topk)
        all_recovery.extend(recovery)
        all_curves.extend(curves)
        cohorts.append(meta)
        np.savez_compressed(OUT / f"bootstrap_{cohort}.npz", weights=weights,
                            bacteria=part.bacteria.to_numpy(str),
                            clusters=part.host_cluster.to_numpy(int))
    metrics = pd.DataFrame(all_metric_rows)
    differences = pd.DataFrame(all_difference_rows)
    topk = pd.DataFrame(all_topk)
    recovery = pd.DataFrame(all_recovery)
    curves = pd.DataFrame(all_curves)
    historical = metrics[(metrics.cohort == "historical_100") & (metrics.method == "historical_product")].iloc[0]
    original = json.loads(INPUTS["historical_summary"].read_text())["primary_cocktail"]
    np.testing.assert_allclose(historical.auroc, original["score_independence"]["roc_auc"], atol=1e-14, rtol=0)
    np.testing.assert_allclose(historical.average_precision, original["score_independence"]["average_precision"], atol=1e-14, rtol=0)
    np.testing.assert_allclose([historical.auroc_ci_low, historical.auroc_ci_high], original["cluster_bootstrap"]["roc_auc_95_ci"], atol=1e-13, rtol=0)
    np.testing.assert_allclose([historical.average_precision_ci_low, historical.average_precision_ci_high], original["cluster_bootstrap"]["average_precision_95_ci"], atol=1e-13, rtol=0)
    shuffle = pd.read_csv(INPUTS["composition_shuffle"])
    assert len(shuffle) == 1000 and set(shuffle) == {"auroc", "average_precision"}
    shuffle_summary = {key: {"median": float(shuffle[key].median()),
                              "diagnostic_95_range": shuffle[key].quantile([.025, .975]).tolist()}
                       for key in shuffle}
    tables = {"host_scores.csv": frame, "training_phylogroup_reference.csv": phylo_reference,
              "method_metrics.csv": metrics, "paired_method_differences.csv": differences,
              "triage_top_k.csv": topk, "triage_recovery.csv": recovery,
              "triage_ranked_records.csv": curves}
    for name, table in tables.items():
        table.to_csv(OUT / name, index=False)
    figure_paths = figures(frame, metrics, curves)
    hashes_after = {name: sha256(path) for name, path in INPUTS.items()}
    if hashes_before != hashes_after:
        raise RuntimeError("Source files changed during the reporting analysis")
    checks.update({"historical_point_metrics_reproduced": True,
                   "historical_2000_draw_confidence_intervals_reproduced": True,
                   "frozen_prediction_hash_verified": True,
                   "training_prevalences_recomputed": True,
                   "all_source_hashes_unchanged": True,
                   "method_rows": len(metrics), "paired_difference_rows": len(differences)})
    limitations = [
        "These are retrospective controls on already inspected external outcomes, not a new frozen validation or preregistration.",
        "Historical product, maximum and mean scores inherit normalization across the 100-host external cohort; subgroup sensitivity preserves those original scores rather than renormalizing.",
        "The outcome is the recorded positive assay of each externally selected historical triple. Naming the source column productive in a legacy CSV does not establish productive replication or clinical efficacy.",
        "One analyzed triple per host ranks records across the cohort; it does not evaluate alternative cocktails for the same host. Source measured pairs/control triples are handled in the separate source-comparison analysis.",
        "All uncertainty conditions on fixed scores, bank and historical host grouping; three missing-distance hosts receive singleton groups and the historical AN03 name is ambiguous.",
        "Seven methods and multiple reporting contrasts are exploratory; no multiplicity-adjusted superiority hypothesis or permutation p-value is asserted.",
        "Composition shuffle is reused as a diagnostic distribution; exchangeability sufficient for a confirmatory null test has not been established.",
        "Triage ranks use descending score and ascending exact host identifier to break ties. Results concern the already selected tested cohort, not prospective thresholds, laboratory cost savings or within-host selection.",
        "Random recovery expectation r*(N+1)/(P+1) is the expected rank of the rth positive record under a uniform permutation of this fixed labelled cohort; it is not a real-world cost estimate.",
    ]
    report = {"analysis": "Saved-score historical cocktail reporting and controls",
              "status": "completed retrospective descriptive analysis",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "code_sha256": sha256(Path(__file__)), "seed": SEED,
              "bootstrap_replicates_per_cohort": REPLICATES,
              "cohorts": cohorts, "methods": METHODS, "training_baselines": baseline,
              "historical_composition_shuffle": {"replicates": len(shuffle), "scope": "Historical 100-host cohort only; diagnostic, no p-value", **shuffle_summary},
              "checks": checks, "limitations": limitations,
              "input_hashes": {name: {"path": str(INPUTS[name].relative_to(ROOT)), "sha256": value}
                               for name, value in hashes_before.items()},
              "outputs": list(tables) + figure_paths + [f"bootstrap_{c['cohort']}.npz" for c in cohorts],
              "versions": {"python": platform.python_version(), "numpy": np.__version__,
                           "pandas": pd.__version__, "scipy": scipy.__version__,
                           "scikit_learn": sklearn.__version__, "matplotlib": matplotlib.__version__}}
    lines = ["# Historical cocktail reporting controls", "",
             "Retrospective analysis of saved scores and already known assay outcomes. No model was fitted.", "",
             "## Discrimination and prevalence", "",
             "Intervals resample the same historical host clusters for all methods. AP minus prevalence is recalculated within each draw.", "",
             "| Cohort | Method | AUROC (95% interval) | AP | Prevalence | AP minus prevalence (95% interval) |",
             "| --- | --- | --- | --- | --- | --- |"]
    for row in metrics.itertuples():
        lines.append(f"| {row.cohort} | {row.label} | {row.auroc:.3f} ({row.auroc_ci_low:.3f}–{row.auroc_ci_high:.3f}) | {row.average_precision:.3f} | {row.positive_prevalence:.3f} | {row.ap_minus_prevalence:.3f} ({row.ap_minus_prevalence_ci_low:.3f}–{row.ap_minus_prevalence_ci_high:.3f}) |")
    lines += ["", "## Historical-score cohort triage", "",
              "These are rank positions among already tested host–cocktail records; they are not validated laboratory cost savings.", "",
              "| Cohort | Target recall | Positive records required | Score-ranked records screened | Expected random-order rank |",
              "| --- | --- | --- | --- | --- |"]
    for row in recovery[recovery.method.eq("historical_product")].itertuples():
        lines.append(f"| {row.cohort} | {row.target_recall:.0%} | {row.positive_records_required} | {row.observed_rank_to_recover} | {row.expected_random_rank_to_recover:.2f} |")
    lines += ["", "## Limits", ""] + ["- " + text for text in limitations]
    lines += ["", "Exact per-method intervals, all 21 paired contrasts in each cohort, top-k results, rank-by-rank records, bootstrap weights, source hashes and training baseline counts are saved alongside this report.", ""]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    dump_json(OUT / "summary.json", report)
    print(metrics[["cohort", "method", "auroc", "average_precision", "positive_prevalence", "ap_minus_prevalence"]].to_string(index=False))
    print(json.dumps({"checks": checks, "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
