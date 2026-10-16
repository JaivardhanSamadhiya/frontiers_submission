"""Retrospective, training-only audit of exact-identifier degree controls.

This extends rather than replaces the original count-degree audit. No sequence
model is fit, and no score is selected from these held-out outcomes.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from precisionphage.splits import combined_unseen_folds, leave_one_group_out

OUT = ROOT / "data/results_v2/frontiers_revision_2026_10_15_followup/degree_control"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def score_controls(train: pd.DataFrame, test: pd.DataFrame):
    """Score using only training identifiers and labels; test labels unused.

    Legacy positive counts retain their original zero fallback and maximum
    normalization. Positive fractions use all assayed rows as denominators.
    Exact IDs absent from training receive the training-set label prevalence,
    rather than being mistaken for assayed-but-inactive entities.
    """
    if len(train) == 0:
        raise ValueError("At least one training pair is required")
    prior = float(train["label"].mean())
    details = pd.DataFrame(index=np.arange(len(test)))
    values = {}
    for axis in ["phage", "host"]:
        total = train.groupby(axis).size()
        positive = train.groupby(axis)["label"].sum().reindex(total.index)
        n = test[axis].map(total).fillna(0).to_numpy(float)
        pos = test[axis].map(positive).fillna(0).to_numpy(float)
        fraction = np.full(len(test), prior, dtype=float)
        np.divide(pos, n, out=fraction, where=n > 0)
        details[axis + "_train_observed"] = n.astype(int)
        details[axis + "_train_positive"] = pos.astype(int)
        details[axis + "_seen"] = n > 0
        details[axis + "_positive_fraction"] = fraction
        values[axis + "_positive_degree"] = pos
        values[axis + "_positive_fraction"] = fraction
        values[axis + "_normalised_degree"] = pos / max(1, float(positive.max()))
    values["additive_positive_degree"] = (
        values.pop("phage_normalised_degree") + values.pop("host_normalised_degree")) / 2
    values["additive_positive_fraction"] = (
        values["phage_positive_fraction"] + values["host_positive_fraction"]) / 2
    # Diagnostic: within each fold this is only a shift of the additive score.
    # It exposes effects of pooling distinct training prevalences, not a tuned
    # second predictor. The two-axis average has prior as its natural centre.
    values["centred_additive_positive_fraction_diagnostic"] = (
        values["additive_positive_fraction"] - prior)
    values["training_prior_diagnostic"] = np.full(len(test), prior)
    details["train_prior"] = prior
    return values, details


def metric_row(labels, score):
    labels, score = np.asarray(labels), np.asarray(score)
    return {"pairs": len(labels), "positive": int(labels.sum()),
            "negative": int((labels == 0).sum()), "prevalence": float(labels.mean()),
            "auroc": float(roc_auc_score(labels, score)),
            "average_precision": float(average_precision_score(labels, score))}


def _write_csv(frame, name):
    frame.to_csv(OUT / name, index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model_path = ROOT / "data/interim_v2/interactions_modelable.csv"
    cluster_path = ROOT / "data/cache_v2/clusters_k21_d0.05.json"
    old_path = ROOT / "data/results_v2/frontiers_revision_2026_10_14/reporting_audit/degree_predictions.csv"
    pooled_path = ROOT / "data/results_v2/significance_pooled_preds.npz"
    plan = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "Retrospective diagnostic prompted by review; original outcomes already known.",
        "no_fitting_or_tuning": True,
        "regimes": ["loso_species", "logo_genus", "host_cluster", "phage_cluster", "combined_unseen"],
        "legacy_methods": "Phage, host and equally weighted maximum-normalised positive training counts; unknown zero.",
        "added_controls": "Phage, host and equally weighted positive training fractions: positive count / all observed training pairs per exact entity; unseen entity falls back to global training prevalence.",
        "diagnostics": "Fold training prior alone, and additive fraction minus that prior, separate pooled-prior shifts from within-fold rankings. Neither is selected as a winner.",
        "metrics": "AUROC/AP/prevalence pooled; AUROC/AP by fold; macro fold AUROC and pair-count weighted within-fold AUROC. Score distributions and seen/unseen status by outcome are diagnostic, not additional predictors.",
        "within_fold_weight": "n_positive * n_negative; pooling only comparable positive-negative pairs from the same fold removes cross-fold ranking pairs.",
        "independence": "No independent-row p-values, confidence intervals or causal-mechanism claims.",
        "test_labels": "Used for reporting only; score_controls does not inspect them.",
        "inputs": {p.relative_to(ROOT).as_posix(): sha(p) for p in [model_path, cluster_path, old_path, pooled_path]},
    }
    plan_path = OUT / "retrospective_plan.json"
    if plan_path.exists():
        retained = json.loads(plan_path.read_text())
        # The saved Windows plan is retained verbatim; canonicalise separators
        # only for verification so the same release can be rerun on POSIX.
        retained_inputs = {k.replace("\\", "/"): v for k, v in retained["inputs"].items()}
        if retained_inputs != plan["inputs"]:
            raise RuntimeError("Inputs changed after saved retrospective plan")
    else:
        plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")

    cov = pd.read_csv(model_path)
    if len(cov) != 1947:
        raise ValueError("Unexpected preserved sequence-covered cohort")
    clusters = json.loads(cluster_path.read_text())
    for axis in ["phage", "host"]:
        cov[axis + "_cluster"] = cov[axis].map(clusters[axis])
    assert not cov[["phage_cluster", "host_cluster"]].isna().any().any()
    design = {
        "loso_species": list(leave_one_group_out(cov, "host_species", "loso", 3, 3)),
        "logo_genus": list(leave_one_group_out(cov, "host_genus", "logo", 3, 3)),
        "host_cluster": list(leave_one_group_out(cov, "host_cluster", "host_cluster", 3, 3)),
        "phage_cluster": list(leave_one_group_out(cov, "phage_cluster", "phage_cluster", 3, 3)),
        "combined_unseen": list(combined_unseen_folds(cov, "phage_cluster", "host_cluster", 5, 42, 3, 3)),
    }
    predictions, fold_rows, row_details, entity_rows = [], [], [], []
    saved = np.load(pooled_path, allow_pickle=False)
    for regime, folds in design.items():
        offset = 0
        for fold in folds:
            train, test = cov.iloc[fold.train_idx], cov.iloc[fold.test_idx]
            values, details = score_controls(train, test)
            labels = test.label.to_numpy()
            if regime + "__GBM" in saved:
                np.testing.assert_array_equal(saved[regime + "__y"][offset:offset + len(test)], labels)
                values["saved_full_gbm_reference"] = saved[regime + "__GBM"][offset:offset + len(test)]
            offset += len(test)
            base = test[["phage", "host", "label"]].reset_index(drop=True)
            base["row_index"] = fold.test_idx
            base["regime"], base["fold"] = regime, fold.name
            row_details.append(pd.concat([base, details], axis=1))
            for method, score in values.items():
                block = base.copy()
                block["method"], block["score"] = method, score
                predictions.append(block)
                fold_rows.append({"regime": regime, "fold": fold.name, "method": method,
                                  **metric_row(labels, score), "train_prior": details.train_prior.iloc[0],
                                  "unique_scores": len(np.unique(score))})
            for axis in ["phage", "host"]:
                table = train.groupby(axis).label.agg(["size", "sum"])
                for entity, subset in test.groupby(axis):
                    n = int(table.loc[entity, "size"]) if entity in table.index else 0
                    positive = int(table.loc[entity, "sum"]) if entity in table.index else 0
                    entity_rows.append({"regime": regime, "fold": fold.name, "axis": axis,
                                        "entity": entity, "train_observed": n, "train_positive": positive,
                                        "test_observed": len(subset), "test_positive": int(subset.label.sum())})
    predictions = pd.concat(predictions, ignore_index=True)
    fold_metrics = pd.DataFrame(fold_rows)
    details = pd.concat(row_details, ignore_index=True)
    # Exact reproduction guard: no old score or held-out cohort has been changed.
    historical = pd.read_csv(old_path)
    key = ["regime", "fold", "method", "row_index"]
    match = historical.merge(predictions, on=key, suffixes=("_old", "_new"), validate="one_to_one")
    assert len(match) == len(historical)
    np.testing.assert_array_equal(match.label_old, match.label_new)
    np.testing.assert_allclose(match.score_old, match.score_new, rtol=0, atol=1e-15)
    summary_rows = []
    for (regime, method), group in predictions.groupby(["regime", "method"]):
        folds = fold_metrics[(fold_metrics.regime == regime) & (fold_metrics.method == method)]
        weights = folds.positive * folds.negative
        summary_rows.append({"regime": regime, "method": method, **metric_row(group.label, group.score),
                             "folds": len(folds), "macro_fold_auroc": float(folds.auroc.mean()),
                             "within_fold_pair_weighted_auroc": float(np.average(folds.auroc, weights=weights)),
                             "min_fold_auroc": float(folds.auroc.min()), "max_fold_auroc": float(folds.auroc.max()),
                             "folds_below_half": int((folds.auroc < .5).sum()),
                             "folds_equal_zero": int((folds.auroc == 0).sum()),
                             "folds_equal_half": int((folds.auroc == .5).sum())})
    summary = pd.DataFrame(summary_rows)
    status_rows = []
    for (regime, label), group in details.groupby(["regime", "label"]):
        for axis in ["phage", "host"]:
            seen = group[axis + "_seen"]
            status_rows.append({"regime": regime, "label": int(label), "axis": axis,
                                "pairs": len(group), "seen_pairs": int(seen.sum()),
                                "seen_fraction": float(seen.mean()),
                                "nonzero_positive_degree_pairs": int((group[axis + "_train_positive"] > 0).sum()),
                                "mean_training_observed_degree": float(group[axis + "_train_observed"].mean()),
                                "mean_training_positive_degree": float(group[axis + "_train_positive"].mean()),
                                "mean_positive_fraction_score": float(group[axis + "_positive_fraction"].mean())})
    score_stats = predictions.groupby(["regime", "fold", "method", "label"]).score.agg(
        ["size", "mean", "median", "min", "max"]).reset_index()
    _write_csv(predictions, "predictions.csv")
    _write_csv(fold_metrics, "fold_metrics.csv")
    _write_csv(summary, "summary.csv")
    _write_csv(details, "row_training_degrees.csv")
    _write_csv(pd.DataFrame(entity_rows), "entity_train_test_counts.csv")
    _write_csv(pd.DataFrame(status_rows), "seen_status_by_outcome.csv")
    _write_csv(score_stats, "fold_score_distributions.csv")
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(),
                "legacy_rows_reproduced": len(historical),
                "script_sha256": sha(Path(__file__)),
                "files": {p.name: sha(p) for p in OUT.glob("*.csv")}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
