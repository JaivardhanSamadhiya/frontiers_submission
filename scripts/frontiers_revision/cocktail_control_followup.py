"""Audit the below-chance source-only phage-marginal control.

This is a retrospective diagnostic of fixed scores and known outcomes. It fits
no predictor, changes no historical score, and does not identify a causal reason
for source cocktail assignment. Outputs are separate from the original results.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/results_v2/frontiers_revision_2026_10_15_followup/cocktail_control"
INPUTS = {
    "host_scores": ROOT / "data/results_v2/frontiers_revision_2026_10_14/cocktail_reporting/host_scores.csv",
    "selected_pair_scores": ROOT / "data/results_v2/external_cocktail_validation/heldout_selected_pair_predictions_and_outcomes.csv",
    "training_matrix": ROOT / "external/gaborieau_coli_phage_interactions_2023/data/interactions/interaction_matrix.csv",
    "source_selection_steps": ROOT / "external/gaborieau_coli_phage_interactions_2023/dev/cocktails/data/cocktails_step.csv",
    "source_selection_documentation": ROOT / "external/gaborieau_coli_phage_interactions_2023/README.md",
    "crossed_metrics": ROOT / "data/results_v2/frontiers_revision_2026_10_15/crossed_cocktail/crossed_metrics.csv",
    "gbm_folds": ROOT / "data/results_v2/frontiers_revision_2026_10_14/reporting_audit/verified_gbm_fold_metrics.csv",
    "fold_rows": ROOT / "data/results_v2/frontiers_revision_2026_10_14/reporting_audit/degree_predictions.csv",
    "cohorts": ROOT / "data/results_v2/frontiers_revision_2026_10_14/reporting_audit/evaluation_cohorts.csv",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def composition_key(value: str) -> str:
    names = value.split("|")
    if len(names) != 3 or len(set(names)) != 3:
        raise ValueError("Expected three distinct phages")
    return "|".join(sorted(names))


def main() -> None:
    before = {name: sha256(path) for name, path in INPUTS.items()}
    hosts = pd.read_csv(INPUTS["host_scores"])
    pairs = pd.read_csv(INPUTS["selected_pair_scores"])
    matrix = pd.read_csv(INPUTS["training_matrix"], sep=";", index_col=0)
    steps = pd.read_csv(INPUTS["source_selection_steps"], sep=";")
    intervals = pd.read_csv(INPUTS["crossed_metrics"])
    folds = pd.read_csv(INPUTS["gbm_folds"])
    fold_rows = pd.read_csv(INPUTS["fold_rows"])
    cohorts = pd.read_csv(INPUTS["cohorts"])
    source_readme = INPUTS["source_selection_documentation"].read_text(encoding="utf-8")
    assert "most **generalist phages**" in source_readme and "fewer than three phages are still recommended" in source_readme
    assert hosts.bacteria.is_unique and len(hosts) == 100
    assert len(pairs) == len(steps) == 300
    assert matrix.shape == (402, 96)
    assert int(hosts.productive.sum()) == 82
    assert not pairs.duplicated(["bacteria", "phage"]).any()
    assert not steps.duplicated(["bacteria", "phage"]).any()
    observed = matrix.notna()
    positives = matrix.gt(0) & observed
    reference = pd.DataFrame({
        "training_positive": positives.sum(),
        "training_observed": observed.sum(),
    })
    reference.index.name = "phage"
    reference["training_positive_fraction"] = reference.training_positive / reference.training_observed
    expected_pair_score = reference.training_positive_fraction.loc[pairs.phage].to_numpy(float)
    np.testing.assert_allclose(expected_pair_score, pairs.training_prevalence, atol=1e-15, rtol=0)
    expected_host_score = pairs.groupby("bacteria").training_prevalence.mean().loc[hosts.bacteria]
    np.testing.assert_allclose(expected_host_score, hosts.training_phage_mean, atol=1e-15, rtol=0)
    hosts["composition_key"] = hosts.predicted_constituents.map(composition_key)
    unique_counts = hosts.groupby("composition_key").training_phage_mean.nunique()
    assert unique_counts.eq(1).all()
    assert hosts.composition_key.nunique() == 18

    merged = pairs.merge(steps[["bacteria", "phage", "chosen_at_step"]],
                         on=["bacteria", "phage"], how="left", validate="one_to_one")
    assert merged.chosen_at_step.notna().all()
    patterns = merged.groupby("bacteria").chosen_at_step.apply(lambda x: "|".join(sorted(x)))
    hosts["selection_step_pattern"] = patterns.loc[hosts.bacteria].to_numpy()
    hosts["rescue_constituents"] = hosts.bacteria.map(
        merged.assign(rescue=merged.chosen_at_step.eq("4_rescue")).groupby("bacteria").rescue.sum())

    by_outcome = hosts.groupby("productive").training_phage_mean.agg(
        records="size", mean="mean", median="median", minimum="min", maximum="max").reset_index()
    by_composition = hosts.groupby("composition_key").agg(
        records=("productive", "size"), positive=("productive", "sum"),
        training_phage_mean=("training_phage_mean", "first"),
        unique_scores=("training_phage_mean", "nunique"),
        host_groups=("host_cluster", "nunique"),
    ).reset_index()
    by_composition["negative"] = by_composition.records - by_composition.positive
    by_composition["observed_activity_fraction"] = by_composition.positive / by_composition.records
    by_composition = by_composition.sort_values(["records", "composition_key"], ascending=[False, True])
    by_stage = hosts.groupby("selection_step_pattern").agg(
        records=("productive", "size"), positive=("productive", "sum"),
        mean_training_score=("training_phage_mean", "mean"),
        compositions=("composition_key", "nunique"),
    ).reset_index()
    by_stage["observed_activity_fraction"] = by_stage.positive / by_stage.records
    by_rescue = hosts.groupby("rescue_constituents").agg(
        records=("productive", "size"), positive=("productive", "sum"),
        mean_training_score=("training_phage_mean", "mean"),
    ).reset_index()
    by_rescue["observed_activity_fraction"] = by_rescue.positive / by_rescue.records
    rescue_presence = hosts.assign(any_rescue=hosts.rescue_constituents.gt(0)).groupby("any_rescue").agg(
        records=("productive", "size"), positive=("productive", "sum"),
        mean_training_score=("training_phage_mean", "mean")).reset_index()
    rescue_presence["observed_activity_fraction"] = rescue_presence.positive / rescue_presence.records
    zero_folds = folds.loc[folds.auroc.eq(0)].copy()
    zero_folds["negative"] = zero_folds.n - zero_folds.positive
    zero_folds["n_positive_negative_pairs"] = zero_folds.positive * zero_folds.negative
    zero_row_sets = []
    for row in zero_folds.itertuples():
        part = fold_rows.loc[(fold_rows.regime == row.regime) &
                             (fold_rows.fold.astype(str) == str(row.fold)) &
                             (fold_rows.method == "phage_positive_degree")]
        assert len(part) == 9 and int(part.label.sum()) == 5
        assert set(part.host) == {"leuconostoc pseudomesenteroides"}
        zero_row_sets.append(set(part.row_index))
    assert len(zero_row_sets) == 2 and zero_row_sets[0] == zero_row_sets[1]

    control = intervals.loc[(intervals.cohort == "historical_100") &
                            (intervals.method == "bank_phage_prevalence_mean")].iloc[0]
    auroc = float(roc_auc_score(hosts.productive, hosts.training_phage_mean))
    average_precision = float(average_precision_score(hosts.productive, hosts.training_phage_mean))
    np.testing.assert_allclose(auroc, control.auroc, atol=1e-14, rtol=0)
    np.testing.assert_allclose(average_precision, control.average_precision, atol=1e-14, rtol=0)
    assert (folds.positive >= 3).all() and ((folds.n - folds.positive) >= 3).all()
    only_negative_compositions = int(((by_composition.negative > 0) & (by_composition.positive == 0)).sum())
    after = {name: sha256(path) for name, path in INPUTS.items()}
    assert before == after
    report = {
        "status": "completed retrospective diagnostic; no model fitted",
        "score_definition": "Arithmetic mean over the three chosen phages of each phage's positive fraction among non-missing individual Picard assays; zero source cocktail outcomes or target-host outcomes are used to calculate this score.",
        "source_only_score_verified": True,
        "records": len(hosts), "positive": int(hosts.productive.sum()),
        "compositions": int(hosts.composition_key.nunique()),
        "composition_constant_score_verified": bool(unique_counts.eq(1).all()),
        "negative_bearing_compositions": int(by_composition.negative.gt(0).sum()),
        "only_negative_compositions": only_negative_compositions,
        "auroc": auroc,
        "auroc_crossed_ci_low": float(control.auroc_ci_low),
        "auroc_crossed_ci_high": float(control.auroc_ci_high),
        "average_precision": average_precision,
        "prevalence": float(hosts.productive.mean()),
        "ap_minus_prevalence": average_precision - float(hosts.productive.mean()),
        "score_reversal_not_used": "The control direction was fixed as greater source susceptibility = greater predicted activity. It is not flipped after seeing the target result.",
        "interpretation": "In this selected cohort, inactive cocktails have higher mean training susceptibility than active cocktails. The control has only one value per composition and no host-specific term. Its crossed interval includes chance; a robust inverse relation is not established.",
        "causal_limit": "The released records describe selection stages and repeated compositions, but do not identify the causal reason for the inverse association. They do not establish that broadly active phages were deliberately assigned to harder hosts.",
        "documented_selection_recipe": "The upstream README section 2.6 states that the rescue step fills unfilled slots using generalist phages with high overall lytic coverage after earlier stages recommend fewer than three phages. This is documentary evidence of the selection recipe, not evidence that the target hosts were intrinsically harder or that the recipe caused the observed inverse score association.",
        "any_rescue_observed_summary": rescue_presence.to_dict(orient="records"),
        "orphan_genus_row": {
            "purpose": "Additional eligible-genus count audit generated in submission_audit.py, not a performance regime with saved predictions.",
            "recommendation": "Remove this row from main Table 1 to match Methods, Table 2 and Figure 1A; preserve count audit in Data S1.",
            "counts": cohorts.loc[cohorts.regime.eq("logo_genus")].to_dict(orient="records"),
        },
        "zero_auroc_folds": zero_folds.to_dict(orient="records"),
        "zero_fold_explanation": "The LOSO and host-cluster zero-AUROC folds each contain the same nine Leuconostoc pseudomesenteroides rows (5 positive, 4 negative); all 20 positive-negative score comparisons are reversed. Small eligible folds permit extreme AUROC values; eligibility requires only three observations of each class. This is not a missing-class calculation.",
        "input_hashes_unchanged": before == after,
        "input_hashes": {name: {"path": str(INPUTS[name].relative_to(ROOT)), "sha256": value}
                         for name, value in before.items()},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scikit_learn": sklearn.__version__},
    }
    OUT.mkdir(parents=True, exist_ok=True)
    for name, frame in {
        "host_diagnostic.csv": hosts,
        "training_phage_reference.csv": reference.reset_index(),
        "score_by_outcome.csv": by_outcome,
        "composition_diagnostic.csv": by_composition,
        "selection_stage_diagnostic.csv": by_stage,
        "rescue_count_diagnostic.csv": by_rescue,
        "rescue_presence_diagnostic.csv": rescue_presence,
        "zero_auroc_folds.csv": zero_folds,
    }.items():
        frame.to_csv(OUT / name, index=False)
    (OUT / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    lines = [
        "# Training phage susceptibility control audit", "",
        report["score_definition"], "",
        f"The 100 records contain 82 active cocktails and 18 distinct compositions. The score is constant within each composition. AUROC is {auroc:.6f}, with crossed host-cluster by composition bootstrap interval {control.auroc_ci_low:.6f}–{control.auroc_ci_high:.6f}. The interval includes 0.5.", "",
        "## Observed score distribution", "",
        "| Observed cocktail outcome | Records | Mean score | Median score |",
        "| --- | --- | --- | --- |",
    ]
    for row in by_outcome.itertuples():
        lines.append(f"| {int(row.productive)} | {int(row.records)} | {row.mean:.6f} | {row.median:.6f} |")
    lines += ["", report["interpretation"], "", report["documented_selection_recipe"], "", report["causal_limit"], "",
              "No control direction is reversed or retuned using the target outcomes. Stage and menu summaries are descriptive, not adjustment models or causal analyses.", "",
              "## Related reporting checks", "",
              report["orphan_genus_row"]["recommendation"], "",
              report["zero_fold_explanation"], "",
              "All input file SHA-256 hashes were unchanged. Exact training counts, menu summaries, released selection-stage patterns, and zero-AUROC folds accompany this report.", ""]
    (OUT / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({key: report[key] for key in ["status", "auroc", "auroc_crossed_ci_low", "auroc_crossed_ci_high", "compositions", "negative_bearing_compositions", "input_hashes_unchanged"]}, indent=2))


if __name__ == "__main__":
    main()
