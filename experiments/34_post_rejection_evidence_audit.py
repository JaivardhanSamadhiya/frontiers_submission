"""Retrospective audit of the frozen cocktail evidence; never tunes a model.

All external outcomes were previously inspected. These analyses are diagnostic,
not a new confirmatory external validation. Frozen inputs are read-only.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/results_v2/external_cocktail_validation"
OUT = ROOT / "data/results_v2/post_rejection_audit"
UPSTREAM = ROOT / "external/gaborieau_coli_phage_interactions_2023"
SEED = 20261005


def metric(y, p, weights=None):
    return {
        "auroc": float(roc_auc_score(y, p, sample_weight=weights)),
        "average_precision": float(average_precision_score(y, p, sample_weight=weights)),
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cocktails = pd.read_csv(SOURCE / "heldout_cocktail_predictions_and_outcomes.csv")
    pairs = pd.read_csv(SOURCE / "heldout_selected_pair_predictions_and_outcomes.csv")
    all_pairs = pd.read_csv(SOURCE / "frozen_pair_predictions.csv")
    summary = json.loads((SOURCE / "external_cocktail_validation_summary.json").read_text())
    assert len(cocktails) == 100 and len(pairs) == 300 and len(all_pairs) == 9600
    assert hashlib.sha256((SOURCE / "frozen_pair_predictions.csv").read_bytes()).hexdigest() == summary["frozen_prediction_sha256"]
    y = cocktails.productive.to_numpy()
    names = cocktails.bacteria.tolist()
    pivot = all_pairs.pivot(index="bacteria", columns="phage", values="score_ensemble").loc[names]
    host_mean = pivot.mean(axis=1).to_numpy()
    scores = {
        "frozen_cocktail": cocktails.score_independence.to_numpy(),
        "max_constituent": cocktails.score_max.to_numpy(),
        "training_generalist": cocktails.score_generalist.to_numpy(),
        "host_mean_all_96_phages": host_mean,
    }
    for key in ["score_logistic", "score_knn"]:
        selected = pairs.pivot(index="bacteria", columns="phage", values=key).loc[names]
        scores[key + "_cocktail"] = 1 - np.nanprod(1 - selected.to_numpy(), axis=1)

    distance = pd.read_csv(UPSTREAM / "dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv", sep="\t", index_col=0)
    available = [n for n in names if n in distance.index]
    _, component = connected_components(csr_matrix(distance.loc[available, available].to_numpy() <= 1e-4), directed=False)
    mapping = dict(zip(available, component))
    nxt = int(max(component)) + 1
    for name in names:
        if name not in mapping:
            mapping[name] = nxt
            nxt += 1
    groups = np.array([mapping[n] for n in names])
    unique, gi = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(SEED)
    replicates = {k: [] for k in scores}
    attempts = 0
    while len(replicates["frozen_cocktail"]) < 2000:
        attempts += 1
        assert attempts < 100000
        counts = np.bincount(rng.integers(len(unique), size=len(unique)), minlength=len(unique))
        w = counts[gi]
        if np.unique(y[w > 0]).size < 2:
            continue
        for key, score in scores.items():
            m = metric(y, score, w)
            replicates[key].append([m["auroc"], m["average_precision"]])
    rows = []
    base = np.asarray(replicates["frozen_cocktail"])
    for key, score in scores.items():
        boot = np.asarray(replicates[key])
        lo, hi = np.quantile(boot, [0.025, 0.975], axis=0)
        delta = base - boot
        dlo, dhi = np.quantile(delta, [0.025, 0.975], axis=0)
        m = metric(y, score)
        rows.append({"score": key, **m, "auroc_ci_low": lo[0], "auroc_ci_high": hi[0],
                     "ap_ci_low": lo[1], "ap_ci_high": hi[1],
                     "frozen_minus_this_auroc": metric(y, scores["frozen_cocktail"])["auroc"] - m["auroc"],
                     "paired_auroc_delta_ci_low": dlo[0], "paired_auroc_delta_ci_high": dhi[0],
                     "paired_ap_delta_ci_low": dlo[1], "paired_ap_delta_ci_high": dhi[1]})

    # Break the chosen host--composition association, while retaining each
    # host's full predicted susceptibility profile and the observed outcomes.
    # This is a diagnostic distribution, not an exchangeability-based p-value.
    idx = {p: i for i, p in enumerate(pivot.columns)}
    composition = np.array([[idx[p] for p in s.split("|")] for s in cocktails.predicted_constituents])
    matrix = pivot.to_numpy()
    shuffled = []
    for _ in range(1000):
        shuffled_composition = composition[rng.permutation(len(names))]
        score = 1 - np.prod(1 - np.take_along_axis(matrix, shuffled_composition, axis=1), axis=1)
        shuffled.append(metric(y, score))

    # A host-specific choice requires within-host ordering, beyond an easy-host
    # effect. Only three externally selected constituents per host were assayed.
    within = []
    for name, part in pairs.groupby("bacteria", sort=False):
        if part.productive.nunique() == 2:
            m = metric(part.productive, part.score_ensemble)
            within.append({"bacteria": name, "n": len(part), **m})
    within = pd.DataFrame(within)
    rows = pd.DataFrame(rows)
    rows.to_csv(OUT / "cocktail_paired_controls.csv", index=False)
    within.to_csv(OUT / "selected_constituent_within_host.csv", index=False)
    pd.DataFrame(shuffled).to_csv(OUT / "composition_shuffle_diagnostic.csv", index=False)
    report = {
        "status": "retrospective diagnostic; external outcomes already known",
        "seed": SEED, "n_cocktails": len(y), "positive_prevalence": float(y.mean()),
        "n_host_clusters": len(unique), "bootstrap_replicates": 2000,
        "bootstrap_scope": "conditional on the fixed bank; paired host-cluster resampling",
        "average_precision_null_reference": float(y.mean()),
        "controls": rows.to_dict("records"),
        "composition_shuffle": {key: {
            "median": float(np.median([r[key] for r in shuffled])),
            "diagnostic_95_range": np.quantile([r[key] for r in shuffled], [0.025, 0.975]).tolist()
        } for key in ["auroc", "average_precision"]},
        "within_host_selected_constituents": {
            "evaluable_hosts": len(within), "macro_auroc": float(within.auroc.mean()),
            "interpretation": "restricted to the three published selections; not full-bank or new-cocktail validation"
        },
        "deployment_issue": "Per-phage percentile normalization uses the other test hosts. A singleton batch gives each constituent 0.5 and a three-phage score 0.875. The frozen rank score is cohort-dependent, not an individual-host predictor.",
        "input_checksums": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [
            SOURCE / "frozen_pair_predictions.csv", SOURCE / "heldout_cocktail_predictions_and_outcomes.csv",
            SOURCE / "heldout_selected_pair_predictions_and_outcomes.csv"]},
    }
    (OUT / "cocktail_evidence_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(rows.round(4).to_string(index=False))
    print(json.dumps({k:report[k] for k in ["composition_shuffle", "within_host_selected_constituents", "deployment_issue"]}, indent=2))


if __name__ == "__main__":
    main()
