"""Repair cohort-dependent ranking using Picard-only reference distributions.

External comparison is retrospective, not replacement of the frozen result.
The feature/model parameters are inherited without external tuning.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from precisionphage.cocktail.scoring import DevelopmentRankNormalizer, combined_ranking_score

spec = importlib.util.spec_from_file_location("frozen_helper", ROOT / "experiments/29_freeze_external_cocktail_predictions.py")
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
OUT = ROOT / "data/results_v2/post_rejection_audit"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = ROOT / "data/results_v2/external_cocktail_validation"
    settings = json.loads((source / "frozen_prediction_manifest.json").read_text())["selected"]
    matrix = pd.read_csv(helper.INTERACTIONS, sep=";", index_col=0)
    names = matrix.index.astype(str).tolist()
    phages = matrix.columns.astype(str).tolist()
    meta = helper.clean_features(pd.read_csv(helper.TRAIN_META, sep=";")).set_index("bacteria").loc[names].reset_index()
    distance = pd.read_csv(helper.DISTANCES, sep="\t", index_col=0)
    groups = helper.build_groups(distance, names)
    y = (matrix.fillna(0).to_numpy() > 0).astype(int)
    observed = ~matrix.isna().to_numpy()
    log_ref = np.full(y.shape, np.nan)
    knn_ref = np.full(y.shape, np.nan)
    for fold, (tr, te) in enumerate(GroupKFold(10).split(meta, groups=groups)):
        for j in range(len(phages)):
            use = tr[observed[tr,j]]
            log_ref[te,j] = helper.constant_or_logistic(meta.iloc[use], y[use,j], meta.iloc[te], settings["logistic_C"])
            knn_ref[te,j] = helper.knn_predict([names[i] for i in use], y[use,j], [names[i] for i in te], distance, settings["knn_k"])
        print(f"development reference fold {fold} complete",flush=True)
    assert np.isfinite(log_ref).all() and np.isfinite(knn_ref).all()
    log_norm = DevelopmentRankNormalizer().fit(log_ref)
    knn_norm = DevelopmentRankNormalizer().fit(knn_ref)
    frozen = pd.read_csv(source / "frozen_pair_predictions.csv")
    external = pd.read_csv(source / "heldout_cocktail_predictions_and_outcomes.csv")
    test_names = external.bacteria.tolist()
    raw_log = frozen.pivot(index="bacteria", columns="phage", values="prob_logistic_raw").loc[test_names,phages].to_numpy()
    raw_knn = frozen.pivot(index="bacteria", columns="phage", values="prob_knn_raw").loc[test_names,phages].to_numpy()
    new = .5 * log_norm.transform(raw_log) + .5 * knn_norm.transform(raw_knn)
    # Every single-host output must equal the corresponding batch output.
    for i in range(len(test_names)):
        np.testing.assert_array_equal(new[i], (.5 * log_norm.transform(raw_log[i:i+1]) + .5 * knn_norm.transform(raw_knn[i:i+1]))[0])
    lookup = {p:i for i,p in enumerate(phages)}
    idx = np.array([[lookup[p] for p in s.split("|")] for s in external.predicted_constituents])
    external["development_reference_score"] = combined_ranking_score(np.take_along_axis(new,idx,axis=1))
    external["raw_mean_score"] = combined_ranking_score(np.take_along_axis(.5*raw_log+.5*raw_knn,idx,axis=1))
    results={}
    for column in ["score_independence", "development_reference_score", "raw_mean_score"]:
        results[column]={"auroc":float(roc_auc_score(external.productive,external[column])),
            "average_precision":float(average_precision_score(external.productive,external[column]))}
    np.savez_compressed(OUT / "development_rank_reference.npz", logistic=log_ref, knn=knn_ref,
                        phages=np.asarray(phages), hosts=np.asarray(names))
    external.to_csv(OUT / "batch_invariant_cocktail_diagnostic.csv",index=False)
    (OUT / "batch_invariant_scoring_summary.json").write_text(json.dumps({
        "status":"retrospective diagnostic; outcomes already inspected",
        "parameters":settings,"reference":"10-fold host-grouped Picard OOF raw scores",
        "reference_hyperparameters":"selected earlier using all Picard CV; not a nested unbiased CV estimate",
        "singleton_batch_consistency":"all 100 hosts passed exact equality",
        "scores_are_probabilities":False,"results":results},indent=2),encoding="utf-8")
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
