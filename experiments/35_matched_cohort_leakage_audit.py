"""Fixed-test, class-count-matched retrospective split comparison on CPU."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from precisionphage.features.genomic import edge_features_from_spectra, kmer_spectrum
from precisionphage.splits.cluster import combined_unseen_folds

OUT = ROOT / "data/results_v2/post_rejection_audit"
SEEDS = [20261005, 20261006, 20261007]


def load_design():
    df = pd.read_csv(ROOT / "data/interim_v2/interactions_modelable.csv")
    assert len(df) == 1947 and df.label.sum() == 1488
    phages, hosts = sorted(df.phage.unique()), sorted(df.host.unique())
    name_hash = hashlib.sha256(("\n".join(phages) + "\n--HOSTS--\n" + "\n".join(hosts)).encode()).hexdigest()
    meta = json.loads((ROOT / "data/cache_v2/node_features.meta.json").read_text())
    assert name_hash == meta["names_sha256"]
    z = np.load(ROOT / "data/cache_v2/node_features.npz")
    pi = df.phage.map(dict(zip(phages, range(len(phages))))).to_numpy()
    hi = df.host.map(dict(zip(hosts, range(len(hosts))))).to_numpy()
    pp, hh = z["phage"][pi], z["host"][hi]
    kdim = len(kmer_spectrum("ACGT" * 4, 4))
    local = np.array([edge_features_from_spectra(p[:kdim], h[:kdim]) for p, h in zip(pp, hh)])
    extra = pd.read_csv(ROOT / "data/interim_v2/seq_pair_features.csv")
    merged = df[["phage", "host"]].merge(extra, on=["phage", "host"], how="left", validate="one_to_one")
    extra_array = merged.drop(columns=["phage", "host"]).to_numpy(np.float32)
    assert np.isfinite(extra_array).all()
    xx = np.hstack([pp, hh, local, df[["k3dist", "k6dist", "GCdiff", "Homology"]].to_numpy(), extra_array]).astype(np.float32)
    clusters = json.loads((ROOT / "data/cache_v2/clusters_k21_d0.05.json").read_text())
    df["phage_cluster"] = df.phage.map(clusters["phage"])
    df["host_cluster"] = df.host.map(clusters["host"])
    assert not df[["phage_cluster", "host_cluster"]].isna().any().any()
    return df, xx


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, xx = load_design()
    yy = df.label.to_numpy()
    folds = list(combined_unseen_folds(df, "phage_cluster", "host_cluster", 5, 42))
    assert len(folds) == 5
    predictions, folds_report = [], []
    regimes = ["shared_entities_allowed", "host_cluster_excluded", "phage_cluster_excluded", "both_clusters_excluded"]
    for fi, fold in enumerate(folds):
        te = fold.test_idx
        host_overlap = df.host_cluster.isin(df.iloc[te].host_cluster.unique()).to_numpy()
        phage_overlap = df.phage_cluster.isin(df.iloc[te].phage_cluster.unique()).to_numpy()
        not_test = ~np.isin(np.arange(len(df)), te)
        masks = [not_test, not_test & ~host_overlap, not_test & ~phage_overlap,
                 not_test & ~host_overlap & ~phage_overlap]
        counts = {c: min(int((mask & (yy == c)).sum()) for mask in masks) for c in [0, 1]}
        for seed in SEEDS:
            for ri, (regime, mask) in enumerate(zip(regimes, masks)):
                rng = np.random.default_rng(seed + fi * 100 + ri)
                tr = np.concatenate([rng.choice(np.flatnonzero(mask & (yy == c)), counts[c], replace=False) for c in [0, 1]])
                assert not np.intersect1d(tr, te).size
                if ri in [1, 3]:
                    assert set(df.iloc[tr].host_cluster).isdisjoint(df.iloc[te].host_cluster)
                if ri in [2, 3]:
                    assert set(df.iloc[tr].phage_cluster).isdisjoint(df.iloc[te].phage_cluster)
                scaler = StandardScaler().fit(xx[tr])
                model = XGBClassifier(n_estimators=400, max_depth=5, learning_rate=0.05,
                    subsample=0.8, colsample_bytree=0.8, random_state=seed,
                    tree_method="hist", n_jobs=2, eval_metric="logloss", verbosity=0)
                model.fit(scaler.transform(xx[tr]), yy[tr])
                score = model.predict_proba(scaler.transform(xx[te]))[:, 1]
                frame = df.iloc[te][["phage", "host", "label", "phage_cluster", "host_cluster"]].copy()
                frame["row_index"], frame["fold"], frame["seed"], frame["regime"], frame["score"] = te, fi, seed, regime, score
                predictions.append(frame)
                folds_report.append({"fold": fi, "seed": seed, "regime": regime,
                    "test_n": len(te), "test_positive": int(yy[te].sum()),
                    "train_n": len(tr), "train_positive": counts[1],
                    "auroc": float(roc_auc_score(yy[te], score)),
                    "average_precision": float(average_precision_score(yy[te], score)),
                    "test_host_cluster_seen_fraction": float(df.iloc[te].host_cluster.isin(df.iloc[tr].host_cluster).mean()),
                    "test_phage_cluster_seen_fraction": float(df.iloc[te].phage_cluster.isin(df.iloc[tr].phage_cluster).mean())})
            print(f"block {fi} seed {seed} complete", flush=True)
    pred = pd.concat(predictions, ignore_index=True)
    result = []
    for (regime, seed), part in pred.groupby(["regime", "seed"]):
        assert not part.row_index.duplicated().any()
        result.append({"regime": regime, "seed": int(seed), "n": len(part),
            "positive": int(part.label.sum()), "auroc": float(roc_auc_score(part.label, part.score)),
            "average_precision": float(average_precision_score(part.label, part.score))})
    pooled = pd.DataFrame(result)
    pred.to_csv(OUT / "matched_cohort_predictions.csv", index=False)
    pd.DataFrame(folds_report).to_csv(OUT / "matched_cohort_folds.csv", index=False)
    pooled.to_csv(OUT / "matched_cohort_summary.csv", index=False)
    input_paths = [ROOT / "data/interim_v2/interactions_modelable.csv", ROOT / "data/cache_v2/node_features.npz",
        ROOT / "data/cache_v2/node_features.meta.json", ROOT / "data/interim_v2/seq_pair_features.csv",
        ROOT / "data/cache_v2/clusters_k21_d0.05.json", ROOT / "POST_REJECTION_ANALYSIS_PROTOCOL.md"]
    (OUT / "matched_cohort_manifest.json").write_text(json.dumps({
        "status": "retrospective; shared test rows and matched training class counts",
        "feature_dimension": xx.shape[1], "seeds": SEEDS, "xgboost_trees": 400,
        "checksums": {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in input_paths},
        "limitations": ["five blocks", "seeds are sensitivity runs, not independent samples",
            "fixed configuration inherited from prior explored analysis", "original forward-strand clusters retained"]
    }, indent=2), encoding="utf-8")
    print(pooled.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
