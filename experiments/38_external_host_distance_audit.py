"""Retrospective external sensitivity to related development hosts.

No tuning or replacement of frozen scores. Direct-distance exclusion and
connected-component exclusion answer different questions; report both.
"""
from pathlib import Path
import json
import hashlib

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "external/gaborieau_coli_phage_interactions_2023"
OUT = ROOT / "data/results_v2/post_rejection_audit"


def estimates(frame, groups):
    y = frame.productive.to_numpy()
    if len(frame) < 2 or len(np.unique(y)) < 2:
        return {"status": "not estimable: fewer than two outcome classes"}
    unique, inverse = np.unique(groups, return_inverse=True)
    result = {}
    for key in ["score_independence", "development_reference_score"]:
        score = frame[key].to_numpy()
        boot = []
        rng = np.random.default_rng(20261005)
        for _ in range(20000):
            weight = np.bincount(rng.integers(len(unique), size=len(unique)), minlength=len(unique))[inverse]
            if len(np.unique(y[weight > 0])) < 2:
                continue
            boot.append([roc_auc_score(y, score, sample_weight=weight),
                         average_precision_score(y, score, sample_weight=weight)])
            if len(boot) == 2000:
                break
        result[key] = {"auroc": float(roc_auc_score(y, score)),
                       "average_precision": float(average_precision_score(y, score)),
                       "valid_bootstrap_replicates": len(boot)}
        if boot:
            result[key]["auroc_ci"] = np.quantile(np.array(boot)[:, 0], [.025, .975]).tolist()
            result[key]["ap_ci"] = np.quantile(np.array(boot)[:, 1], [.025, .975]).tolist()
    return result


def main():
    paths = [UPSTREAM / "data/interactions/interaction_matrix.csv",
             UPSTREAM / "dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv",
             OUT / "batch_invariant_cocktail_diagnostic.csv"]
    training = pd.read_csv(paths[0], sep=";", index_col=0).index.astype(str).tolist()
    distance = pd.read_csv(paths[1], sep="\t", index_col=0)
    frame = pd.read_csv(paths[2])
    external = frame.bacteria.astype(str).tolist()
    collisions = sorted(set(training) & set(external))
    # Short names are not guaranteed globally unique. Do not silently map a
    # different external Gembase record to a development host distance row.
    training_meta = pd.read_csv(UPSTREAM / "data/genomics/bacteria/picard_collection.csv", sep=";").set_index("bacteria")
    test_meta = pd.read_csv(UPSTREAM / "dev/cocktails/data/test_collection.csv", sep=";").set_index("bacteria")
    collision_records = [{"name": name, "development_genome": str(training_meta.loc[name, "Gembase"]),
                          "external_genome": str(test_meta.loc[name, "Gembase"])} for name in collisions]
    available = [name for name in external if name in distance.index and name in distance.columns]
    unavailable = sorted(set(external) - set(available))
    names = training + [name for name in available if name not in set(training)]
    dd = distance.loc[names, names].to_numpy(float)
    assert np.isfinite(dd).all() and np.allclose(dd, dd.T)
    nearest = np.array([float(distance.loc[name, training].min()) if name in available else np.nan for name in external])
    frame["nearest_development_host_distance"] = nearest
    unambiguous = ~frame.bacteria.isin(collisions).to_numpy()
    all_external_distances = distance.loc[available, available].to_numpy(float)
    def groups_at(threshold):
        _, labels = connected_components(csr_matrix(all_external_distances <= threshold), directed=False)
        lookup = dict(zip(available, labels))
        for i, name in enumerate(unavailable):
            lookup[name] = int(labels.max()) + 1 + i
        return np.array([lookup[name] for name in external])
    original_groups = groups_at(1e-4)
    summaries = [{"threshold": 1e-4, "exclusion_rule": "exclude_identifier_collisions_only",
                  "n_hosts": int(unambiguous.sum()), "n_positive": int(frame.loc[unambiguous].productive.sum()),
                  "n_negative": int((1 - frame.loc[unambiguous].productive).sum()),
                  "n_host_groups": int(len(np.unique(original_groups[unambiguous]))),
                  "results": estimates(frame.loc[unambiguous], original_groups[unambiguous])}]
    for threshold in [1e-4, 1e-3, 1e-2]:
        _, joint_groups = connected_components(csr_matrix(dd <= threshold), directed=False)
        external_groups = groups_at(threshold)
        lookup = dict(zip(names, joint_groups))
        shared = np.isin([lookup.get(name, -1) for name in external], joint_groups[:len(training)])
        frame[f"shares_development_component_{threshold:g}"] = shared
        for rule, eligible in [("nearest_distance_greater_than_threshold", (nearest > threshold) & unambiguous),
                               ("no_shared_development_component", ~shared & unambiguous & np.isfinite(nearest))]:
            subset = frame.loc[eligible].copy()
            summaries.append({"threshold": threshold, "exclusion_rule": rule,
                              "n_hosts": len(subset), "n_positive": int(subset.productive.sum()),
                              "n_negative": int((1 - subset.productive).sum()),
                              "n_host_groups": int(len(np.unique(external_groups[eligible]))),
                              "results": estimates(subset, external_groups[eligible])})
    frame.to_csv(OUT / "external_host_distance_diagnostic.csv", index=False)
    report = {"status": "retrospective diagnostic; external outcomes previously observed",
              "distance": "published core-genome host distance; not ANI",
              "identifier_collisions": collision_records,
              "external_hosts_without_distance": unavailable,
              "missing_distance_policy": "retain in collision-only analysis as singleton bootstrap units; exclude from distance sensitivity",
              "collision_policy": "retain historical score; exclude ambiguous identifiers in every new sensitivity",
              "nearest_distance_quantiles": dict(zip(["min", "q25", "median", "q75", "max"],
                                                      np.quantile(nearest[np.isfinite(nearest)], [0, .25, .5, .75, 1]).tolist())),
              "summaries": summaries,
              "input_checksums": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (OUT / "external_host_distance_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
