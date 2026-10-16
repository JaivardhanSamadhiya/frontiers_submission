"""Retrospective nested evaluation of raw kNN selection and simple breadth.

This evaluates observed constituent coverage, not measured mixture efficacy.
Family choice was informed by previous development analyses. Original outputs
are not changed. Missing selected assays are reported with coverage bounds.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[1]
UPSTREAM = ROOT / "external/gaborieau_coli_phage_interactions_2023"
OUT = ROOT / "data/results_v2/post_rejection_audit"
GRID = [5, 10, 20, 40]
THRESHOLDS = [1e-4, 1e-3, 1e-2]
SEED = 20261005


def nearest_scores(tr, te, y, observed, distance, k):
    result = np.empty((len(te), y.shape[1]))
    dd = distance[np.ix_(te, tr)]
    for j in range(y.shape[1]):
        keep = observed[tr, j]
        prior = float(y[tr[keep], j].mean()) if keep.any() else 0.0
        if not keep.any():
            result[:, j] = prior
            continue
        values = y[tr[keep], j]
        d = dd[:, keep]
        order = np.argsort(d, axis=1, kind="stable")[:, :min(k, len(values))]
        chosen_d = np.take_along_axis(d, order, axis=1)
        finite = np.isfinite(chosen_d)
        weights = np.where(finite, 1 / np.maximum(chosen_d, 1e-6), 0)
        denominator = weights.sum(axis=1)
        result[:, j] = np.divide((weights * values[order]).sum(axis=1), denominator,
                                out=np.full(len(te), prior), where=denominator > 0)
    return result


def selection(scores, genus, isolation):
    selected = []
    for j in np.argsort(-scores, kind="stable"):
        if not any(genus[j] == genus[a] and isolation[j] == isolation[a] for a in selected):
            selected.append(j)
        if len(selected) == 3:
            return np.array(selected)
    raise RuntimeError("Insufficient diverse candidate phages")


def outcome(host, idx, y, observed):
    known_positive = bool((y[host, idx] * observed[host, idx]).any())
    complete = bool(observed[host, idx].all())
    return int(known_positive), int(known_positive or not complete), float((y[host, idx] * observed[host, idx]).sum() / 3), complete


def tune(tr, y, observed, groups, distance, genus, isolation):
    inner = GroupKFold(n_splits=min(3, len(np.unique(groups[tr]))))
    results = {k:[] for k in GRID}
    for ti, vi in inner.split(tr, groups=groups[tr]):
        inner_tr, inner_te = tr[ti], tr[vi]
        for k in GRID:
            scores = nearest_scores(inner_tr, inner_te, y, observed, distance, k)
            for local, host in enumerate(inner_te):
                chosen = selection(scores[local], genus, isolation)
                lower, upper, frac, complete = outcome(host, chosen, y, observed)
                results[k].append((lower, frac))
    # Optimize the conservative coverage lower bound over ALL inner held-out
    # hosts. Unknown assays never disappear from the denominator.
    return max(GRID, key=lambda k:(np.mean([v[0] for v in results[k]]),
                                  np.mean([v[1] for v in results[k]]), -k)), {
        k:{"lower_coverage":float(np.mean([v[0] for v in values])),
           "positive_fraction_lower":float(np.mean([v[1] for v in values]))}
        for k,values in results.items()}


def paired_bootstrap(frame, n=2000):
    unique, gi = np.unique(frame.group, return_inverse=True)
    rng = np.random.default_rng(SEED)
    columns = ["design_lower", "design_upper", "baseline_lower", "baseline_upper"]
    values = frame[columns].to_numpy()
    boot = []
    for _ in range(n):
        counts = np.bincount(rng.integers(len(unique),size=len(unique)),minlength=len(unique))
        mean = np.average(values,axis=0,weights=counts[gi])
        boot.append([mean[0]-mean[2],mean[0]-mean[3],mean[1]-mean[2]])
    arr=np.asarray(boot)
    return {"paired_lower_bound_difference_ci":np.quantile(arr[:,0],[.025,.975]).tolist(),
            "missingness_robust_difference_interval_with_sampling":
                [float(np.quantile(arr[:,1],.025)),float(np.quantile(arr[:,2],.975))],
            "replicates":n,"unit":"host sequence component; conditional on fixed bank"}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    paths=[UPSTREAM / "data/interactions/interaction_matrix.csv",
           UPSTREAM / "data/genomics/phages/guelin_collection.csv",
           UPSTREAM / "dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv"]
    matrix=pd.read_csv(paths[0],sep=";",index_col=0)
    metadata=pd.read_csv(paths[1],sep=";").set_index("phage").loc[matrix.columns]
    names=matrix.index.astype(str).tolist()
    distances=pd.read_csv(paths[2],sep="\t",index_col=0).loc[names,names].to_numpy(float)
    assert np.allclose(distances,distances.T,equal_nan=True)
    assert np.isfinite(distances).all()
    y=(matrix.fillna(0).to_numpy()>0).astype(int)
    observed=~matrix.isna().to_numpy()
    genus=metadata.Genus.fillna("Unknown").to_numpy()
    isolation=metadata.Phage_host.fillna("Unknown").to_numpy()
    rows,summaries,settings=[],[],[]
    for threshold in THRESHOLDS:
        _,groups=connected_components(csr_matrix(distances<=threshold),directed=False)
        count=len(np.unique(groups))
        if count<5:
            summaries.append({"threshold":threshold,"n_groups":count,"status":"skipped: fewer than five groups"})
            continue
        outer=GroupKFold(n_splits=5)
        threshold_rows=[]
        for fold,(tr,te) in enumerate(outer.split(y,groups=groups)):
            if len(np.unique(groups[tr]))<2:
                raise RuntimeError("Insufficient inner training groups")
            k,inner=tune(tr,y,observed,groups,distances,genus,isolation)
            scores=nearest_scores(tr,te,y,observed,distances,k)
            numerator=(y[tr]*observed[tr]).sum(axis=0)
            denominator=observed[tr].sum(axis=0)
            prevalence=np.divide(numerator,denominator,out=np.zeros(y.shape[1]),where=denominator>0)
            base=selection(prevalence,genus,isolation)
            settings.append({"threshold":threshold,"fold":fold,"k":k,"inner":inner})
            assert set(groups[tr]).isdisjoint(groups[te])
            for local,host in enumerate(te):
                design=selection(scores[local],genus,isolation)
                dl,du,df,dc=outcome(host,design,y,observed)
                bl,bu,bf,bc=outcome(host,base,y,observed)
                threshold_rows.append({"threshold":threshold,"fold":fold,"group":int(groups[host]),
                    "bacteria":names[host],"k":k,"design_phages":"|".join(matrix.columns[design]),
                    "baseline_phages":"|".join(matrix.columns[base]),
                    "design_lower":dl,"design_upper":du,"baseline_lower":bl,"baseline_upper":bu,
                    "design_positive_fraction_lower":df,"baseline_positive_fraction_lower":bf,
                    "design_complete":dc,"baseline_complete":bc,
                    "nearest_train_distance":float(distances[host,tr].min())})
            print(f"threshold {threshold} fold {fold} k {k} done",flush=True)
        frame=pd.DataFrame(threshold_rows)
        assert len(frame)==402 and not frame.bacteria.duplicated().any()
        rows.extend(threshold_rows)
        summary={"threshold":threshold,"n_groups":count,"n_hosts":len(frame),
            "status":"completed exploratory nested evaluation",
            **{c:float(frame[c].mean()) for c in ["design_lower","design_upper","baseline_lower","baseline_upper",
                "design_positive_fraction_lower","baseline_positive_fraction_lower"]},
            "design_incomplete_n":int((~frame.design_complete).sum()),
            "baseline_incomplete_n":int((~frame.baseline_complete).sum()),
            "paired_bootstrap":paired_bootstrap(frame)}
        summaries.append(summary)
        print(json.dumps(summary,indent=2),flush=True)
    pd.DataFrame(rows).to_csv(OUT / "corrected_nested_selector_hosts.csv",index=False)
    report={"status":"retrospective; family chosen after prior development audit",
        "outcome":"at least one active constituent; no mixture efficacy measured",
        "settings":settings,"summaries":summaries,
        "input_checksums":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}
    (OUT / "corrected_nested_selector_summary.json").write_text(json.dumps(report,indent=2),encoding="utf-8")


if __name__=="__main__":
    main()
