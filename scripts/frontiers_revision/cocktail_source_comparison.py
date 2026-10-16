"""Retrospective controls for the Frontiers revision; never refit a predictor.

The published selected and fixed-baseline triplets have measured outcomes. The
released source recommender also records 233 of 300 selected constituent scores.
Its 67 unscored rescue selections are missing, not zero. All comparative metrics
are computed on explicitly matched hosts. Historical and later repaired scores
remain separate. Intervals resample whole host sequence groups and condition on
the fitted predictors, fixed phage bank, and externally selected compositions.
"""
from __future__ import annotations

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[2]
UPSTREAM = ROOT / "external/gaborieau_coli_phage_interactions_2023"
FROZEN = ROOT / "data/results_v2/external_cocktail_validation"
AUDIT = ROOT / "data/results_v2/post_rejection_audit"
OUT = ROOT / "data/results_v2/frontiers_revision_2026_10_14/cocktail_controls"
SEED = 20261014
REPLICATES = 2000
THRESHOLD = 1e-4
ALIASES = {"BCH857": "BCH857-72"}
COLLISIONS = {"AN03"}
BASELINE = tuple(sorted(["LF82_P8", "CLB_P2", "LF73_P1"]))
PHAGE_ALIASES = {"CLB_P2": "55989_P2"}
BASELINE_BANK = tuple(sorted([PHAGE_ALIASES.get(p, p) for p in BASELINE]))
SOURCE_COMMIT = "99656d8c020fcd1dbbc30a02ee10274a07655562"
SOURCE_DOI = "https://doi.org/10.1038/s41564-024-01832-5"
SOURCE_MANUSCRIPT = (
    "https://hal.science/hal-04781608v1/file/"
    "Gaborieau_et_al_Accepted%20Manuscript.pdf"
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_names(frame):
    frame = frame.copy()
    frame["bacteria"] = frame.bacteria.replace(ALIASES)
    return frame


def product_score(values):
    values = np.asarray(values, float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("Supply a finite two-dimensional constituent-score table")
    if ((values < 0) | (values > 1)).any():
        raise ValueError("Constituent scores must be in [0, 1]")
    return 1 - np.prod(1 - values, axis=1)


def verify_baseline_alias():
    """Require exact genome identity, rather than inferring an alias by name."""
    url = "https://www.ebi.ac.uk/ena/browser/api/fasta/OZ035781.1?download=false"
    local = UPSTREAM / "data/genomics/phages/FNA/55989_P2.fna"
    cached = OUT / "ENA_OZ035781_1.fasta"
    if not cached.exists():
        request = Request(url, headers={"User-Agent": "PrecisionPhage-reproducibility-audit/1.0"})
        with urlopen(request, timeout=90) as response:
            data = response.read()
        if not data.startswith(b">"):
            raise ValueError("ENA did not return a FASTA record")
        cached.write_bytes(data)

    def sequence(path):
        return "".join(line.strip().upper() for line in path.read_text().splitlines()
                       if not line.startswith(">"))

    known, candidate = sequence(cached), sequence(local)
    exact = known == candidate
    prefix_redundancy = (len(candidate) == len(known) + 60 and candidate[:len(known)] == known
                         and candidate[len(known):] == known[:60])
    if not (exact or prefix_redundancy):
        raise ValueError("CLB_P2 accession identity could not be established for bank 55989_P2")
    return {"published_name": "CLB_P2", "bank_identifier": "55989_P2",
            "primary_accession": "OZ035781.1", "primary_registry_url": url,
            "published_base_pairs": len(known), "local_base_pairs": len(candidate),
            "exact_whole_sequence_identity": exact,
            "complete_accession_sequence_identical": True,
            "local_terminal_duplicate_of_accession_prefix_bp": 60 if prefix_redundancy else 0,
            "accession_sequence_sha256": hashlib.sha256(known.encode()).hexdigest(),
            "local_sequence_sha256": hashlib.sha256(candidate.encode()).hexdigest(),
            "downloaded_record_sha256": sha(cached), "local_fasta_sha256": sha(local),
            "downloaded_record_path": str(cached), "local_fasta_path": str(local)}


class RankedMetric:
    """Fast tied-score, weighted AUROC/AP; checked against scikit-learn below."""

    def __init__(self, labels, scores):
        labels, scores = np.asarray(labels, int), np.asarray(scores, float)
        if labels.shape != scores.shape or not np.isfinite(scores).all():
            raise ValueError("Invalid score/label arrays")
        self.order = np.argsort(-scores, kind="stable")
        self.labels = labels[self.order]
        sorted_score = scores[self.order]
        self.ends = np.r_[np.flatnonzero(np.diff(sorted_score)), len(labels) - 1]

    def __call__(self, weights):
        weights = np.asarray(weights, float)[self.order]
        tp = np.cumsum(weights * self.labels)[self.ends]
        fp = np.cumsum(weights * (1 - self.labels))[self.ends]
        if tp[-1] <= 0 or fp[-1] <= 0:
            return None
        recall, fpr = tp / tp[-1], fp / fp[-1]
        precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=tp + fp > 0)
        ap = np.sum(np.diff(np.r_[0.0, recall]) * precision)
        auc = np.trapezoid(np.r_[0.0, recall], np.r_[0.0, fpr])
        prevalence = tp[-1] / (tp[-1] + fp[-1])
        return np.array([auc, ap, prevalence, ap - prevalence])


def verify_metric_implementation():
    rng = np.random.default_rng(42)
    checked = 0
    for _ in range(200):
        y = rng.integers(0, 2, 40)
        score = rng.integers(0, 8, 40) / 7
        weights = rng.integers(0, 4, 40)
        if (weights[y == 1].sum() == 0) or (weights[y == 0].sum() == 0):
            continue
        observed = RankedMetric(y, score)(weights)
        expected = [roc_auc_score(y, score, sample_weight=weights),
                    average_precision_score(y, score, sample_weight=weights)]
        np.testing.assert_allclose(observed[:2], expected, atol=1e-14, rtol=1e-14)
        checked += 1
    constant = RankedMetric([0, 1, 1], [.5, .5, .5])(np.ones(3))
    np.testing.assert_allclose(constant, [.5, 2 / 3, 2 / 3, 0])
    if RankedMetric([0, 1], [.1, .9])([0, 1]) is not None:
        raise AssertionError("A one-class resample was not rejected")
    return {"weighted_tied_score_checks_against_sklearn": checked,
            "constant_score_and_single_class_checks": "passed"}


def historical_groups(hosts, distances):
    """Reconstruct the exact historical host-group rule; save assignments."""
    available = [host for host in hosts if host in distances.index]
    square = distances.loc[available, available].to_numpy(float)
    _, labels = connected_components(csr_matrix(square <= THRESHOLD), directed=False)
    mapping = dict(zip(available, labels))
    next_id = int(max(labels)) + 1
    for host in hosts:
        if host not in mapping:
            mapping[host] = next_id
            next_id += 1
    return np.asarray([mapping[host] for host in hosts], int)


def bootstrap_weights(groups, count=REPLICATES, seed=SEED):
    unique, inverse = np.unique(groups, return_inverse=True)
    rng = np.random.default_rng(seed)
    for _ in range(count):
        counts = np.bincount(rng.integers(len(unique), size=len(unique)),
                             minlength=len(unique))
        yield counts[inverse].astype(float)


def intervals(values):
    return np.quantile(np.asarray(values), [.025, .975], axis=0)


def summarize_cohort(frame, score_columns, cohort, reference="pp_historical_product"):
    y = frame.observed_lytic.to_numpy(int)
    groups = frame.host_group.to_numpy(int)
    evaluators = {key: RankedMetric(y, frame[key]) for key in score_columns}
    unit = np.ones(len(frame))
    points = {key: evaluator(unit) for key, evaluator in evaluators.items()}
    if any(value is None for value in points.values()):
        raise ValueError(f"Cohort {cohort} has fewer than two outcome classes")
    replicates = {key: [] for key in score_columns}
    valid = 0
    attempts = 0
    # Use identical resampled host groups for all methods within a cohort.
    for weights in bootstrap_weights(groups, REPLICATES * 50):
        attempts += 1
        if len(np.unique(y[weights > 0])) < 2:
            continue
        for key, evaluator in evaluators.items():
            replicates[key].append(evaluator(weights))
        valid += 1
        if valid == REPLICATES:
            break
    if valid != REPLICATES:
        raise RuntimeError(f"Only {valid} valid bootstrap resamples for {cohort}")
    metric_names = ["auroc", "average_precision", "prevalence", "ap_minus_prevalence"]
    rows = []
    for key in score_columns:
        low, high = intervals(replicates[key])
        delta = np.asarray(replicates[reference]) - np.asarray(replicates[key])
        dlow, dhigh = intervals(delta)
        repaired_key = "pp_development_reference_product"
        repaired_delta = np.asarray(replicates[repaired_key]) - np.asarray(replicates[key])
        repaired_low, repaired_high = intervals(repaired_delta)
        row = {"cohort": cohort, "score": key, "hosts": len(frame),
               "positive": int(y.sum()), "negative": int((1-y).sum()),
               "host_groups": len(np.unique(groups)), "bootstrap_valid": valid,
               "bootstrap_attempts": attempts}
        for j, name in enumerate(metric_names):
            row[name] = float(points[key][j])
            row[name + "_ci_low"] = float(low[j])
            row[name + "_ci_high"] = float(high[j])
            row["historical_minus_this_" + name] = float(points[reference][j] - points[key][j])
            row["historical_minus_this_" + name + "_ci_low"] = float(dlow[j])
            row["historical_minus_this_" + name + "_ci_high"] = float(dhigh[j])
            row["development_reference_minus_this_" + name] = float(points[repaired_key][j] - points[key][j])
            row["development_reference_minus_this_" + name + "_ci_low"] = float(repaired_low[j])
            row["development_reference_minus_this_" + name + "_ci_high"] = float(repaired_high[j])
        rows.append(row)
    return rows


def reference_percentiles(raw_scores, references):
    raw_scores, references = np.asarray(raw_scores), np.asarray(references)
    if raw_scores.shape[1] != references.shape[1]:
        raise ValueError("Reference phage order differs")
    result = np.empty_like(raw_scores)
    n = references.shape[0]
    for j in range(raw_scores.shape[1]):
        ref = np.sort(references[:, j])
        left = np.searchsorted(ref, raw_scores[:, j], side="left")
        right = np.searchsorted(ref, raw_scores[:, j], side="right")
        result[:, j] = (left + right + 1) / (2 * (n + 1))
    return result


def preference_summary(frame, score_pairs, cohort):
    selected_y = frame.selected_lytic.to_numpy(float)
    baseline_y = frame.baseline_lytic.to_numpy(float)
    discordant = selected_y != baseline_y
    selected_mlc = frame.selected_mlc.to_numpy(float)
    baseline_mlc = frame.baseline_mlc.to_numpy(float)
    ordinal_differs = selected_mlc != baseline_mlc
    groups = frame.host_group.to_numpy(int)
    rows, host_rows = [], []
    for key, (selected_column, baseline_column) in score_pairs.items():
        sd = frame[selected_column].to_numpy(float) - frame[baseline_column].to_numpy(float)
        choice = np.where(sd > 0, 1.0, np.where(sd < 0, 0.0, .5))
        observed = choice * selected_y + (1 - choice) * baseline_y
        correct = np.where(selected_y > baseline_y, choice, 1 - choice)
        ordinal_correct = np.where(selected_mlc > baseline_mlc, choice, 1 - choice)

        def summarize(weights):
            denom = weights.sum()
            relevant = np.sum(weights[discordant])
            ordinal_n = np.sum(weights[ordinal_differs])
            policy = np.sum(weights * observed) / denom
            return np.array([
                np.sum(weights[discordant] * correct[discordant]) / relevant
                if relevant else np.nan,
                policy,
                policy - np.sum(weights * selected_y) / denom,
                policy - np.sum(weights * baseline_y) / denom,
                policy - np.sum(weights * (selected_y + baseline_y) / 2) / denom,
                np.sum(weights[ordinal_differs] * ordinal_correct[ordinal_differs]) / ordinal_n
                if ordinal_n else np.nan,
            ])

        point = summarize(np.ones(len(frame)))
        samples = np.array([summarize(w) for w in bootstrap_weights(groups)])
        low, high = np.nanquantile(samples, [.025, .975], axis=0)
        names = ["binary_discordant_preference_accuracy", "policy_observed_lytic_fraction",
                 "policy_minus_always_selected", "policy_minus_fixed_baseline",
                 "policy_minus_uniform_random", "ordinal_discordant_preference_accuracy"]
        row = {"cohort": cohort, "score": key, "hosts": len(frame),
               "binary_discordant_hosts": int(discordant.sum()),
               "ordinal_discordant_hosts": int(ordinal_differs.sum()),
               "host_groups": int(len(np.unique(groups))),
               "selected_only_lytic": int(((selected_y == 1) & (baseline_y == 0)).sum()),
               "baseline_only_lytic": int(((selected_y == 0) & (baseline_y == 1)).sum()),
               "both_lytic": int(((selected_y == 1) & (baseline_y == 1)).sum()),
               "neither_lytic": int(((selected_y == 0) & (baseline_y == 0)).sum()),
               "always_selected_lytic_fraction": float(selected_y.mean()),
               "fixed_baseline_lytic_fraction": float(baseline_y.mean()),
               "oracle_lytic_fraction": float(np.maximum(selected_y, baseline_y).mean()),
               "uniform_random_lytic_fraction": float(((selected_y + baseline_y) / 2).mean()),
               "ties": int((sd == 0).sum()), "bootstrap_replicates": REPLICATES,
               "bootstrap_valid_preference_replicates": int(np.isfinite(samples[:, 0]).sum())}
        for j, name in enumerate(names):
            row[name] = float(point[j])
            row[name + "_ci_low"] = float(low[j])
            row[name + "_ci_high"] = float(high[j])
        rows.append(row)
        for i, host in enumerate(frame.bacteria):
            host_rows.append({"cohort": cohort, "bacteria": host, "score": key,
                              "host_group": int(groups[i]), "score_difference": float(sd[i]),
                              "binary_discordant": bool(discordant[i]),
                              "selected_lytic": int(selected_y[i]),
                              "baseline_lytic": int(baseline_y[i]),
                              "probability_select_selected_if_tie_random": float(choice[i]),
                              "policy_observed_lytic": float(observed[i])})
    return rows, host_rows


def markdown_table(rows, columns):
    result = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    for row in rows:
        result.append("| " + " | ".join(str(row[column]) for column in columns) + " |")
    return "\n".join(result)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    baseline_alias_evidence = verify_baseline_alias()
    source_paths = {
        "source_step_scores": UPSTREAM / "dev/cocktails/data/cocktails_step.csv",
        "source_composition": UPSTREAM / "dev/cocktails/data/cocktails_composition.csv",
        "source_outcomes": UPSTREAM / "dev/cocktails/data/cml_cocktails.csv",
        "source_combined_selected_table": UPSTREAM / "dev/cocktails/data/cocktails.csv",
        "source_distances": UPSTREAM / "dev/cocktails/data/picard+test_collection_phylogenetic_distances.tsv",
        "source_train_metadata": UPSTREAM / "data/genomics/bacteria/picard_collection.csv",
        "source_test_metadata": UPSTREAM / "dev/cocktails/data/test_collection.csv",
        "source_interactions": UPSTREAM / "data/interactions/interaction_matrix.csv",
        "frozen_pair_predictions": FROZEN / "frozen_pair_predictions.csv",
        "frozen_cocktail_predictions": FROZEN / "heldout_cocktail_predictions_and_outcomes.csv",
        "frozen_manifest": FROZEN / "frozen_prediction_manifest.json",
        "historical_evaluation": FROZEN / "external_cocktail_validation_summary.json",
        "development_reference": AUDIT / "development_rank_reference.npz",
        "repaired_scores": AUDIT / "batch_invariant_cocktail_diagnostic.csv",
        "previous_controls": AUDIT / "cocktail_evidence_audit.json",
        "baseline_alias_published_genome": OUT / "ENA_OZ035781_1.fasta",
        "baseline_alias_local_genome": UPSTREAM / "data/genomics/phages/FNA/55989_P2.fna",
        "analysis_script": Path(__file__),
    }
    hashes_before = {key: sha(path) for key, path in source_paths.items()}
    verification = verify_metric_implementation()
    history = json.loads(source_paths["historical_evaluation"].read_text())
    manifest = json.loads(source_paths["frozen_manifest"].read_text())
    expected_hash = manifest["checksums"]["frozen_pair_predictions"]
    assert hashes_before["frozen_pair_predictions"] == expected_hash
    assert history["frozen_prediction_sha256"] == expected_hash

    source = canonical_names(pd.read_csv(source_paths["source_step_scores"], sep=";"))
    composition = canonical_names(pd.read_csv(source_paths["source_composition"], sep=";"))
    outcomes = canonical_names(pd.read_csv(source_paths["source_outcomes"], sep=";"))
    published = canonical_names(pd.read_csv(source_paths["source_combined_selected_table"], sep=";"))
    cocktail = pd.read_csv(source_paths["frozen_cocktail_predictions"])
    pair = pd.read_csv(source_paths["frozen_pair_predictions"])
    assert len(pair) == 9600 and pair[["bacteria", "phage"]].duplicated().sum() == 0
    assert len(source) == 300 and source[["bacteria", "phage"]].duplicated().sum() == 0
    assert source.bacteria.nunique() == 100 and source.groupby("bacteria").size().eq(3).all()
    assert len(composition) == len(outcomes) == len(cocktail) == 100
    assert all(x.bacteria.is_unique for x in [composition, outcomes, cocktail, published])
    assert set(composition.bacteria) == set(outcomes.bacteria) == set(pair.bacteria)
    composition = composition.set_index("bacteria").loc[cocktail.bacteria].reset_index()
    outcomes = outcomes.set_index("bacteria").loc[cocktail.bacteria].reset_index()
    published = published.set_index("bacteria").loc[cocktail.bacteria].reset_index()
    np.testing.assert_array_equal(published["p1+p2+p3_MLC"], outcomes["p1+p2+p3"])
    np.testing.assert_array_equal(cocktail.productive, (outcomes["p1+p2+p3"] > 0).astype(int))
    selected_names = composition[["p1", "p2", "p3"]].to_numpy(str)
    assert all(len(set(row)) == 3 for row in selected_names)
    for host, phages in zip(cocktail.bacteria, selected_names):
        assert set(source.loc[source.bacteria == host, "phage"]) == set(phages)

    distance = pd.read_csv(source_paths["source_distances"], sep="\t", index_col=0)
    groups = historical_groups(cocktail.bacteria.tolist(), distance)
    previous = json.loads(source_paths["previous_controls"].read_text())
    assert len(np.unique(groups)) == previous["n_host_clusters"] == 81
    group_table = pd.DataFrame({"bacteria": cocktail.bacteria, "host_group": groups,
                                "distance_record_available": cocktail.bacteria.isin(distance.index),
                                "identifier_ambiguous": cocktail.bacteria.isin(COLLISIONS)})
    group_table.to_csv(OUT / "host_group_assignments.csv", index=False)

    references = np.load(source_paths["development_reference"], allow_pickle=False)
    phages = references["phages"].astype(str).tolist()
    host_names = cocktail.bacteria.tolist()
    raw_log = pair.pivot(index="bacteria", columns="phage", values="prob_logistic_raw").loc[host_names, phages].to_numpy()
    raw_knn = pair.pivot(index="bacteria", columns="phage", values="prob_knn_raw").loc[host_names, phages].to_numpy()
    historical = pair.pivot(index="bacteria", columns="phage", values="score_ensemble").loc[host_names, phages].to_numpy()
    prevalence = pair.pivot(index="bacteria", columns="phage", values="training_prevalence").loc[host_names, phages].to_numpy()
    repaired = .5 * reference_percentiles(raw_log, references["logistic"]) + .5 * reference_percentiles(raw_knn, references["knn"])
    lookup = {phage: j for j, phage in enumerate(phages)}
    selected_indices = np.array([[lookup[p] for p in row] for row in selected_names])
    baseline_indices = np.tile([lookup[p] for p in BASELINE_BANK], (100, 1))
    take = lambda matrix, indices: np.take_along_axis(matrix, indices, axis=1)
    selected_historical = take(historical, selected_indices)
    pp_historical = product_score(selected_historical)
    pp_repaired = product_score(take(repaired, selected_indices))
    np.testing.assert_allclose(pp_historical, cocktail.score_independence, rtol=0, atol=2e-15)
    saved_repaired = pd.read_csv(source_paths["repaired_scores"]).set_index("bacteria").loc[host_names]
    np.testing.assert_allclose(pp_repaired, saved_repaired.development_reference_score, rtol=0, atol=2e-15)
    # Existing references alone determine every query; check singleton equality.
    for i in range(100):
        one = (.5 * reference_percentiles(raw_log[i:i+1], references["logistic"])
               + .5 * reference_percentiles(raw_knn[i:i+1], references["knn"]))
        np.testing.assert_array_equal(one[0], repaired[i])
    verification["repaired_singleton_batch_equal_hosts"] = 100

    train_meta = pd.read_csv(source_paths["source_train_metadata"], sep=";").set_index("bacteria")
    test_meta = pd.read_csv(source_paths["source_test_metadata"], sep=";").set_index("bacteria").loc[host_names]
    matrix = pd.read_csv(source_paths["source_interactions"], sep=";", index_col=0)
    assert set(matrix.columns) == set(phages)
    phenotype = matrix.gt(0).where(matrix.notna())
    prior = float(np.nansum(phenotype.to_numpy()) / np.isfinite(phenotype.to_numpy(float)).sum())
    group_susceptibility = {}
    for phylogroup, names in train_meta.groupby("Clermont_Phylo").groups.items():
        names = sorted(set(names).intersection(matrix.index))
        values = phenotype.loc[names].to_numpy(float)
        group_susceptibility[str(phylogroup)] = float(np.nansum(values) / np.isfinite(values).sum())
    phylo = test_meta.Clermont_Phylo.astype(str).map(group_susceptibility).fillna(prior).to_numpy(float)
    # This phylogroup control uses only Picard labels, averaged over the bank.
    # Both candidate triplets receive the same host-only score by construction.
    source_pivot = source.pivot(index="bacteria", columns="phage", values="y_pred_proba").reindex(index=host_names, columns=phages)
    original_selected = take(source_pivot.to_numpy(float), selected_indices)
    available = np.isfinite(original_selected).all(axis=1)
    original_complete_score = np.full(100, np.nan)
    original_complete_score[available] = product_score(original_selected[available])

    selected = pd.DataFrame({
        "bacteria": host_names, "host_group": groups,
        "observed_lytic": (outcomes["p1+p2+p3"] > 0).astype(int),
        "mlc_score": outcomes["p1+p2+p3"], "pp_historical_product": pp_historical,
        "pp_development_reference_product": pp_repaired,
        "pp_historical_max": selected_historical.max(axis=1),
        "pp_historical_mean": selected_historical.mean(axis=1),
        "pp_logistic_raw_product": product_score(take(raw_log, selected_indices)),
        "pp_knn_raw_product": product_score(take(raw_knn, selected_indices)),
        "bank_phage_prevalence_mean": take(prevalence, selected_indices).mean(axis=1),
        "bank_phage_prevalence_product": product_score(take(prevalence, selected_indices)),
        "training_phylogroup_susceptibility": phylo,
        "source_released_product": original_complete_score,
        "source_constituent_scores_available": np.isfinite(original_selected).sum(axis=1),
        "source_product_lower_per_host": 1 - np.prod(1 - np.nan_to_num(original_selected, nan=0.0), axis=1),
        "source_product_upper_per_host": 1 - np.prod(1 - np.nan_to_num(original_selected, nan=1.0), axis=1),
        "composition_key": ["|".join(sorted(row)) for row in selected_names],
        "identifier_ambiguous": [name in COLLISIONS for name in host_names],
    })
    assert int(available.sum()) == 64
    selected.to_csv(OUT / "selected_triplet_scores.csv", index=False)
    source.to_csv(OUT / "source_released_constituent_scores.csv", index=False)

    # Identification bounds, deliberately outcome-informed: unknown constituent
    # scores may be anywhere in [0,1]. These are not generated predictions.
    bound_rows = []
    for exclude_collision in [False, True]:
        bound_frame = selected.loc[~selected.identifier_ambiguous] if exclude_collision else selected
        by = bound_frame.observed_lytic.to_numpy(int)
        low = bound_frame.source_product_lower_per_host.to_numpy(float)
        high = bound_frame.source_product_upper_per_host.to_numpy(float)
        assert (low <= high).all()
        for case, score in [("worst", np.where(by == 1, low, high)),
                            ("best", np.where(by == 1, high, low))]:
            value = RankedMetric(by, score)(np.ones(len(by)))
            bound_rows.append({"cohort": "exclude_AN03_99" if exclude_collision else "historical_100",
                               "case": case, "hosts": len(by), "positive": int(by.sum()),
                               "auroc": value[0], "average_precision": value[1],
                               "prevalence": value[2], "ap_minus_prevalence": value[3],
                               "interpretation": "outcome-informed identification bound; not a model prediction or confidence interval"})
    pd.DataFrame(bound_rows).to_csv(OUT / "source_full_cohort_identification_bounds.csv", index=False)
    bound_checks = 0
    bound_rng = np.random.default_rng(SEED + 1)
    for _ in range(200):
        possible = original_selected.copy()
        missing = ~np.isfinite(possible)
        possible[missing] = bound_rng.random(missing.sum())
        possible = product_score(possible)
        for excluded in [False, True]:
            keep = ~selected.identifier_ambiguous.to_numpy() if excluded else np.ones(100, bool)
            yy = selected.observed_lytic.to_numpy()[keep]
            value = RankedMetric(yy, possible[keep])(np.ones(keep.sum()))
            cohort = "exclude_AN03_99" if excluded else "historical_100"
            bounds = [r for r in bound_rows if r["cohort"] == cohort]
            lower, upper = bounds[0], bounds[1]
            assert lower["auroc"] - 1e-14 <= value[0] <= upper["auroc"] + 1e-14
            assert lower["average_precision"] - 1e-14 <= value[1] <= upper["average_precision"] + 1e-14
            bound_checks += 1
    verification["random_feasible_imputation_checks_inside_identification_bounds"] = bound_checks

    base_scores = ["pp_historical_product", "pp_development_reference_product", "pp_historical_max",
                   "pp_historical_mean", "pp_logistic_raw_product", "pp_knn_raw_product",
                   "bank_phage_prevalence_mean", "bank_phage_prevalence_product",
                   "training_phylogroup_susceptibility"]
    discrimination_rows = []
    for exclude_collision in [False, True]:
        keep = ~selected.identifier_ambiguous if exclude_collision else np.ones(100, bool)
        suffix = "exclude_AN03" if exclude_collision else "historical_100"
        full = selected.loc[keep]
        discrimination_rows += summarize_cohort(full, base_scores, "selected_" + suffix)
        common = full.loc[full.source_constituent_scores_available == 3]
        discrimination_rows += summarize_cohort(common, base_scores + ["source_released_product"], "source_complete_" + suffix)
    discrimination = pd.DataFrame(discrimination_rows)
    discrimination.to_csv(OUT / "discrimination_and_paired_controls.csv", index=False)

    paired = pd.DataFrame({"bacteria": host_names, "host_group": groups,
                           "selected_lytic": selected.observed_lytic,
                           "baseline_lytic": (outcomes["b1+b2+b3"] > 0).astype(int),
                           "selected_mlc": outcomes["p1+p2+p3"], "baseline_mlc": outcomes["b1+b2+b3"],
                           "selected_composition_key": selected.composition_key,
                           "baseline_composition_key": "|".join(BASELINE_BANK)})
    assert (paired.selected_composition_key != paired.baseline_composition_key).all()
    score_pairs = {}
    for name, values in [("pp_historical_product", historical), ("pp_development_reference_product", repaired),
                         ("pp_logistic_raw_product", raw_log), ("pp_knn_raw_product", raw_knn),
                         ("bank_phage_prevalence_product", prevalence)]:
        sel_name, base_name = name + "_selected", name + "_baseline"
        paired[sel_name] = product_score(take(values, selected_indices))
        paired[base_name] = product_score(take(values, baseline_indices))
        score_pairs[name] = (sel_name, base_name)
    paired["training_phylogroup_susceptibility_selected"] = phylo
    paired["training_phylogroup_susceptibility_baseline"] = phylo
    score_pairs["training_phylogroup_susceptibility"] = (
        "training_phylogroup_susceptibility_selected", "training_phylogroup_susceptibility_baseline")
    paired.to_csv(OUT / "paired_selected_baseline_triplets.csv", index=False)
    preferences, preference_hosts = [], []
    for excluded in [False, True]:
        subset = paired.loc[~paired.bacteria.isin(COLLISIONS)] if excluded else paired
        label = "exclude_AN03_99" if excluded else "historical_100"
        rows, host_rows = preference_summary(subset, score_pairs, label)
        preferences += rows
        preference_hosts += host_rows
    preference_table = pd.DataFrame(preferences)
    preference_table.to_csv(OUT / "paired_preference_and_policy_summary.csv", index=False)
    pd.DataFrame(preference_hosts).to_csv(OUT / "paired_preference_host_values.csv", index=False)

    measured = pd.concat([
        pd.DataFrame({"bacteria": host_names, "kind": "externally_selected",
                      "composition_key": selected.composition_key, "observed_lytic": selected.observed_lytic}),
        pd.DataFrame({"bacteria": host_names, "kind": "fixed_baseline", "composition_key": "|".join(BASELINE_BANK),
                      "observed_lytic": paired.baseline_lytic}),
    ], ignore_index=True)
    measured = measured.drop_duplicates(["bacteria", "composition_key"])
    assert len(measured) == 200
    composition_counts = measured.groupby(["kind", "composition_key"]).agg(
        hosts=("bacteria", "size"), positive=("observed_lytic", "sum")).reset_index()
    composition_counts["prevalence"] = composition_counts.positive / composition_counts.hosts
    composition_counts.to_csv(OUT / "unique_composition_counts.csv", index=False)

    # These small-cohort summaries support review, not model selection.
    limitations = [
        "Every analysis is retrospective; external outcomes were known before these controls were defined.",
        "No predictor was refitted. None of these comparisons is a new independent frozen validation.",
        "The original source records 233 of 300 selected constituent decision scores; 67 rescue scores are missing, on 36 hosts.",
        "The 64 complete hosts are enriched for selections made before rescue. Complete-case source comparison may not represent all 100 hosts.",
        "Source y_pred_proba comes from several recommender stages. Its product aggregation is imposed here for comparison and is not a published end-to-end cocktail probability.",
        "Per-host lower/upper source-score columns bound missing constituent scores. Full-cohort identification bounds set unknown scores to outcome-dependent extremes solely to bound AUROC/AP. These bounds are not model predictions, confidence intervals, or substitutes for missing released scores.",
        "Historical PrecisionPhage percentiles depend on the external query cohort; the later development-reference repair is reported separately.",
        "Intervals resample entire historical host sequence groups, conditional on existing predictions and the fixed bank/composition menu. They do not generalize over new phages, banks, or selection procedures.",
        "AN03 labels refer to different development/external genome identifiers, while the published distance table has one AN03. The 99-host and 63-complete-host sensitivity excludes that ambiguity; it does not silently replace historical estimates.",
        "Only 17 hosts have discordant binary outcomes between the two measured triplets. Within-host binary preference accuracy is estimated on those hosts with ties worth one half.",
        "Repeated compositions across different hosts are legitimate distinct assays and are retained. Canonicalization removes ordering differences; only exact duplicate host/composition records are deduplicated.",
        "MLC > 0 is observed lytic activity; high-titre clearing does not establish productive phage replication, synergy, clinical benefit, or resistance suppression.",
        "Only the unordered fixed-baseline triplet is assigned. Published CLB_P2 maps to bank 55989_P2 by exact identity of the complete ENA OZ035781.1 sequence plus a 60-base terminal duplication of its prefix in the local assembly. No identity is inferred for b1/b2/b3 individual phages or two-phage baseline mixtures.",
        "AP-minus-prevalence intervals recompute the prevalence baseline in each identical bootstrap replicate. A finite-sample permutation AP reference can differ from prevalence.",
        "No permutation significance test, exchangeability claim, or confirmatory p-value is supplied.",
        "The ordinal MLC preference diagnostic treats the source scale as ordered, not equally spaced biological efficacy.",
    ]
    report = {
        "status": "completed retrospective saved-score analysis",
        "utc_created": datetime.now(timezone.utc).isoformat(),
        "seed": SEED, "bootstrap_replicates": REPLICATES,
        "source_commit": SOURCE_COMMIT,
        "source_citation": "Gaborieau et al. (2024), Nature Microbiology 9:2847-2861",
        "source_doi": SOURCE_DOI, "source_accepted_manuscript": SOURCE_MANUSCRIPT,
        "baseline_identity_evidence": {
            "unordered_phages": list(BASELINE), "source": "Assembling phages in cocktails, source accepted manuscript Methods",
            "bank_identifiers": list(BASELINE_BANK),
            "exact_genome_alias_verification": baseline_alias_evidence,
            "verification": "Primary manuscript text retrieved in indexed search; it names LF82_P8, CLB_P2, LF73_P1 and confirms triplet, pair and single assays. Direct HAL download was access-controlled; no bypass attempted.",
            "baseline_assay_column": "b1+b2+b3",
            "individual_b_columns_assigned": False,
        },
        "aliases": ALIASES,
        "phage_aliases": PHAGE_ALIASES,
        "source_score_coverage": {
            "constituent_rows": 300, "available": int(source.y_pred_proba.notna().sum()),
            "missing": int(source.y_pred_proba.isna().sum()),
            "complete_hosts": int(available.sum()),
            "complete_host_positives": int(selected.loc[available, "observed_lytic"].sum()),
            "complete_host_negatives": int((1-selected.loc[available, "observed_lytic"]).sum()),
            "missing_scores_by_stage": source.loc[source.y_pred_proba.isna(), "chosen_at_step"].value_counts().to_dict(),
            "all_stage_counts": source.chosen_at_step.value_counts().to_dict(),
        },
        "source_full_cohort_identification_bounds": bound_rows,
        "assay_composition_counts": {"selected_hosts": 100, "selected_distinct_compositions": 18,
                                     "baseline_distinct_compositions": 1,
                                     "host_composition_records_after_deduplication": len(measured),
                                     "distinct_compositions_overall": int(measured.composition_key.nunique()),
                                     "exact_duplicate_host_composition_rows_removed": 0},
        "historical_groups": {"groups": 81, "rule": "connected components of published core-genome distances <= 1e-4; missing records singleton",
                              "reconstructed_from_frozen_source": True, "historical_group_count_verified": True},
        "verification": verification,
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "scipy": scipy.__version__, "scikit_learn": sklearn.__version__},
        "sources": {key: {"path": str(path), "sha256": hashes_before[key]} for key, path in source_paths.items()},
        "limitations": limitations,
    }
    for key, path in source_paths.items():
        if sha(path) != hashes_before[key]:
            raise RuntimeError(f"Source changed during this analysis: {key}")
    report["verification"]["all_source_hashes_unchanged_after_analysis"] = True
    report["verification"]["historical_and_repaired_scores_match_saved_artifacts"] = True
    (OUT / "analysis_manifest.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    common_rows = []
    for row in discrimination.loc[discrimination.cohort.str.startswith("source_complete")].itertuples():
        if row.score not in ["pp_historical_product", "pp_development_reference_product", "source_released_product"]:
            continue
        common_rows.append({"cohort": row.cohort, "score": row.score, "hosts": row.hosts,
                            "AUROC": f"{row.auroc:.3f} [{row.auroc_ci_low:.3f}, {row.auroc_ci_high:.3f}]",
                            "AP": f"{row.average_precision:.3f}", "prevalence": f"{row.prevalence:.3f}",
                            "AP minus prevalence": f"{row.ap_minus_prevalence:.3f} [{row.ap_minus_prevalence_ci_low:.3f}, {row.ap_minus_prevalence_ci_high:.3f}]"})
    pref_rows = []
    for row in preference_table.itertuples():
        if row.score not in ["pp_historical_product", "pp_development_reference_product", "bank_phage_prevalence_product"]:
            continue
        pref_rows.append({"cohort": row.cohort, "score": row.score, "discordant": row.binary_discordant_hosts,
                          "preference accuracy": f"{row.binary_discordant_preference_accuracy:.3f} [{row.binary_discordant_preference_accuracy_ci_low:.3f}, {row.binary_discordant_preference_accuracy_ci_high:.3f}]",
                          "chosen activity": f"{row.policy_observed_lytic_fraction:.3f}",
                          "change vs baseline": f"{row.policy_minus_fixed_baseline:.3f} [{row.policy_minus_fixed_baseline_ci_low:.3f}, {row.policy_minus_fixed_baseline_ci_high:.3f}]"})
    text = "\n\n".join([
        "# Retrospective cocktail controls for the Frontiers revision",
        "This analysis reuses saved scores and measured outcomes. No predictor was refitted. Results and intervals condition on the fixed phage bank and the externally chosen composition menu.",
        "## Released original-source scores on exactly common hosts",
        "The source file contains 233 numeric decision scores and 67 missing rescue scores. Exactly 64 hosts (56 lytic and 8 non-lytic selected cocktails) have three numeric original scores; excluding ambiguous AN03 leaves 63. The same product formula is applied to the source decision scores, without pretending that the source publication reported a cocktail probability. All methods below use identical hosts.",
        markdown_table(common_rows, ["cohort", "score", "hosts", "AUROC", "AP", "prevalence", "AP minus prevalence"]),
        "On the 64 complete hosts, the historical score's paired AUROC difference from the released-source aggregate is 0.258 (95% host-group interval -0.042 to 0.563); the AP difference is 0.074 (interval -0.001 to 0.179). Both intervals include zero. The matched point estimates therefore do not establish superiority, and missing rescue scores prevent an equivalent point comparison over the whole cohort.",
        "All-cohort identification bounds are also supplied: unknown rescue scores are allowed anywhere in [0,1]. The worst assignment takes positive cases to their lowest feasible cocktail score and negative cases to their highest; the best assignment reverses those choices. Because these assignments use known outcomes, they are bounds on an incompletely observed source comparator, not predicted scores or confidence intervals.",
        markdown_table([{**r, "AUROC": f"{r['auroc']:.3f}", "AP": f"{r['average_precision']:.3f}"} for r in bound_rows], ["cohort", "case", "hosts", "AUROC", "AP"]),
        "## Choice between two experimentally measured triplets for a host",
        "The source measured both the selected triplet and the fixed baseline LF82_P8/CLB_P2/LF73_P1 for all 100 hosts. Selected outcomes are 82/100 lytic and baseline outcomes 81/100: 73 both lytic, 10 neither, 9 selected-only and 8 baseline-only. Preference accuracy uses only the 17 binary-discordant hosts; score ties count one half. Chosen activity uses all hosts, with uniform random choice for score ties. These are retrospective comparisons within this two-option menu.",
        markdown_table(pref_rows, ["cohort", "score", "discordant", "preference accuracy", "chosen activity", "change vs baseline"]),
        "There are 18 unique selected triplets and one unique baseline triplet. Canonical host/composition keys yield 200 unique measured triplet records. Repeated compositions on different hosts are retained; their outcomes are different assays. Single and double mixture outcomes are deliberately unused because the individual baseline b1/b2/b3 order has not been independently established.",
        "## Interpretation and reproducibility",
        "Host-group percentile intervals use 2,000 resamples. The AP-minus-prevalence reference is recomputed within each resample. Both historical test-batch percentiles and the later development-reference repair are reported, rather than replacing the historical result. The complete-source subset is selection-biased toward non-rescue recommendations; it cannot support an unqualified whole-cohort superiority claim.",
        f"Source: [Gaborieau et al. (2024)]({SOURCE_DOI}), Methods ‘Assembling phages in cocktails’ and [public code](https://github.com/mdmparis/coli_phage_interactions_2023/tree/{SOURCE_COMMIT}). Endpoint: observed lytic activity (MLC > 0), which is not proof of productive replication at every titre.",
        "\n".join("- " + value for value in limitations),
        "Run with: `python scripts/frontiers_revision/cocktail_source_comparison.py`. Exact source/script checksums, aliases, software versions, reconstructed historical groups, matched scores, all controls, and verification results are supplied alongside this report.",
    ])
    (OUT / "ANALYSIS_SUMMARY.md").write_text(text + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT), "verification": verification,
                      "common_source_comparison": common_rows,
                      "paired_preference": pref_rows}, indent=2), flush=True)


if __name__ == "__main__":
    main()
