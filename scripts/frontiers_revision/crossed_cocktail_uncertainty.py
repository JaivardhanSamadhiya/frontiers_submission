"""Retrospective crossed-cluster uncertainty for saved cocktail scores only.

This is a generic saved-score statistical audit. It does not fit a predictor,
resolve genomes, select a biological composition, or change historical outputs.
The dated analysis plan is written before any new bootstrap metrics are run.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
OLD = ROOT / 'data/results_v2/frontiers_revision_2026_10_14'
OUT = ROOT / 'data/results_v2/frontiers_revision_2026_10_15/crossed_cocktail'
SCORES = OLD / 'cocktail_controls/selected_triplet_scores.csv'
SEED = 2026101501
REPLICATES = 2000
METRICS = ['auroc', 'average_precision', 'prevalence', 'ap_minus_prevalence']
METHODS = [
    'pp_historical_product', 'pp_development_reference_product',
    'pp_historical_max', 'pp_historical_mean', 'pp_logistic_raw_product',
    'pp_knn_raw_product', 'bank_phage_prevalence_mean',
    'bank_phage_prevalence_product', 'training_phylogroup_susceptibility',
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def canonical_composition(value):
    parts = str(value).split('|')
    if len(parts) != 3 or len(set(parts)) != 3 or any(not p for p in parts):
        raise ValueError('Expected exactly three distinct saved composition identifiers')
    return '|'.join(sorted(parts))


def weighted_metrics_many(y, scores, weights):
    """Vectorized tied-score AUROC/AP, retaining original row-weighted estimand."""
    y, scores = np.asarray(y, int), np.asarray(scores, float)
    weights = np.atleast_2d(np.asarray(weights, float))
    if y.ndim != 1 or scores.shape != y.shape or weights.shape[1] != len(y):
        raise ValueError('Mismatched shapes')
    if not np.isin(y, [0, 1]).all() or not np.isfinite(scores).all():
        raise ValueError('Invalid outcome or score')
    if not np.isfinite(weights).all() or (weights < 0).any():
        raise ValueError('Invalid weights')
    order = np.argsort(-scores, kind='stable')
    starts = np.r_[0, np.flatnonzero(np.diff(scores[order])) + 1]
    pos = np.add.reduceat(weights[:, order] * y[order], starts, axis=1)
    neg = np.add.reduceat(weights[:, order] * (1-y[order]), starts, axis=1)
    p, n = pos.sum(axis=1), neg.sum(axis=1)
    if (p <= 0).any() or (n <= 0).any():
        raise ValueError('Both classes required')
    cp, cn = pos.cumsum(axis=1), neg.cumsum(axis=1)
    auc = np.sum(pos * (n[:, None] - cn + .5*neg), axis=1) / (p*n)
    precision = np.divide(cp, cp+cn, out=np.zeros_like(cp), where=cp+cn > 0)
    ap = (pos*precision).sum(axis=1) / p
    prevalence = p/(p+n)
    return np.column_stack([auc, ap, prevalence, ap-prevalence])


def crossed_weights(y, host_groups, composition_groups, replicates=REPLICATES, seed=SEED):
    """Pigeonhole-style bootstrap; product of independent sampled group counts.

    Draw H host-cluster labels from H labels and C composition labels from C
    labels independently with replacement. A row receives the product of its
    two multiplicities. All metrics/methods in a cohort share these row weights.
    Reject only draws lacking either outcome class. Never rebalance labels or
    invent an unobserved host/composition outcome.
    """
    y = np.asarray(y, int)
    hg, hi = np.unique(host_groups, return_inverse=True)
    cg, ci = np.unique(composition_groups, return_inverse=True)
    if len(y) != len(hi) or len(y) != len(ci) or not np.isin(y, [0, 1]).all():
        raise ValueError('Invalid row grouping or outcome')
    if np.unique(y).size != 2:
        raise ValueError('Both classes required')
    rng = np.random.default_rng(seed)
    weights, host_counts, composition_counts = [], [], []
    attempts = 0
    while len(weights) < replicates and attempts < 100*replicates:
        attempts += 1
        hc = np.bincount(rng.integers(len(hg), size=len(hg)), minlength=len(hg))
        cc = np.bincount(rng.integers(len(cg), size=len(cg)), minlength=len(cg))
        w = hc[hi]*cc[ci]
        if np.unique(y[w > 0]).size != 2:
            continue
        weights.append(w)
        host_counts.append(hc)
        composition_counts.append(cc)
    if len(weights) != replicates:
        raise RuntimeError('Insufficient valid draws')
    return (np.asarray(weights), np.asarray(host_counts), np.asarray(composition_counts),
            hg, cg, hi, ci, attempts)


def plan():
    """Specify dated protocol before calculation; never claim preregistration."""
    OUT.mkdir(parents=True, exist_ok=True)
    spec = {
        'analysis': 'Post-review crossed host-cluster and composition uncertainty sensitivity',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'status': 'Retrospective; external labels already inspected before this plan',
        'seed': SEED, 'valid_replicates_per_cohort': REPLICATES,
        'input_score_sha256': sha(SCORES), 'script_sha256': sha(__file__),
        'algorithm': 'Independently resample H historical host-cluster labels H times and C canonical composition labels C times with replacement; original observed row weight is product of its two counts. No Cartesian completion or unobserved assay creation. Reject draws with zero positive or zero negative total weight. Shared weights pair every method contrast. Percentile 2.5/97.5 intervals; AP minus reweighted prevalence within each draw.',
        'cohorts': ['historical_100', 'exclude_AN03_99', 'source_complete_64', 'source_complete_exclude_AN03_63'],
        'methods': METHODS + ['source_released_product on complete cases only'],
        'shuffle_audit': 'Reconstruct the original 20261005 RNG state after its 2000 valid host-only bootstrap draws; permute complete saved composition rows 1000 times. Verify each permutation preserves every composition multiplicity and reproduces the saved diagnostic metrics. Do not create a new formal permutation null or p-value.',
        'interpretation': 'Additional post hoc dependence sensitivity, not a replacement for or original prespecified historical CI. Crossed resampling may be conservative/unstable with sparse, non-randomly selected cells and 18 compositions. It does not establish generalization to new banks, phages, selection procedures or composition populations. Shared phages across distinct compositions remain outside this grouping model. No multiplicity-adjusted confirmatory superiority test.',
    }
    dest = OUT / 'analysis_plan.json'
    if dest.exists():
        old = json.loads(dest.read_text())
        if old['input_score_sha256'] != spec['input_score_sha256'] or old['script_sha256'] != spec['script_sha256']:
            raise RuntimeError('Existing plan differs; preserve it rather than overwrite dated specifications')
    else:
        write_json(dest, spec)
    return json.loads(dest.read_text())


def cohort_metrics(frame, name, seed):
    y = frame.observed_lytic.to_numpy(int)
    methods = METHODS + (['source_released_product'] if name.startswith('source_complete') else [])
    w, hc, cc, hg, cg, hi, ci, attempts = crossed_weights(
        y, frame.host_group.to_numpy(), frame.composition_key.to_numpy(), seed=seed)
    assert np.array_equal(w, hc[:, hi]*cc[:, ci])
    assert (hc.sum(axis=1) == len(hg)).all() and (cc.sum(axis=1) == len(cg)).all()
    point, boots, rows, differences = {}, {}, [], []
    for method in methods:
        s = frame[method].to_numpy(float)
        point[method] = weighted_metrics_many(y, s, np.ones(len(y)))[0]
        np.testing.assert_allclose(point[method][:2], [roc_auc_score(y, s), average_precision_score(y, s)], rtol=0, atol=1e-13)
        boots[method] = weighted_metrics_many(y, s, w)
        low, high = np.quantile(boots[method], [.025, .975], axis=0)
        row = {'cohort': name, 'method': method, 'records': len(frame), 'positive': int(y.sum()),
               'host_groups': len(hg), 'compositions': len(cg), 'valid_replicates': len(w), 'attempts': attempts, 'seed': seed}
        for j, metric in enumerate(METRICS):
            row.update({metric: point[method][j], metric+'_ci_low': low[j], metric+'_ci_high': high[j]})
        rows.append(row)
    for a, b in combinations(methods, 2):
        low, high = np.quantile(boots[a]-boots[b], [.025, .975], axis=0)
        row = {'cohort': name, 'method_a': a, 'method_b': b, 'definition': 'a minus b'}
        for j, metric in enumerate(METRICS):
            row.update({metric+'_difference': point[a][j]-point[b][j],
                        metric+'_difference_ci_low': low[j], metric+'_difference_ci_high': high[j]})
        differences.append(row)
    np.savez_compressed(OUT / f'bootstrap_{name}.npz', weights=w, host_counts=hc,
                        composition_counts=cc, host_labels=hg, composition_labels=cg,
                        host_inverse=hi, composition_inverse=ci,
                        bacteria=frame.bacteria.to_numpy(str), labels=y,
                        methods=np.array(methods), metric_names=np.array(METRICS),
                        metric_replicates=np.stack([boots[m] for m in methods]))
    return rows, differences


def audit_original_shuffle(frame):
    """Whole-menu permutation already preserves repetition; verify exact output."""
    frozen = ROOT / 'data/results_v2/external_cocktail_validation'
    old_shuffle = ROOT / 'data/results_v2/post_rejection_audit/composition_shuffle_diagnostic.csv'
    original_cocktails = pd.read_csv(frozen / 'heldout_cocktail_predictions_and_outcomes.csv')
    names = original_cocktails.bacteria.tolist()
    part = frame.set_index('bacteria').loc[names]
    y = part.observed_lytic.to_numpy(int)
    hg, hi = np.unique(part.host_group.to_numpy(), return_inverse=True)
    rng = np.random.default_rng(20261005)
    valid, attempts = 0, 0
    while valid < 2000:
        attempts += 1
        counts = np.bincount(rng.integers(len(hg), size=len(hg)), minlength=len(hg))
        if np.unique(y[counts[hi] > 0]).size == 2:
            valid += 1
    all_pairs = pd.read_csv(frozen / 'frozen_pair_predictions.csv')
    matrix = all_pairs.pivot(index='bacteria', columns='phage', values='score_ensemble').loc[names]
    index = {value: i for i, value in enumerate(matrix.columns)}
    compositions = np.array([[index[p] for p in s.split('|')] for s in original_cocktails.predicted_constituents])
    keys = np.asarray([canonical_composition(s) for s in original_cocktails.predicted_constituents])
    expected_counts = Counter(keys)
    expected = pd.read_csv(old_shuffle)
    metrics, permutations = [], []
    for _ in range(len(expected)):
        perm = rng.permutation(len(names))
        assert Counter(keys[perm]) == expected_counts
        score = 1-np.prod(1-np.take_along_axis(matrix.to_numpy(), compositions[perm], axis=1), axis=1)
        metrics.append(weighted_metrics_many(y, score, np.ones(len(y)))[0][:2])
        permutations.append(perm)
    actual = np.asarray(metrics)
    np.testing.assert_allclose(actual, expected[['auroc', 'average_precision']].to_numpy(), rtol=0, atol=1e-13)
    np.savez_compressed(OUT / 'original_shuffle_verified.npz', permutations=np.asarray(permutations),
                        bacteria=np.asarray(names), canonical_composition=keys,
                        diagnostic_metrics=actual, labels=y)
    return {
        'existing_draws': len(expected), 'distinct_compositions': len(expected_counts),
        'bootstrap_rng_draws_consumed_before_shuffle': attempts,
        'all_draws_preserve_whole_compositions_and_exact_counts': True,
        'all_saved_metrics_reproduced_at_absolute_tolerance': 1e-13,
        'composition_counts': dict(sorted(expected_counts.items())),
        'auroc_median': float(np.median(actual[:, 0])),
        'auroc_central_95_range': np.quantile(actual[:, 0], [.025, .975]).tolist(),
        'average_precision_median': float(np.median(actual[:, 1])),
        'average_precision_central_95_range': np.quantile(actual[:, 1], [.025, .975]).tolist(),
        'formal_permutation_p_value': None,
        'scope': 'Fixed host profiles, original observed labels, fixed bank and observed multiset of complete triplets; changes only which host receives each already observed triplet. Assay outcomes for reassigned compositions are not known; this is a score-association diagnostic, not simulated measured outcomes or a confirmatory exchangeability test.',
        'sources': {str(p.relative_to(ROOT)): sha(p) for p in [old_shuffle, frozen / 'frozen_pair_predictions.csv', frozen / 'heldout_cocktail_predictions_and_outcomes.csv', ROOT / 'experiments/34_post_rejection_evidence_audit.py']},
    }


def main():
    protocol = plan()
    before = sha(SCORES)
    frame = pd.read_csv(SCORES)
    assert len(frame) == 100 and frame.bacteria.is_unique
    assert int(frame.observed_lytic.sum()) == 82 and frame.host_group.nunique() == 81
    frame['composition_key'] = frame.composition_key.map(canonical_composition)
    assert frame.composition_key.nunique() == 18
    counts = frame.groupby('composition_key').agg(records=('bacteria', 'size'),
        positive=('observed_lytic', 'sum'), host_groups=('host_group', 'nunique')).reset_index()
    counts['negative'] = counts.records-counts.positive
    counts.to_csv(OUT / 'composition_groups.csv', index=False)
    cohorts = [('historical_100', frame), ('exclude_AN03_99', frame[frame.bacteria.ne('AN03')]),
               ('source_complete_64', frame[frame.source_constituent_scores_available.eq(3)]),
               ('source_complete_exclude_AN03_63', frame[frame.source_constituent_scores_available.eq(3) & frame.bacteria.ne('AN03')])]
    rows, deltas = [], []
    for offset, (name, part) in enumerate(cohorts):
        r, d = cohort_metrics(part, name, SEED+offset)
        rows += r
        deltas += d
    result, differences = pd.DataFrame(rows), pd.DataFrame(deltas)
    result.to_csv(OUT / 'crossed_metrics.csv', index=False)
    differences.to_csv(OUT / 'crossed_paired_differences.csv', index=False)
    shuffle = audit_original_shuffle(frame)
    write_json(OUT / 'original_shuffle_audit.json', shuffle)
    if sha(SCORES) != before:
        raise RuntimeError('Input changed')
    manifest = {'status': 'Completed retrospective saved-score dependence sensitivity',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'analysis_plan_sha256': sha(OUT/'analysis_plan.json'),
        'script_sha256': sha(__file__), 'input_score_sha256': before,
        'historical_files_modified': False, 'original_shuffle_verified': True,
        'metric_rows': len(rows), 'paired_contrast_rows': len(deltas),
        'versions': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__, 'sklearn': sklearn.__version__},
        'cautions': protocol['interpretation'],
        'outputs': {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'manifest.json'}}
    lines = ['# Retrospective crossed host-cluster × composition sensitivity', '',
        'The historical host-only intervals remain unchanged. This additional analysis was planned after outcomes were known and is not preregistration.', '',
        'Each draw independently resamples host-cluster and canonical-triplet labels; each original record receives the product of the two multiplicities. The same weights pair all methods. AP-minus-prevalence recomputes the prevalence in each draw. Exactly 2,000 two-class draws were retained per cohort.', '',
        '| Cohort | Score | N | Compositions | AUROC [95% interval] | AP [95% interval] | AP minus prevalence [95% interval] |',
        '| --- | --- | --- | --- | --- | --- | --- |']
    for row in result.itertuples():
        lines.append(f'| {row.cohort} | {row.method} | {row.records} | {row.compositions} | {row.auroc:.3f} [{row.auroc_ci_low:.3f}, {row.auroc_ci_high:.3f}] | {row.average_precision:.3f} [{row.average_precision_ci_low:.3f}, {row.average_precision_ci_high:.3f}] | {row.ap_minus_prevalence:.3f} [{row.ap_minus_prevalence_ci_low:.3f}, {row.ap_minus_prevalence_ci_high:.3f}] |')
    lines += ['', '## Audit of the existing composition shuffle', '',
        f'All {shuffle["existing_draws"]:,} saved draws were independently reproduced. Each permutes 100 **whole triplet rows**, preserving every one of the 18 composition counts; no constituent is shuffled independently. The AUROC diagnostic central 95% range remains {shuffle["auroc_central_95_range"][0]:.3f}–{shuffle["auroc_central_95_range"][1]:.3f}. This is not a formal permutation null or p-value.', '',
        '## Limits', '', protocol['interpretation'], '',
        'Only five of the 18 selected compositions have any negative outcomes, and three repeated compositions account for 72 of the 100 records. The non-random externally selected menu and sparse observed host/composition cells limit any composition-population interpretation. Distinct menus can share phages; grouping complete triplets does not eliminate all phage-level dependence. Host-only intervals are conditional on the saved menu, while these crossed intervals are an exploratory sensitivity to an additional grouping assumption. Neither is a guarantee of calibrated population coverage.', '',
        'Exact source inputs, bootstrap row weights, factor multiplicities, per-method replicates, paired contrasts, full permutation indices and checksums are retained. No historical result has been overwritten.', '']
    (OUT/'REPORT.md').write_text('\n'.join(lines), encoding='utf-8')
    manifest['outputs']['REPORT.md'] = sha(OUT/'REPORT.md')
    write_json(OUT/'manifest.json', manifest)
    print(result[result.method.isin(['pp_historical_product', 'pp_development_reference_product', 'source_released_product'])].to_string(index=False))
    print(json.dumps({'output': str(OUT), 'metric_rows': len(rows), 'contrasts': len(deltas), 'original_shuffle_verified': True}, indent=2))


if __name__ == '__main__':
    main()
