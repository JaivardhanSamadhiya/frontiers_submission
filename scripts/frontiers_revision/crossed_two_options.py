"""Post-review crossed-cluster sensitivity of saved two-option assay records.

No biological model is trained and no unmeasured assay is generated. Whole
hosts and selected-menu groups are reweighted, retaining paired observed
selected/fixed-baseline outcomes and conditioning on the single baseline menu.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'data/results_v2/frontiers_revision_2026_10_14/cocktail_controls/paired_selected_baseline_triplets.csv'
HISTORICAL = ROOT / 'data/results_v2/frontiers_revision_2026_10_14/cocktail_controls/paired_preference_and_policy_summary.csv'
OUT = ROOT / 'data/results_v2/frontiers_revision_2026_10_15/crossed_two_options'
SEED = 2026101505
REPLICATES = 2000
METHODS = ['pp_historical_product', 'pp_development_reference_product',
           'pp_logistic_raw_product', 'pp_knn_raw_product',
           'bank_phage_prevalence_product', 'training_phylogroup_susceptibility']
METRICS = ['binary_discordant_preference_accuracy', 'policy_observed_lytic_fraction',
           'always_selected_lytic_fraction', 'fixed_baseline_lytic_fraction',
           'uniform_random_lytic_fraction', 'oracle_lytic_fraction',
           'policy_minus_always_selected', 'policy_minus_fixed_baseline',
           'policy_minus_uniform_random', 'policy_minus_oracle',
           'ordinal_discordant_preference_accuracy']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def make_plan():
    OUT.mkdir(parents=True, exist_ok=True)
    spec = {'analysis': 'Retrospective amendment: two-measured-options crossed dependence sensitivity',
        'specified_utc_before_computation': datetime.now(timezone.utc).isoformat(),
        'status': 'Post hoc; original external outcomes and earlier analyses already inspected',
        'script_sha256': sha(__file__), 'input_sha256': sha(SOURCE),
        'historical_summary_sha256': sha(HISTORICAL), 'base_seed': SEED,
        'cohorts': ['historical_100', 'exclude_AN03_99'], 'valid_replicates_per_cohort': REPLICATES,
        'algorithm': 'Independently sample H host-cluster labels H times and C selected-composition labels C times with replacement. Original host record weight is product of its two multiplicities. Preserve both measured outcomes and scores within each host. Baseline composition is fixed, identical for every host and not sampled as an identifiable population. Retain first 2000 draws with positive total weight and positive binary/ordinal-discordant weight. No outcome balance or result-based selection. All methods/metrics/contrasts share the same retained weights. Report attempts and exclusions. Percentile intervals at 2.5/97.5.',
        'tie_rule': 'Choice weight is 1 if selected score exceeds baseline, 0 if below and 0.5 if exactly equal. Chosen binary activity for a tie is the actual mean of the two observed outcomes, not unconditional activity 0.5. Discordant binary/ordinal preference ties receive half credit.',
        'metrics': METRICS, 'methods': METHODS,
        'interpretation': 'Sensitivity to repeated selected menus, conditional on the existing fixed baseline menu, fixed predictions and observed paired assay records. Not a new frozen validation, prospective laboratory policy, selection-mechanism test or composition-population confidence guarantee. Only 17 hosts have discordant binary outcomes. Sparse crossed cells, 18 non-random selected menus and phages shared across distinct menus limit inference. No confirmatory p-value, multiplicity-adjusted superiority claim or therapeutic interpretation.',
    }
    path = OUT/'amendment_protocol.json'
    if path.exists():
        old = json.loads(path.read_text())
        for key in ['script_sha256', 'input_sha256', 'historical_summary_sha256']:
            if old[key] != spec[key]:
                raise RuntimeError('Existing dated amendment does not match source/code; do not overwrite')
        return old
    dump(path, spec)
    return spec


def sample_weights(hosts, compositions, binary_discordant, ordinal_discordant,
                   replicates=REPLICATES, seed=SEED):
    h, hi = np.unique(hosts, return_inverse=True)
    c, ci = np.unique(compositions, return_inverse=True)
    binary_discordant, ordinal_discordant = np.asarray(binary_discordant, bool), np.asarray(ordinal_discordant, bool)
    if len(hi) != len(ci) or len(hi) != len(binary_discordant) or len(hi) != len(ordinal_discordant):
        raise ValueError('Grouping shape differs')
    if not binary_discordant.any() or not ordinal_discordant.any():
        raise ValueError('No discordant records')
    rng = np.random.default_rng(seed)
    weights, hcounts, ccounts = [], [], []
    reject = {'empty': 0, 'no_binary_discordance': 0, 'no_ordinal_discordance': 0}
    attempts = 0
    while len(weights) < replicates and attempts < 100*replicates:
        attempts += 1
        hc = np.bincount(rng.integers(len(h), size=len(h)), minlength=len(h))
        cc = np.bincount(rng.integers(len(c), size=len(c)), minlength=len(c))
        w = hc[hi]*cc[ci]
        if w.sum() == 0:
            reject['empty'] += 1
            continue
        if w[binary_discordant].sum() == 0:
            reject['no_binary_discordance'] += 1
            continue
        if w[ordinal_discordant].sum() == 0:
            reject['no_ordinal_discordance'] += 1
            continue
        weights.append(w)
        hcounts.append(hc)
        ccounts.append(cc)
    if len(weights) != replicates:
        raise RuntimeError('Insufficient valid draws')
    return (np.asarray(weights), np.asarray(hcounts), np.asarray(ccounts), h, c, hi, ci,
            {'attempts': attempts, 'valid': len(weights), 'rejections': reject, 'seed': seed})


def policy_values(selected_y, baseline_y, selected_mlc, baseline_mlc, selected_score, baseline_score, weights):
    sy, by = np.asarray(selected_y, float), np.asarray(baseline_y, float)
    sm, bm = np.asarray(selected_mlc, float), np.asarray(baseline_mlc, float)
    ss, bs = np.asarray(selected_score, float), np.asarray(baseline_score, float)
    w = np.atleast_2d(np.asarray(weights, float))
    if any(arr.shape != sy.shape for arr in [by, sm, bm, ss, bs]) or w.shape[1] != len(sy):
        raise ValueError('Shape mismatch')
    if any(not np.isfinite(a).all() for a in [sy, by, sm, bm, ss, bs, w]):
        raise ValueError('Nonfinite input')
    if (w < 0).any() or not np.isin(sy, [0, 1]).all() or not np.isin(by, [0, 1]).all():
        raise ValueError('Invalid labels or weights')
    choice = np.where(ss > bs, 1., np.where(ss < bs, 0., .5))
    dy, dm = sy != by, sm != bm
    total, binary_n, ordinal_n = w.sum(axis=1), w[:, dy].sum(axis=1), w[:, dm].sum(axis=1)
    if (total <= 0).any() or (binary_n <= 0).any() or (ordinal_n <= 0).any():
        raise ValueError('Positive total and discordant weights required')
    observed = choice*sy+(1-choice)*by
    correct = np.where(sy > by, choice, 1-choice)
    ordinal_correct = np.where(sm > bm, choice, 1-choice)
    def mean(values):
        return w @ values / total
    activity = mean(observed)
    selected, baseline = mean(sy), mean(by)
    random, oracle = mean((sy+by)/2), mean(np.maximum(sy, by))
    values = np.column_stack([
        w[:, dy] @ correct[dy] / binary_n, activity, selected, baseline,
        random, oracle, activity-selected, activity-baseline, activity-random,
        activity-oracle, w[:, dm] @ ordinal_correct[dm] / ordinal_n])
    return values, choice, observed


def analyze_cohort(frame, cohort, seed):
    sy, by = frame.selected_lytic.to_numpy(), frame.baseline_lytic.to_numpy()
    sm, bm = frame.selected_mlc.to_numpy(), frame.baseline_mlc.to_numpy()
    binary_d, ordinal_d = sy != by, sm != bm
    w, hc, cc, h, c, hi, ci, meta = sample_weights(
        frame.host_group, frame.selected_composition_key, binary_d, ordinal_d, seed=seed)
    assert np.array_equal(w, hc[:, hi]*cc[:, ci])
    points, boots, rows, per_host, differences = {}, {}, [], [], []
    for method in METHODS:
        ss, bs = frame[method+'_selected'].to_numpy(), frame[method+'_baseline'].to_numpy()
        point, choice, observed = policy_values(sy, by, sm, bm, ss, bs, np.ones(len(frame)))
        boot, _, _ = policy_values(sy, by, sm, bm, ss, bs, w)
        points[method], boots[method] = point[0], boot
        low, high = np.quantile(boot, [.025, .975], axis=0)
        row = {'cohort': cohort, 'method': method, 'hosts': len(frame),
            'host_groups': len(h), 'selected_compositions': len(c), 'baseline_compositions': 1,
            'binary_discordant_hosts': int(binary_d.sum()), 'ordinal_discordant_hosts': int(ordinal_d.sum()),
            'selected_only_active': int(((sy == 1) & (by == 0)).sum()),
            'baseline_only_active': int(((by == 1) & (sy == 0)).sum()),
            'both_active': int(((sy == 1) & (by == 1)).sum()), 'neither_active': int(((sy == 0) & (by == 0)).sum()),
            'score_ties': int((ss == bs).sum()), **{k: meta[k] for k in ['seed', 'attempts', 'valid']}}
        for j, metric in enumerate(METRICS):
            row.update({metric: point[0, j], metric+'_ci_low': low[j], metric+'_ci_high': high[j]})
        rows.append(row)
        for i in range(len(frame)):
            per_host.append({'cohort': cohort, 'method': method, 'bacteria': frame.bacteria.iloc[i],
                'host_group': frame.host_group.iloc[i], 'selected_composition': frame.selected_composition_key.iloc[i],
                'selected_score': ss[i], 'baseline_score': bs[i], 'selected_outcome': sy[i], 'baseline_outcome': by[i],
                'probability_choose_selected_with_tie_randomization': choice[i],
                'observed_activity_if_score_policy': observed[i], 'binary_discordant': binary_d[i]})
    for a, b in combinations(METHODS, 2):
        low, high = np.quantile(boots[a]-boots[b], [.025, .975], axis=0)
        row = {'cohort': cohort, 'method_a': a, 'method_b': b, 'definition': 'method_a minus method_b'}
        for j, metric in enumerate(METRICS):
            row.update({metric+'_difference': points[a][j]-points[b][j],
                        metric+'_difference_ci_low': low[j], metric+'_difference_ci_high': high[j]})
        differences.append(row)
    np.savez_compressed(OUT/f'bootstrap_{cohort}.npz', weights=w, host_counts=hc, composition_counts=cc,
        host_labels=h, composition_labels=c, host_inverse=hi, composition_inverse=ci,
        bacteria=frame.bacteria.to_numpy(str), methods=np.asarray(METHODS), metric_names=np.asarray(METRICS),
        metric_replicates=np.stack([boots[m] for m in METHODS]))
    return rows, differences, per_host, {'cohort': cohort, **meta}


def main():
    spec = make_plan()
    hashes_before = {str(p): sha(p) for p in [SOURCE, HISTORICAL]}
    frame = pd.read_csv(SOURCE)
    assert len(frame) == 100 and frame.bacteria.is_unique
    assert frame.baseline_composition_key.nunique() == 1 and frame.selected_composition_key.nunique() == 18
    assert frame.host_group.nunique() == 81
    assert int((frame.selected_lytic != frame.baseline_lytic).sum()) == 17
    rows, differences, host_rows, metas = [], [], [], []
    for offset, (name, part) in enumerate([('historical_100', frame), ('exclude_AN03_99', frame[frame.bacteria.ne('AN03')])]):
        r, d, h, m = analyze_cohort(part, name, SEED+offset)
        rows += r
        differences += d
        host_rows += h
        metas.append(m)
    result = pd.DataFrame(rows)
    old = pd.read_csv(HISTORICAL)
    verified = 0
    for row in result.itertuples():
        prev = old[old.cohort.eq(row.cohort) & old.score.eq(row.method)].iloc[0]
        for metric in METRICS:
            if metric in prev.index:
                np.testing.assert_allclose(getattr(row, metric), prev[metric], atol=1e-13, rtol=0)
                verified += 1
    result.to_csv(OUT/'crossed_policy_metrics.csv', index=False)
    pd.DataFrame(differences).to_csv(OUT/'crossed_policy_paired_differences.csv', index=False)
    pd.DataFrame(host_rows).to_csv(OUT/'policy_host_values.csv', index=False)
    summary = {'status': 'Completed retrospective amendment; historical output preserved',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'protocol_sha256': sha(OUT/'amendment_protocol.json'),
        'cohorts': metas, 'point_estimates_reproduced_from_historical_summary': verified,
        'software': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__},
        'input_hashes': hashes_before, 'source_files_unchanged': all(sha(p) == digest for p, digest in hashes_before.items()),
        'limits': spec['interpretation']}
    if not summary['source_files_unchanged']:
        raise RuntimeError('Input changed')
    lines = ['# Two-measured-options crossed dependence sensitivity', '',
        'This is a post hoc uncertainty amendment, not a new outcome-locked analysis. No saved source or historical interval was overwritten.', '',
        'Each replicate independently resamples whole host-cluster labels and whole selected-menu labels. Product multiplicities reweight paired measured outcomes for each host. The single fixed baseline menu is conditioned upon. All metrics and methods share 2,000 valid draws. Preference accuracy uses weighted discordant hosts only; policy activity and gains use all reweighted hosts. Tied choices use the actual mean of the two observed outcomes.', '',
        '| Cohort | Method | Preference on 17 discordant hosts [95% interval] | Chosen activity [95% interval] | Fixed-baseline activity [95% interval] | Change versus fixed baseline [95% interval] |',
        '| --- | --- | --- | --- | --- | --- |']
    for row in result.itertuples():
        parts = [row.cohort, row.method]
        for m in ['binary_discordant_preference_accuracy', 'policy_observed_lytic_fraction', 'fixed_baseline_lytic_fraction', 'policy_minus_fixed_baseline']:
            parts.append(f'{getattr(row,m):.3f} [{getattr(row,m+"_ci_low"):.3f}, {getattr(row,m+"_ci_high"):.3f}]')
        lines.append('| '+' | '.join(parts)+' |')
    lines += ['', '## Interpretation', '', spec['interpretation'], '',
        'The historical and repaired point estimates are unchanged. These wider dependence sensitivities do not establish improved selection over the fixed baseline. The fixed baseline is one composition used across every host; a baseline-composition population effect cannot be estimated from it. Selected-menu resampling does not capture phage-level dependence across distinct selected menus. Rejected zero-discordance draws mean intervals condition on an estimable preference denominator; counts are retained in the manifest.', '']
    (OUT/'REPORT.md').write_text('\n'.join(lines), encoding='utf-8')
    summary['output_hashes'] = {p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.name != 'manifest.json'}
    dump(OUT/'manifest.json', summary)
    print(result[result.method.isin(METHODS[:2])].to_string(index=False))
    print(json.dumps({'output': str(OUT), 'cohorts': metas, 'historical_points_verified': verified}, indent=2))


if __name__ == '__main__':
    main()
