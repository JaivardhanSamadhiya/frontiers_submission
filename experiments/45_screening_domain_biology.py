"""Retrospective, fixed-definition diagnostics of saved OOF/external scores.

Does not access BASEL-ECOR outcomes and does not train or select a model.
See NEXT_STAGE_VALIDATION_PROTOCOL_2026_10_05.md.
"""
from pathlib import Path
import hashlib
import json
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.benchmarks.screening import (
    screening_curve, random_curve, cluster_mean_ci, distance_bin,
)
OUT = ROOT/'data/results_v2/screening_domain_biology'
BANK = ROOT/'data/results_v2/additional_banks'
SOURCE = ROOT/'external/gaborieau_coli_phage_interactions_2023'


def discrimination(y, score):
    y, score = np.asarray(y), np.asarray(score)
    return {'pairs': len(y), 'positive': int(y.sum()),
            'negative': int((y == 0).sum()),
            'auroc': float(roc_auc_score(y, score)) if len(np.unique(y)) == 2 else None,
            'average_precision': float(average_precision_score(y, score)) if y.sum() else None}


def screening(pairs, hosts):
    rows = []
    for (bank, grouping, host), part in pairs.groupby(['bank', 'grouping', 'host'], sort=True):
        meta = hosts.loc[(bank, grouping, host)]
        actual = screening_curve(part.label, part.score, part.phage)
        breadth = screening_curve(part.label, part.breadth_score, part.phage)
        random = random_curve(actual['candidates'], actual['positive'])
        for method, curve in [('model', actual), ('training_prevalence', breadth), ('random_expected', random)]:
            common = {'bank': bank, 'grouping': grouping, 'host': host,
                      'group': meta['group'], 'fold': meta['fold'], 'method': method,
                      'nearest_training_distance': meta.nearest_training_distance,
                      'positive_phages': actual['positive'], 'bank_size': actual['candidates'],
                      'first_active_rank_actionable_only': curve['first_hit'],
                      'normalized_hit_area': curve['normalized_hit_area']}
            for k in [1, 3, 5, 10]:
                if k > actual['candidates']: continue
                rows.append({**common, 'budget': k, 'hit': float(curve['hit'][k-1]),
                             'precision': float(curve['precision'][k-1]),
                             'active_recall_actionable_only': float(curve['recall'][k-1])})
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT/'screening_host_metrics.csv', index=False)
    summaries = []
    for (bank, grouping, budget), part in frame.groupby(['bank', 'grouping', 'budget'], sort=True):
        for metric in ['hit', 'precision', 'active_recall_actionable_only',
                       'first_active_rank_actionable_only', 'normalized_hit_area']:
            wide = part.pivot(index=['host', 'group'], columns='method', values=metric)
            groups = wide.index.get_level_values('group').to_numpy()
            for method in wide.columns:
                summaries.append({'bank': bank, 'grouping': grouping, 'budget': budget,
                                  'metric': metric, 'comparison': method,
                                  **cluster_mean_ci(wide[method], groups)})
            for baseline in ['training_prevalence', 'random_expected']:
                summaries.append({'bank': bank, 'grouping': grouping, 'budget': budget,
                                  'metric': metric, 'comparison': f'model_minus_{baseline}',
                                  **cluster_mean_ci(wide['model']-wide[baseline], groups)})
    pd.DataFrame(summaries).to_csv(OUT/'screening_summary.csv', index=False)
    return frame


def domain(pairs, hosts, screening_hosts):
    rows, associations, selective = [], [], []
    for (bank, grouping), part in pairs.groupby(['bank', 'grouping'], sort=True):
        meta = hosts.loc[(bank, grouping)].copy()
        q = part.merge(meta[['group', 'nearest_training_distance']], left_on='host', right_index=True,
                       how='left', validate='many_to_one')
        if q['group'].isna().any(): raise ValueError('Unresolved host groups')
        q['distance_bin'] = q.nearest_training_distance.map(distance_bin)
        s = screening_hosts[(screening_hosts.bank == bank) & (screening_hosts.grouping == grouping)
                            & (screening_hosts.budget == 3) & (screening_hosts.method == 'model')].copy()
        s['distance_bin'] = s.nearest_training_distance.map(distance_bin)
        for label in ['<=0.0001', '(0.0001,0.001]', '(0.001,0.01]', '>0.01', 'missing']:
            z = q[q.distance_bin == label]; h = s[s.distance_bin == label]
            if not len(z):
                rows.append({'bank': bank, 'grouping': grouping, 'distance_bin': label, 'hosts': 0, 'pairs': 0})
                continue
            stats = discrimination(z.label, z.score)
            coverage = cluster_mean_ci(h.hit, h.group)
            rows.append({'bank': bank, 'grouping': grouping, 'distance_bin': label,
                         'hosts': z.host.nunique(), **stats,
                         **{f'budget3_hit_{k}': v for k, v in coverage.items()}})
        # Within-source distances only; Spearman rho is descriptive, no rowwise p-value.
        error = q.assign(squared_score_error=(q.score-q.label)**2).groupby('host').squared_score_error.mean()
        context = meta.join(error).dropna(subset=['nearest_training_distance', 'squared_score_error'])
        rho = (float(spearmanr(context.nearest_training_distance, context.squared_score_error).statistic)
               if context.nearest_training_distance.nunique() > 1 and context.squared_score_error.nunique() > 1 else None)
        associations.append({'bank': bank, 'grouping': grouping, 'hosts': len(context),
                             'spearman_distance_vs_squared_score_error': rho})
        # Fixed domain cutoffs: no best cutoff is chosen from these results.
        for cutoff in [1e-4, 1e-3, 1e-2, float('inf')]:
            h = s[s.nearest_training_distance <= cutoff]
            selective.append({'bank': bank, 'grouping': grouping,
                              'maximum_distance': 'all_finite' if np.isinf(cutoff) else cutoff,
                              'host_fraction_retained': len(h)/len(s),
                              **cluster_mean_ci(h.hit, h.group)})
    pd.DataFrame(rows).to_csv(OUT/'distance_strata.csv', index=False)
    pd.DataFrame(associations).to_csv(OUT/'distance_error_associations.csv', index=False)
    pd.DataFrame(selective).to_csv(OUT/'distance_selective_screening.csv', index=False)


def biology():
    pair_path = ROOT/'data/results_v2/external_cocktail_validation/heldout_selected_pair_predictions_and_outcomes.csv'
    pair = pd.read_csv(pair_path)
    metadata = pd.read_csv(SOURCE/'dev/cocktails/data/test_collection.csv', sep=';', dtype=str)
    phage_metadata = pd.read_csv(SOURCE/'data/genomics/phages/guelin_collection.csv', sep=';', dtype=str)
    if pair.duplicated(['bacteria', 'phage']).any(): raise ValueError('Duplicate measured pair')
    # Previously documented ambiguous cross-collection identity; preserve raw artifact.
    excluded = pair[pair.bacteria == 'AN03'].copy()
    excluded.to_csv(OUT/'biology_excluded_ambiguous_AN03.csv', index=False)
    pair = pair[pair.bacteria != 'AN03'].copy()
    distances = pd.read_csv(ROOT/'data/results_v2/post_rejection_audit/external_host_distance_diagnostic.csv')
    pair = pair.merge(distances[['bacteria', 'nearest_development_host_distance']],
                      on='bacteria', how='left', validate='many_to_one')
    pair = pair.merge(metadata, on='bacteria', how='left', validate='many_to_one', indicator=True)
    if (pair._merge != 'both').any(): raise ValueError('Unresolved host metadata')
    pair.drop(columns='_merge', inplace=True)
    pair = pair.merge(phage_metadata[['phage', 'Genus', 'Phage_host_phylo']],
                      on='phage', how='left', validate='many_to_one', indicator=True)
    if (pair._merge != 'both').any(): raise ValueError('Unresolved phage metadata')
    pair.drop(columns='_merge', inplace=True)
    comparable = pair.phage_genus.notna() & pair.Genus.notna()
    if not (pair.loc[comparable, 'phage_genus'] == pair.loc[comparable, 'Genus']).all():
        raise ValueError('Phage genus disagreement')
    pair['predicted_at_0_5'] = (pair.score_ensemble >= .5).astype(int)
    for name, truth, pred in [('tp', 1, 1), ('tn', 0, 0), ('fp', 0, 1), ('fn', 1, 0)]:
        pair[name] = ((pair.productive == truth) & (pair.predicted_at_0_5 == pred)).astype(int)
    pair['wrong_at_0_5'] = (pair.productive != pair.predicted_at_0_5).astype(float)
    pair.to_csv(OUT/'biology_pair_annotations.csv', index=False)
    summaries = []
    for field in ['Clermont_Phylo', 'ST_Warwick', 'LPS_type', 'O-type', 'H-type', 'Genus', 'Phage_host_phylo']:
        for value, z in pair.groupby(pair[field].fillna('MISSING'), sort=True):
            ci = cluster_mean_ci(z.wrong_at_0_5, z.bacteria)
            stats = discrimination(z.productive, z.score_ensemble)
            tp, tn, fp, fn = [int(z[c].sum()) for c in ['tp', 'tn', 'fp', 'fn']]
            summaries.append({'field': field, 'value': value, 'hosts': z.bacteria.nunique(),
                              'phages': z.phage.nunique(), 'small_group_under_10_hosts': z.bacteria.nunique() < 10,
                              'positive_fraction': float(z.productive.mean()),
                              'mean_training_phage_prevalence': float(z.training_prevalence.mean()),
                              'hosts_missing_distance': z.loc[z.nearest_development_host_distance.isna(), 'bacteria'].nunique(),
                              'median_host_distance': z.drop_duplicates('bacteria').nearest_development_host_distance.median(),
                              **stats, 'tp': tp, 'tn': tn, 'fp': fp, 'fn': fn,
                              'fpr': fp/(fp+tn) if fp+tn else None,
                              'fnr': fn/(fn+tp) if fn+tp else None,
                              **{f'error_{k}': v for k, v in ci.items()}})
    pd.DataFrame(summaries).to_csv(OUT/'biology_subgroup_errors.csv', index=False)
    return {'pairs': len(pair), 'hosts': pair.bacteria.nunique(), 'phages': pair.phage.nunique(),
            'excluded_ambiguous_pairs': len(excluded),
            **{c: int(pair[c].sum()) for c in ['tp', 'tn', 'fp', 'fn']}}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = pd.read_csv(BANK/'pair_predictions.csv', dtype={'grouping': str})
    hosts = pd.read_csv(BANK/'host_selections.csv', dtype={'grouping': str})
    if hosts.duplicated(['bank', 'grouping', 'host']).any(): raise ValueError('Duplicate host selections')
    if pairs.duplicated(['bank', 'grouping', 'host', 'phage']).any(): raise ValueError('Duplicate OOF pairs')
    if not np.isin(pairs.label, [0, 1]).all(): raise ValueError('Incomplete assay matrix')
    for key, z in pairs.groupby(['bank', 'grouping']):
        if len(z) != z.host.nunique()*z.phage.nunique(): raise ValueError(f'Incomplete bank {key}')
    hosts = hosts.set_index(['bank', 'grouping', 'host']).sort_index()
    frame = screening(pairs, hosts)
    domain(pairs, hosts, frame)
    bio = biology()
    inputs = [BANK/'pair_predictions.csv', BANK/'host_selections.csv',
              ROOT/'data/results_v2/external_cocktail_validation/heldout_selected_pair_predictions_and_outcomes.csv',
              ROOT/'data/results_v2/post_rejection_audit/external_host_distance_diagnostic.csv',
              SOURCE/'dev/cocktails/data/test_collection.csv', SOURCE/'data/genomics/phages/guelin_collection.csv',
              ROOT/'NEXT_STAGE_VALIDATION_PROTOCOL_2026_10_05.md', Path(__file__),
              ROOT/'src/precisionphage/benchmarks/screening.py']
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(), 'biology': bio,
                'sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                'limitations': [
                    'Exploratory retrospective diagnostics; no model refitting or external outcome tuning.',
                    'Additional banks were trained separately. This is not frozen cross-bank transfer.',
                    'Hit means at least one active individual phage; it does not measure cocktail efficacy.',
                    'Bootstrap conditions on saved predictions and does not include training uncertainty.',
                    'Bank uncertainty resamples held-out host groups; biology resamples hosts, not independent pairs.',
                    'Shared phages and genetically related external hosts may leave additional dependence.',
                    'Distance bins are source-specific patristic units, not universally equivalent genomic distances.',
                    'Single-phage scores are not calibrated probabilities; 0.5 is a descriptive fixed threshold.',
                    'Biology pairs were selected by external experimenters; subgroup differences do not establish mechanisms.',
                    'No correction or significance claim across exploratory subgroups; small groups are retained and flagged.',
                ]}
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'status': 'complete', 'bank_groupings': pairs.groupby(['bank', 'grouping']).ngroups,
                      'biology': bio, 'output': str(OUT)}, indent=2))


if __name__ == '__main__':
    main()
