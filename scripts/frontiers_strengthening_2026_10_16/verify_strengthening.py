"""Independent checks of saved strengthening predictions and split constraints.

Does not fit a model or modify historical results. Incomplete work is reported
as incomplete, never accepted from a stale status file.
"""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import joblib
from sklearn.metrics import roc_auc_score, average_precision_score
from prediction_portability_policy import assert_saved_predictions

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/results_v2/frontiers_strengthening_2026_10_16'
OUT = ROOT / 'output/frontiers_strengthening_2026_10_17'


def check_metric(rows, reported):
    np.testing.assert_allclose(roc_auc_score(rows.label, rows.score), reported.auroc, atol=1e-10)
    np.testing.assert_allclose(average_precision_score(rows.label, rows.score), reported.average_precision, atol=1e-10)
    assert len(rows) == reported.n and int(rows.label.sum()) == reported.positive


def audit_models(directory, nested=False, portable_float32=False):
    complete = directory / 'completion.json'
    if not complete.exists():
        return {'complete': False, 'reason': 'No completion artifact'}
    state = json.loads(complete.read_text())
    assert state['complete']
    protocol = json.loads((directory / 'retrospective_protocol.json').read_text())
    assert hashlib.sha256((directory / 'retrospective_protocol.json').read_bytes()).hexdigest() == state['protocol_sha256']
    for name, digest in protocol['inputs'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    rows = pd.read_csv(directory / 'input_rows.csv')
    features = np.load(directory / 'input_features.npz')['X']
    # Preserve exact serialized double ties; default parser may perturb a tie.
    prediction = pd.read_csv(directory / 'predictions.csv',float_precision='round_trip')
    split_records = json.loads((directory / 'splits.json').read_text())
    checked_reload = checked_splits = 0
    for split in split_records:
        tr, te = np.asarray(split['train']), np.asarray(split['test'])
        assert len(set(te)) == len(te) and not np.intersect1d(tr, te).size
        regime = split['regime']
        if nested:
            axes = ['host'] if regime == 'host_cluster' else ['host', 'phage']
        else:
            axes = [] if regime == 'shared_entities_allowed' else (['host', 'phage'] if regime == 'both_clusters_excluded' else [regime.split('_')[0]])
        for axis in axes:
            assert set(rows.iloc[tr][axis+'_cluster']).isdisjoint(rows.iloc[te][axis+'_cluster'])
        for inner in split.get('inner', []):
            it, iv = inner['train'], inner['validation']
            assert set(it).issubset(tr) and set(iv).issubset(tr)
            assert set(it).isdisjoint(iv) and set(te).isdisjoint(it + iv)
            for axis in axes:
                assert set(rows.iloc[it][axis+'_cluster']).isdisjoint(rows.iloc[iv][axis+'_cluster'])
        part = prediction.loc[prediction.regime.eq(regime) & prediction.fold.eq(split['fold'])]
        if not nested:
            part = part.loc[part.seed.eq(split['seed'])]
        for method, scored in part.groupby('method'):
            scored = scored.sort_values('row_index')
            assert set(scored.row_index) == set(te)
            assert np.array_equal(rows.iloc[scored.row_index].label, scored.label)
            if nested:
                filename = f'{regime}_{split["fold"]}_{method}.joblib'
                cols = slice(None)
            elif method in ('phage_only', 'pair_only', 'full'):
                filename = f'{split["fold"]}_{split["seed"]}_{regime}_{method}.joblib'
                cols = {'phage_only': slice(0,212), 'pair_only':slice(424,448), 'full':slice(None)}[method]
            else:
                # Recalculate all training-only entity/prior controls independently.
                prior = rows.iloc[tr].label.mean()
                rates = {}
                for axis in ['host','phage']:
                    rate = rows.iloc[tr].groupby(axis).label.mean()
                    rates[axis] = rows.iloc[scored.row_index][axis].map(rate).fillna(prior).to_numpy()
                expected = {'training_prior': np.repeat(prior,len(scored)), 'host_positive_fraction':rates['host'],
                            'phage_positive_fraction':rates['phage'], 'additive_positive_fraction':(rates['host']+rates['phage'])/2}[method]
                np.testing.assert_allclose(expected, scored.score, atol=1e-12)
                continue
            model = joblib.load(directory / 'fitted_models' / filename)
            if nested:
                np.testing.assert_allclose(model['scale'].mean_, features[tr].mean(axis=0,dtype=np.float64), rtol=1e-10, atol=1e-10)
            assert_saved_predictions(scored.score,model.predict_proba(features[scored.row_index, cols])[:,1],model,
                                     portable_float32=portable_float32,strict_rtol=1e-7,strict_atol=1e-7)
            checked_reload += 1
        checked_splits += 1
    if not nested:
        # Every exclusion condition must use the same test rows and class counts.
        records = pd.DataFrame([{'fold':s['fold'],'seed':s['seed'],'test':tuple(s['test']),
                    'positive':int(rows.iloc[s['train']].label.sum()),'negative':int((1-rows.iloc[s['train']].label).sum())} for s in split_records])
        for _, part in records.groupby(['fold','seed']):
            assert part.test.nunique() == part.positive.nunique() == part.negative.nunique() == 1
    pooled = pd.read_csv(directory/'pooled_metrics.csv',float_precision='round_trip')
    for reported in pooled.itertuples():
        part = prediction.loc[prediction.regime.eq(reported.regime) & prediction.method.eq(reported.method)]
        if not nested: part = part.loc[part.seed.eq(reported.seed)]
        assert not part.row_index.duplicated().any()
        check_metric(part,reported)
    return {'complete':True,'checked_split_records':checked_splits,'reloaded_prediction_models':checked_reload,
            'independently_recalculated_pooled_metrics':len(pooled),'prediction_rows':len(prediction)}


def audit_doublets():
    directory = DATA / 'doublets'
    if not (directory/'summary.json').exists():return {'complete':False}
    summary = json.loads((directory/'summary.json').read_text())
    options = pd.read_csv(directory/'doublet_option_scores.csv')
    decisions = pd.read_csv(directory/'doublet_host_decisions.csv')
    # The exact recorded source alias predates this extension. No fuzzy matching.
    source = pd.read_csv(ROOT/'external/gaborieau_coli_phage_interactions_2023/dev/cocktails/data/cml_cocktails.csv',sep=';').replace({'bacteria':{'BCH857':'BCH857-72'}}).set_index('bacteria')
    compositions = pd.read_csv(ROOT/'external/gaborieau_coli_phage_interactions_2023/dev/cocktails/data/cocktails_composition.csv',sep=';').replace({'bacteria':{'BCH857':'BCH857-72'}}).set_index('bacteria')
    for (host, method), part in options.groupby(['bacteria','method']):
        assert set(part.option) == {'p1+p2','p1+p3','p2+p3'}
        expected = np.isclose(part.score,part.score.max(),atol=1e-12,rtol=0).astype(float)
        expected /= expected.sum()
        np.testing.assert_allclose(part.selection_probability,expected,atol=1e-12)
        for row in part.itertuples():
            assert row.assay_value == source.loc[host,row.option]
            assert row.active == int(row.assay_value > 0)
            p1,p2 = row.option.split('+')
            assert row.phage1 == compositions.loc[host,p1] and row.phage2 == compositions.loc[host,p2]
        outcome = np.dot(part.selection_probability,part.active)
        np.testing.assert_allclose(decisions.set_index('bacteria').loc[host,method],outcome)
    np.testing.assert_allclose(decisions.repaired_product.mean(),.78)
    np.testing.assert_allclose(decisions.random_choice_activity.mean(),224/300)
    cond = pd.read_csv(directory/'triplet_composition_metrics.csv')
    info = cond.loc[cond.method.eq('repaired_product') & cond.positive.gt(0) & cond.negative.gt(0)]
    assert len(info)==3 and info.n.sum()==72 and info.positive_negative_pairs.sum()==301
    np.testing.assert_allclose(np.average(info.conditional_auroc,weights=info.positive_negative_pairs),217/301)
    assert summary['informative_triplet_compositions']==3
    return {'complete':True,'options_verified':len(options),'unique_measured_doublets':300,
            'informative_triplet_compositions':3,'conditional_host_pairs':301,'repaired_active_choices':78,
            'uniform_activity':224/300}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    report = {'doublets':audit_doublets(), 'second_source':audit_models(DATA/'models/secondsource'),
              'nested_model_comparison':audit_models(DATA/'models/nested_comparison',nested=True)}
    report['all_three_complete'] = all(value['complete'] for value in report.values() if isinstance(value,dict))
    (OUT/'INDEPENDENT_ANALYSIS_CHECKS.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
