"""Verify the completed higher-budget sensitivity without fitting a model."""
from pathlib import Path
import json
import joblib
import numpy as np
import pandas as pd
from model_common import ROOT,OUT,sha,metric,assert_split,save_json


def verify(write=True):
    base=OUT/'nested_comparison';dest=OUT/'mlp_budget_sensitivity'
    complete=json.loads((dest/'completion.json').read_text());assert complete['complete']
    for name,digest in complete['outputs_sha256'].items():assert sha(dest/name)==digest
    assert sha(Path(__file__).with_name('model_mlp_budget_sensitivity.py'))==complete['executed_script_sha256']
    lock=json.loads((dest/'protocol_lock.json').read_text())
    repair=json.loads((dest/'technical_verifier_repair.json').read_text())
    assert repair['previous_script_sha256']==lock['script_sha256']==sha(Path(__file__).with_name('model_mlp_budget_before_verifier_repair.py'))
    assert repair['revised_script_sha256']==complete['executed_script_sha256']
    assert sha(Path(__file__).with_name('MLP_BUDGET_AMENDMENT.md'))==lock['amendment_sha256']
    for name,digest in lock['inputs_sha256'].items():assert sha(base/name)==digest
    X=np.load(base/'input_features.npz')['X'];df=pd.read_csv(base/'input_rows.csv')
    splits=json.loads((base/'splits.json').read_text())
    scored=pd.read_csv(dest/'predictions.csv',float_precision='round_trip')
    inners=pd.read_csv(dest/'inner_metrics.csv',float_precision='round_trip')
    selections=json.loads((dest/'selections.json').read_text())
    epoch_counts=[]
    for split in splits:
        regime,fold=split['regime'],split['fold'];tr,te=np.array(split['train']),np.array(split['test'])
        axes=['host'] if regime=='host_cluster' else ['host','phage']
        assert_split(df,tr,te,axes)
        for inner in split['inner']:
            assert_split(df,inner['train'],inner['validation'],axes)
            assert set(inner['train']).issubset(tr) and set(inner['validation']).issubset(tr)
            assert set(te).isdisjoint(inner['train']+inner['validation'])
        work=dest/f'{regime}_{fold}';model=joblib.load(work/'model.joblib')
        assert model['model'].max_iter==1000 and model['model'].early_stopping is False
        assert model['model'].tol==.0001 and model['model'].n_iter_no_change==20
        np.testing.assert_allclose(model['scale'].mean_,X[tr].mean(0,dtype=np.float64),rtol=1e-10,atol=1e-10)
        pred=scored.loc[scored.regime.eq(regime)&scored.fold.eq(fold)]
        np.testing.assert_array_equal(pred.row_index,te)
        np.testing.assert_allclose(model.predict_proba(X[te])[:,1],pred.score,rtol=1e-6,atol=1e-7)
        selected=next(s for s in selections if s['regime']==regime and s['fold']==fold)
        inner=inners.loc[inners.regime.eq(regime)&inners.outer_fold.eq(fold)]
        means=inner.groupby('candidate').auroc.mean().to_numpy()
        np.testing.assert_allclose(means,selected['mean_inner_aurocs'])
        assert selected['chosen']==int(np.argmax(means))
        assert selected['parameters']==lock['candidates'][selected['chosen']]
        assert len(inner)==9 and inner.epochs.between(1,1000).all()
        for row in inner.itertuples():
            trace=json.loads((work/f'training_trace_candidate{row.candidate}_inner{row.inner_fold}.json').read_text())
            assert len(trace)==row.epochs and np.isfinite(trace).all()
        epoch_counts.extend(inner.epochs.tolist());epoch_counts.append(model['model'].n_iter_)
        checkpoint=json.loads((work/'complete.json').read_text())
        for name,digest in checkpoint['sha256'].items():assert sha(work/name)==digest
    pooled=pd.read_csv(dest/'pooled_metrics.csv',float_precision='round_trip')
    for regime,pred in scored.groupby('regime'):
        expected=pooled.loc[pooled.regime.eq(regime)&pooled.method.eq('MLP_extended')].iloc[0]
        assert not pred.row_index.duplicated().any()
        for name,value in metric(pred.label,pred.score).items():np.testing.assert_allclose(expected[name],value,rtol=1e-6,atol=1e-7)
    report={'passed':True,'outer_models_reloaded':10,'inner_score_records_checked':len(inners),
            'epochs_min':int(min(epoch_counts)),'epochs_max':int(max(epoch_counts)),
            'fits_reaching_1000':int(sum(e==1000 for e in epoch_counts)),
            'source_sha256':complete['executed_script_sha256'],'completion_sha256':sha(dest/'completion.json'),
            'scope':'Retrospective fixed-fold training-budget sensitivity, not exhaustive family optimization'}
    if write:save_json(dest/'verification.json',report)
    return report


if __name__=='__main__':print(json.dumps(verify(),indent=2))
