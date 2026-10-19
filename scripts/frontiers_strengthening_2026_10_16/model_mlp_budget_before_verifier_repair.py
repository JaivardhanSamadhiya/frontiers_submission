"""Predeclared retrospective optimization-budget sensitivity, no new data."""
import os
for _key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_key]='2'
from pathlib import Path
import hashlib
import json
import time
import warnings
from datetime import datetime,timezone
import joblib
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPClassifier
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
from model_common import ROOT,OUT,sha,metric,save_json,assert_split

BASE=OUT/'nested_comparison'
DEST=OUT/'mlp_budget_sensitivity'
PLAN=Path(__file__).with_name('MLP_BUDGET_AMENDMENT.md')


def estimator(params,seed):
    params=dict(params);params['hidden_layer_sizes']=tuple(params['hidden_layer_sizes'])
    return Pipeline([('scale',StandardScaler()),('model',MLPClassifier(**params,random_state=seed,
        max_iter=1000,tol=.0001,n_iter_no_change=20,early_stopping=False,batch_size=64,
        learning_rate_init=.001,solver='adam'))])


def fit(model,X,y):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always',ConvergenceWarning)
        with threadpool_limits(limits=2):model.fit(X,y)
    return [str(w.message) for w in caught]


def main():
    started=time.monotonic();DEST.mkdir(parents=True,exist_ok=True)
    assert json.loads((BASE/'completion.json').read_text())['complete']
    protocol=json.loads((BASE/'retrospective_protocol.json').read_text())
    candidates=protocol['candidates']['MLP']
    inputs={name:sha(BASE/name) for name in ['input_rows.csv','input_features.npz','splits.json','pooled_metrics.csv','retrospective_protocol.json']}
    plan={'amendment_sha256':sha(PLAN),'script_sha256':sha(Path(__file__)),
          'inputs_sha256':inputs,'created_utc':datetime.now(timezone.utc).isoformat(),
          'max_iter':1000,'training_loss_tol':.0001,'training_loss_patience':20,
          'early_stopping':False,'candidates':candidates,'status':'retrospective outcome-exposed sensitivity'}
    lock=DEST/'protocol_lock.json'
    if lock.exists():
        previous=json.loads(lock.read_text())
        for key in ['amendment_sha256','script_sha256','inputs_sha256','candidates']:
            assert previous[key]==plan[key],key
    else:save_json(lock,plan)
    df=pd.read_csv(BASE/'input_rows.csv');X=np.load(BASE/'input_features.npz')['X'];y=df.label.to_numpy()
    splits=json.loads((BASE/'splits.json').read_text());allpred=[];allinner=[];allouter=[];selections=[]
    for split in splits:
        regime,fold=split['regime'],split['fold'];tr,te=np.array(split['train']),np.array(split['test'])
        axes=['host'] if regime=='host_cluster' else ['phage','host']
        assert_split(df,tr,te,axes)
        work=DEST/f'{regime}_{fold}';work.mkdir(exist_ok=True)
        done=work/'complete.json'
        if done.exists():
            saved=json.loads(done.read_text())
            for name,digest in saved['sha256'].items():assert sha(work/name)==digest
        else:
            innerrows=[];means=[]
            for ci,params in enumerate(candidates):
                values=[]
                for ii,inner in enumerate(split['inner']):
                    it,iv=np.array(inner['train']),np.array(inner['validation'])
                    assert_split(df,it,iv,axes)
                    assert set(it).issubset(tr) and set(iv).issubset(tr)
                    assert not np.intersect1d(te,np.r_[it,iv]).size
                    model=estimator(params,20261016+fold);warn=fit(model,X[it],y[it])
                    np.testing.assert_allclose(model['scale'].mean_,X[it].mean(0,dtype=np.float64),rtol=1e-10,atol=1e-10)
                    met=metric(y[iv],model.predict_proba(X[iv])[:,1]);values.append(met['auroc'])
                    innerrows.append({'regime':regime,'outer_fold':fold,'candidate':ci,'inner_fold':ii,
                        'epochs':model['model'].n_iter_,'warnings':' | '.join(warn),**met})
                    save_json(work/f'training_trace_candidate{ci}_inner{ii}.json',model['model'].loss_curve_)
                means.append(float(np.mean(values)))
            chosen=int(np.argmax(means));model=estimator(candidates[chosen],20261016+fold)
            warn=fit(model,X[tr],y[tr]);p=model.predict_proba(X[te])[:,1]
            np.testing.assert_allclose(model['scale'].mean_,X[tr].mean(0,dtype=np.float64),rtol=1e-10,atol=1e-10)
            joblib.dump(model,work/'model.joblib',compress=3)
            np.testing.assert_allclose(joblib.load(work/'model.joblib').predict_proba(X[te])[:,1],p,rtol=0,atol=0)
            part=df.iloc[te][['phage','host','label','phage_cluster','host_cluster']].copy()
            part['row_index'],part['regime'],part['fold'],part['method'],part['score']=te,regime,fold,'MLP_extended',p
            part.to_csv(work/'predictions.csv',index=False)
            pd.DataFrame(innerrows).to_csv(work/'inner_metrics.csv',index=False)
            save_json(work/'outer_metrics.json',{'regime':regime,'fold':fold,'epochs':model['model'].n_iter_,
                'warnings':' | '.join(warn),**metric(y[te],p)})
            save_json(work/'selection.json',{'chosen':chosen,'parameters':candidates[chosen],'mean_inner_aurocs':means})
            save_json(work/'training_trace_outer.json',model['model'].loss_curve_)
            save_json(done,{'sha256':{p.name:sha(p) for p in work.iterdir() if p.is_file() and p.name!='complete.json'}})
        pred=pd.read_csv(work/'predictions.csv',float_precision='round_trip')
        model=joblib.load(work/'model.joblib')
        np.testing.assert_array_equal(pred.row_index,te)
        np.testing.assert_allclose(model.predict_proba(X[te])[:,1],pred.score,rtol=1e-12,atol=1e-12)
        selected=json.loads((work/'selection.json').read_text());assert selected['chosen']==int(np.argmax(selected['mean_inner_aurocs']))
        inner=pd.read_csv(work/'inner_metrics.csv',float_precision='round_trip')
        np.testing.assert_allclose(inner.groupby('candidate').auroc.mean(),selected['mean_inner_aurocs'])
        allpred.append(pred);allinner.append(inner);allouter.append(json.loads((work/'outer_metrics.json').read_text()))
        selections.append({'regime':regime,'fold':fold,**selected})
        save_json(DEST/'status.json',{'fold':fold,'regime':regime,'completed_outer_folds':len(allouter),
            'updated_utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid()})
        print(f'MLP sensitivity {regime}/{fold} complete; elapsed {time.monotonic()-started:.1f}s',flush=True)
    preds=pd.concat(allpred,ignore_index=True);inner=pd.concat(allinner,ignore_index=True)
    preds.to_csv(DEST/'predictions.csv',index=False);inner.to_csv(DEST/'inner_metrics.csv',index=False)
    pd.DataFrame(allouter).to_csv(DEST/'outer_metrics.csv',index=False);save_json(DEST/'selections.json',selections)
    summary=[]
    old=pd.read_csv(BASE/'pooled_metrics.csv',float_precision='round_trip')
    for regime,part in preds.groupby('regime'):
        assert not part.row_index.duplicated().any()
        folds=[v for v in allouter if v['regime']==regime]
        summary.append({'regime':regime,'method':'MLP_extended',**metric(part.label,part.score),
            'macro_fold_auroc':float(np.mean([f['auroc'] for f in folds])),
            'within_fold_pair_weighted_auroc':float(np.average([f['auroc'] for f in folds],weights=[f['positive']*(f['n']-f['positive']) for f in folds]))})
    pd.concat([old,pd.DataFrame(summary)],ignore_index=True).to_csv(DEST/'pooled_metrics.csv',index=False)
    save_json(DEST/'completion.json',{'complete':True,'outer_folds':10,'outer_models_reloaded':10,
        'inner_fits':len(inner),'elapsed_seconds_this_resume':time.monotonic()-started,
        'amendment_sha256':sha(PLAN),'protocol_lock_sha256':sha(lock),
        'max_epoch_warnings':int(inner.warnings.fillna('').str.contains('Maximum iterations').sum()+sum('Maximum iterations' in f['warnings'] for f in allouter)),
        'scope':'Retrospective training-budget sensitivity; all prior 150-epoch results retained',
        'outputs_sha256':{p.name:sha(p) for p in DEST.iterdir() if p.is_file() and p.name not in ['completion.json','status.json']}})
    print(pd.DataFrame(summary).to_string(index=False),flush=True)


if __name__=='__main__':main()
