"""Fresh, fully logged grouped nested GBM/MLP comparison; retrospective."""
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[_v]='2'
import json
import time
import warnings
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import joblib
import sklearn
import xgboost
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold
from sklearn.neural_network import MLPClassifier
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits
from xgboost import XGBClassifier
from model_common import ROOT,OUT,sha,load_ncbi,metric,save_json,assert_split
from precisionphage.splits.cluster import combined_unseen_folds

DEST=OUT/'nested_comparison'
CANDIDATES={
 'GBM':[{'n_estimators':200,'max_depth':3}, {'n_estimators':400,'max_depth':3}, {'n_estimators':400,'max_depth':5}],
 'MLP':[{'hidden_layer_sizes':[64],'alpha':.001}, {'hidden_layer_sizes':[128,64],'alpha':.001}, {'hidden_layer_sizes':[128,64],'alpha':.01}]
}

def folds_for(df,regime,n,seed):
    if regime=='host_cluster':
        folds=list(GroupKFold(n_splits=n).split(np.zeros(len(df)),groups=df.host_cluster))
    else:
        folds=[(f.train_idx,f.test_idx) for f in combined_unseen_folds(df,'phage_cluster','host_cluster',n,seed,3,3)]
    return [(np.asarray(tr),np.asarray(te)) for tr,te in folds if df.iloc[tr].label.nunique()==2 and df.iloc[te].label.nunique()==2]

def estimator(family,params,seed):
    if family=='GBM':
        model=XGBClassifier(**params,learning_rate=.05,subsample=.8,colsample_bytree=.8,
            random_state=seed,tree_method='hist',n_jobs=2,eval_metric='logloss',verbosity=0)
    else:
        settings=dict(params);settings['hidden_layer_sizes']=tuple(settings['hidden_layer_sizes'])
        model=MLPClassifier(**settings,random_state=seed,max_iter=150,early_stopping=False,
            n_iter_no_change=151,tol=0,batch_size=64,learning_rate_init=.001,solver='adam')
    return Pipeline([('scale',StandardScaler()),('model',model)])

def fit(model,X,y):
    with warnings.catch_warnings(record=True) as records:
        warnings.simplefilter('always',ConvergenceWarning)
        with threadpool_limits(limits=2):model.fit(X,y)
    return [str(w.message) for w in records]

def main():
    started=time.monotonic();DEST.mkdir(parents=True,exist_ok=True)
    plan={'created_utc':datetime.now(timezone.utc).isoformat(),
      'status':'Retrospective fresh model-family comparison; original held-out outcomes already known. Does not reconstruct lost historical tuning or justify a historical model selection.',
      'input':'NCBI_HR original 1947 rows and compatible full448 columns (424node+24pair), preserved historical forward-strand .05 groups',
      'outer':'Five group-size-balanced host-cluster folds; five original diagonal dual-exclusion blocks seed42. Host and dual regimes have different held-out cohorts; compare families only within a regime.',
      'inner':'Three grouped folds of outer training only, matching host-only or dual exclusion; no row-stratified fallback; ineligible one-class folds reported and omitted; require>=2 inner folds.',
      'preprocessing':'StandardScaler separately fitted inside each inner train and outer refit; no outcome-based feature selection, calibration, or random row early stopping.',
      'selection':'Mean inner-fold AUROC; exact tie resolved by candidate list order; outer outcomes never select a model.',
      'budgets':'Three predeclared candidates per family evaluated on identical inner folds. Equal candidate/split budget, not equal wall-clock or parameter count. MLP150 full epochs no early stop; GBM200/400 trees fixed.',
      'candidates':CANDIDATES,'seed':20261016,
      'reporting':'All predictions, train/test indices, inner metrics, chosen candidates, scalers/models, warnings and timing retained. Pooled and macro-fold AUROC/AP descriptive; no independent-row significance tests.',
      'inputs':{p:sha(ROOT/p) for p in ['data/interim_v2/interactions_modelable.csv','data/cache_v2/node_features.npz','data/cache_v2/node_features.meta.json','data/interim_v2/seq_pair_features.csv','data/cache_v2/clusters_k21_d0.05.json']},
      'versions':{'sklearn':sklearn.__version__,'xgboost':xgboost.__version__}}
    pp=DEST/'retrospective_protocol.json'
    if pp.exists():assert json.loads(pp.read_text())['inputs']==plan['inputs']
    else:save_json(pp,plan)
    df,X=load_ncbi();y=df.label.to_numpy();assert X.shape==(1947,448)
    df.to_csv(DEST/'input_rows.csv',index=False);np.savez_compressed(DEST/'input_features.npz',X=X)
    allsplits=[];inner_metrics=[];outer_metrics=[];predictions=[];selection=[]
    for regime in ['host_cluster','dual_cluster']:
      axes=['host'] if regime=='host_cluster' else ['host','phage']
      outer=folds_for(df,regime,5,42)
      for fi,(tr,te) in enumerate(outer):
        assert_split(df,tr,te,axes)
        inner=folds_for(df.iloc[tr].reset_index(drop=True),regime,3,20261016+fi)
        if len(inner)<2:raise RuntimeError(f'{regime} outer{fi} has only {len(inner)} eligible inner folds; no fallback permitted')
        inner_global=[(tr[it],tr[iv]) for it,iv in inner]
        for it,iv in inner_global:
            assert_split(df,it,iv,axes)
            assert set(it).issubset(tr) and set(iv).issubset(tr)
            assert not np.intersect1d(te,np.r_[it,iv]).size
        allsplits.append({'regime':regime,'fold':fi,'train':tr.tolist(),'test':te.tolist(),
                          'inner':[{'train':it.tolist(),'validation':iv.tolist()} for it,iv in inner_global]})
        for family,candidates in CANDIDATES.items():
          candidate_aucs=[]
          for ci,params in enumerate(candidates):
            vals=[]
            for ii,(it,iv) in enumerate(inner_global):
              t=time.monotonic();model=estimator(family,params,20261016+fi)
              warn=fit(model,X[it],y[it]);p=model.predict_proba(X[iv])[:,1]
              np.testing.assert_allclose(model['scale'].mean_,X[it].mean(axis=0,dtype=np.float64),rtol=1e-10,atol=1e-10)
              met=metric(y[iv],p);vals.append(met['auroc'])
              inner_metrics.append({'regime':regime,'outer_fold':fi,'family':family,'candidate':ci,'inner_fold':ii,
                  'train_n':len(it),'elapsed_seconds':time.monotonic()-t,'warnings':' | '.join(warn),**met})
            candidate_aucs.append(float(np.mean(vals)))
          chosen=int(np.argmax(candidate_aucs));t=time.monotonic()
          model=estimator(family,candidates[chosen],20261016+fi)
          warn=fit(model,X[tr],y[tr]);p=model.predict_proba(X[te])[:,1]
          np.testing.assert_allclose(model['scale'].mean_,X[tr].mean(axis=0,dtype=np.float64),rtol=1e-10,atol=1e-10)
          modeldir=DEST/'fitted_models';modeldir.mkdir(exist_ok=True)
          joblib.dump(model,modeldir/f'{regime}_{fi}_{family}.joblib',compress=3)
          part=df.iloc[te][['phage','host','label','phage_cluster','host_cluster']].copy()
          part['row_index'],part['regime'],part['fold'],part['method'],part['score']=te,regime,fi,family,p
          predictions.append(part)
          outer_metrics.append({'regime':regime,'fold':fi,'family':family,'train_n':len(tr),'chosen_candidate':chosen,
              'refit_seconds':time.monotonic()-t,'warnings':' | '.join(warn),**metric(y[te],p)})
          selection.append({'regime':regime,'fold':fi,'family':family,'chosen':chosen,'parameters':candidates[chosen],
                            'mean_inner_aurocs':candidate_aucs,'inner_folds':len(inner)})
          print(f'{regime} outer{fi} {family}: chosen{chosen} AUC{outer_metrics[-1]["auroc"]:.4f}; elapsed {time.monotonic()-started:.1f}s',flush=True)
        # Durable progress after each complete paired outer fold.
        save_json(DEST/'splits.json',allsplits);save_json(DEST/'selections.json',selection)
        pd.DataFrame(inner_metrics).to_csv(DEST/'inner_metrics.csv',index=False)
        pd.DataFrame(outer_metrics).to_csv(DEST/'outer_metrics.csv',index=False)
        pd.concat(predictions,ignore_index=True).to_csv(DEST/'predictions.csv',index=False)
    pred=pd.concat(predictions,ignore_index=True);foldmet=pd.DataFrame(outer_metrics);rows=[]
    for (regime,family),part in pred.groupby(['regime','method']):
        assert not part.row_index.duplicated().any()
        fm=foldmet.query('regime==@regime and family==@family')
        rows.append({'regime':regime,'method':family,**metric(part.label,part.score),'macro_fold_auroc':float(fm.auroc.mean()),
                    'within_fold_pair_weighted_auroc':float(np.average(fm.auroc,weights=fm.positive*(fm.n-fm.positive)))})
    pd.DataFrame(rows).to_csv(DEST/'pooled_metrics.csv',index=False)
    save_json(DEST/'completion.json',{'complete':True,'elapsed_seconds':time.monotonic()-started,'protocol_sha256':sha(pp),
        'inner_fits':len(inner_metrics),'outer_refits':len(outer_metrics),'convergence_warning_count':sum(bool(x['warnings']) for x in inner_metrics+outer_metrics)})
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':main()
