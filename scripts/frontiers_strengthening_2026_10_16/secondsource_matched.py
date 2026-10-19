"""Matched training-size/class-count hierarchy on already accessed StaphStudy."""
import os
for _v in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ[_v] = '2'
from datetime import datetime, timezone
import time
import json
import numpy as np
import pandas as pd
import joblib
from xgboost import XGBClassifier
from model_common import ROOT, OUT, sha, load_staph, metric, save_json, assert_split
from precisionphage.splits.cluster import combined_unseen_folds

DEST = OUT / 'secondsource'
SEEDS = [20261016,20261017,20261018]
REGIMES = [('shared_entities_allowed', []), ('host_cluster_excluded',['host']),
           ('phage_cluster_excluded',['phage']), ('both_clusters_excluded',['phage','host'])]

def fractions(train, test):
    prior = float(train.label.mean())
    values = {axis: test[axis].map(train.groupby(axis).label.mean()).fillna(prior).to_numpy() for axis in ['phage','host']}
    return {'phage_positive_fraction': values['phage'], 'host_positive_fraction': values['host'],
            'additive_positive_fraction': (values['phage']+values['host'])/2,
            'training_prior': np.repeat(prior,len(test))}

def main():
    started = time.monotonic()
    DEST.mkdir(parents=True, exist_ok=True)
    plan = {'created_utc': datetime.now(timezone.utc).isoformat(),
        'status':'Retrospective new analysis; StaphStudy outcomes already accessed in frozen transfer evaluation.',
        'question':'Does within-source dependence on shared entities recur in the dense StaphStudy matrix?',
        'clustering':'sourmash canonical 21-mer, 2048 bottom hashes, seed42, distance<=0.05 single-linkage; no label-guided cutoff selection',
        'test_design':'five diagonal dual-axis blocks seed42, eligibility at least3 positives3 negatives; same eligible test rows in every regime',
        'training_design':'exclude exact test rows; then indicated axes; each fold uses minimum class counts across regimes and three prespecified sampling seeds',
        'features':{'phage_only':212,'pair_only':24,'full':448},
        'model':{'n_estimators':400,'max_depth':5,'learning_rate':.05,'subsample':.8,'colsample_bytree':.8,'tree_method':'hist','n_jobs':2},
        'controls':'training-only positive fraction per exact phage/host, their mean, and global training prior; unknown entity uses training prior',
        'reporting':'all per-fold and pooled AUROC/AP, class prevalence, exposure; no row-wise p-values, seeds not independent replication',
        'limits':'This is within-source retraining in StaphStudy, not renewed external frozen validation; canonical groups differ from historical forward-strand clusters.',
        'seeds':SEEDS,
        'inputs':{p:sha(ROOT/p) for p in ['data/interim_v2/interactions.csv','data/interim_v2/external_staph_seq_pair_features.csv','data/interim_v2/external_staph_seq_pair_features.meta.json']}}
    plan_path = DEST / 'retrospective_protocol.json'
    if not plan_path.exists(): save_json(plan_path,plan)
    else: assert json.loads(plan_path.read_text())['inputs'] == plan['inputs']
    df,X,hashes,groups = load_staph()
    df.to_csv(DEST/'input_rows.csv',index=False)
    np.savez_compressed(DEST/'input_features.npz',X=X)
    pd.DataFrame(groups).to_csv(DEST/'canonical_clusters.csv',index=False)
    save_json(DEST/'genome_checksums.json',hashes)
    y=df.label.to_numpy()
    folds=list(combined_unseen_folds(df,'phage_cluster','host_cluster',5,42,3,3))
    save_json(DEST/'input_audit.json',{'rows':len(df),'positives':int(y.sum()),'phages':df.phage.nunique(),'hosts':df.host.nunique(),
        'phage_clusters':df.phage_cluster.nunique(),'host_clusters':df.host_cluster.nunique(),'eligible_blocks':len(folds),
        'test_rows':sum(len(f.test_idx) for f in folds),'complete_matrix':len(df)==df.phage.nunique()*df.host.nunique()})
    preds, reports, splits=[] ,[],[]
    for fi,fold in enumerate(folds):
        te=fold.test_idx
        nt=~np.isin(np.arange(len(df)),te)
        masks=[]
        for regime,axes in REGIMES:
            mask=nt.copy()
            for axis in axes: mask &= ~df[axis+'_cluster'].isin(df.iloc[te][axis+'_cluster']).to_numpy()
            masks.append(mask)
        counts={c:min(int((mask & (y==c)).sum()) for mask in masks) for c in [0,1]}
        assert min(counts.values())>0
        for seed in SEEDS:
            for ri,((regime,axes),mask) in enumerate(zip(REGIMES,masks)):
                rng=np.random.default_rng(seed+fi*100+ri)
                tr=np.concatenate([rng.choice(np.flatnonzero(mask & (y==c)),counts[c],replace=False) for c in [0,1]])
                assert_split(df,tr,te,axes)
                splits.append({'fold':fi,'seed':seed,'regime':regime,'train':tr.tolist(),'test':te.tolist()})
                scores=fractions(df.iloc[tr],df.iloc[te])
                for name,cols in [('phage_only',slice(0,212)),('pair_only',slice(424,448)),('full',slice(0,448))]:
                    model=XGBClassifier(**plan['model'],random_state=seed,eval_metric='logloss',verbosity=0)
                    model.fit(X[tr,cols],y[tr])
                    scores[name]=model.predict_proba(X[te,cols])[:,1]
                    modeldir=DEST/'fitted_models';modeldir.mkdir(exist_ok=True)
                    joblib.dump(model,modeldir/f'{fi}_{seed}_{regime}_{name}.joblib',compress=3)
                for name,p in scores.items():
                    part=df.iloc[te][['phage','host','label','phage_cluster','host_cluster']].copy()
                    part['row_index'],part['fold'],part['seed'],part['regime'],part['method'],part['score']=te,fi,seed,regime,name,p
                    preds.append(part)
                    reports.append({'fold':fi,'seed':seed,'regime':regime,'method':name,'train_n':len(tr),'train_positive':int(y[tr].sum()),
                        'test_host_seen_fraction':float(df.iloc[te].host_cluster.isin(df.iloc[tr].host_cluster).mean()),
                        'test_phage_seen_fraction':float(df.iloc[te].phage_cluster.isin(df.iloc[tr].phage_cluster).mean()),**metric(y[te],p)})
            print(f'Staph block {fi}, seed {seed} complete; elapsed {time.monotonic()-started:.1f}s',flush=True)
    pred=pd.concat(preds,ignore_index=True); report=pd.DataFrame(reports)
    pred.to_csv(DEST/'predictions.csv',index=False);report.to_csv(DEST/'fold_metrics.csv',index=False)
    save_json(DEST/'splits.json',splits)
    pooled=[]
    for (regime,method,seed),p in pred.groupby(['regime','method','seed']):
        assert not p.row_index.duplicated().any()
        f=report.query('regime==@regime and method==@method and seed==@seed')
        pooled.append({'regime':regime,'method':method,'seed':int(seed),**metric(p.label,p.score),
            'macro_fold_auroc':float(f.auroc.mean()),'within_fold_pair_weighted_auroc':float(np.average(f.auroc,weights=f.positive*(f.n-f.positive)))})
    pd.DataFrame(pooled).to_csv(DEST/'pooled_metrics.csv',index=False)
    save_json(DEST/'completion.json',{'complete':True,'elapsed_seconds':time.monotonic()-started,'protocol_sha256':sha(plan_path),'fits':len(splits)*3})
    print(pd.DataFrame(pooled).groupby(['method','regime']).auroc.mean().to_string(),flush=True)

if __name__=='__main__':main()
