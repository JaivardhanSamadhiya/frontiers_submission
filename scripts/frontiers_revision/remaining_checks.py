"""Retrospective reviewer checks using existing fixed data; CPU only.

No historical result is overwritten. No external outcome selects parameters.
Run one of: cutoffs, ablations, recalibration.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.isotonic import IsotonicRegression
from xgboost import XGBClassifier

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from precisionphage.splits import leave_one_group_out, combined_unseen_folds
OUT=ROOT/'data/results_v2/frontiers_revision_2026_10_14/remaining_checks'
SEEDS=[20261005,20261006,20261007]


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def design():
    spec=importlib.util.spec_from_file_location('matched',ROOT/'experiments/35_matched_cohort_leakage_audit.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    df,x=module.load_design()
    assert x.shape==(1947,448)
    return df,x


def score_model(x,y,tr,te,seed):
    assert len(np.intersect1d(tr,te))==0
    scale=StandardScaler().fit(x[tr])
    model=XGBClassifier(n_estimators=400,max_depth=5,learning_rate=.05,
        subsample=.8,colsample_bytree=.8,random_state=seed,tree_method='hist',
        n_jobs=2,eval_metric='logloss',verbosity=0)
    model.fit(scale.transform(x[tr]),y[tr])
    return model.predict_proba(scale.transform(x[te]))[:,1]


def metrics(y,s):
    bins=np.minimum((np.clip(s,0,1)*10).astype(int),9)
    ece=sum(np.mean(bins==b)*abs(np.mean(y[bins==b])-np.mean(s[bins==b]))
            for b in range(10) if (bins==b).any())
    return {'n':len(y),'positive':int(np.sum(y)),'prevalence':float(np.mean(y)),
            'auroc':float(roc_auc_score(y,s)) if len(np.unique(y))==2 else None,
            'ap':float(average_precision_score(y,s)),
            'brier':float(brier_score_loss(y,s)),'ece':float(ece)}


def finish(name,inputs,details):
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),
        'status':'retrospective reviewer sensitivity; not independently frozen validation',
        'mode':name,'script_sha256':sha(__file__),
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},**details}
    (OUT/(name+'_manifest.json')).write_text(json.dumps(manifest,indent=2))


def cutoffs():
    df,x=design();y=df.label.to_numpy();old=df.copy()
    folds=list(combined_unseen_folds(old,'phage_cluster','host_cluster',5,42))
    assert len(folds)==5
    cache=ROOT/'data/results_v2/canonical_sensitivity'
    distances={axis:np.load(cache/(axis+'_distances_2048.npz'),allow_pickle=False) for axis in ['host','phage']}
    regimes=['shared_entities_allowed','host_cluster_excluded','phage_cluster_excluded','both_clusters_excluded']
    predictions=[];counts=[];train_counts=[]
    # 0.05 is reproduced from the preserved identical-design output, not refitted.
    existing=pd.read_csv(cache/'predictions.csv');existing['cutoff']=.05
    predictions.append(existing)
    for threshold in [.025,.05,.1]:
        for axis,z in distances.items():
            n,labels=connected_components(csr_matrix(z['distance']<=threshold),directed=False)
            df[axis+'_cluster']=df[axis].map(dict(zip(z['names'],labels)))
            assert df[axis+'_cluster'].notna().all()
            counts.append({'cutoff':threshold,'axis':axis,'entities':len(labels),'clusters':n})
        if threshold==.05:continue
        for fi,fold in enumerate(folds):
            te=fold.test_idx
            h=df.host_cluster.isin(df.iloc[te].host_cluster).to_numpy()
            p=df.phage_cluster.isin(df.iloc[te].phage_cluster).to_numpy()
            nt=~np.isin(np.arange(len(df)),te)
            masks=[nt,nt & ~h,nt & ~p,nt & ~h & ~p]
            classes={c:min(int(np.sum(mask & (y==c))) for mask in masks) for c in [0,1]}
            if not all(classes.values()):
                train_counts.append({'cutoff':threshold,'fold':fi,'status':'not estimable: a matched training class is empty',**{f'class_{c}':v for c,v in classes.items()}})
                continue
            for seed in SEEDS:
                for ri,(regime,mask) in enumerate(zip(regimes,masks)):
                    rng=np.random.default_rng(seed+fi*100+ri)
                    tr=np.concatenate([rng.choice(np.flatnonzero(mask & (y==c)),classes[c],replace=False) for c in [0,1]])
                    if ri in [1,3]:assert set(df.iloc[tr].host_cluster).isdisjoint(df.iloc[te].host_cluster)
                    if ri in [2,3]:assert set(df.iloc[tr].phage_cluster).isdisjoint(df.iloc[te].phage_cluster)
                    s=score_model(x,y,tr,te,seed)
                    q=df.iloc[te][['phage','host','label','phage_cluster','host_cluster']].copy()
                    q['row_index']=te;q['fold']=fi;q['seed']=seed;q['regime']=regime;q['score']=s;q['cutoff']=threshold
                    predictions.append(q)
                    train_counts.append({'cutoff':threshold,'fold':fi,'seed':seed,'regime':regime,
                        'status':'completed','train_n':len(tr),'train_positive':classes[1],'test_n':len(te)})
                print(f'cutoff {threshold} block {fi} seed {seed}',flush=True)
    pred=pd.concat(predictions,ignore_index=True)
    results=[]
    for (threshold,regime,seed),q in pred.groupby(['cutoff','regime','seed']):
        assert not q.row_index.duplicated().any()
        results.append({'cutoff':threshold,'regime':regime,'seed':int(seed),**metrics(q.label.to_numpy(),q.score.to_numpy())})
    pred.to_csv(OUT/'cutoff_predictions.csv',index=False)
    pd.DataFrame(results).to_csv(OUT/'cutoff_summary.csv',index=False)
    pd.DataFrame(counts).to_csv(OUT/'cutoff_cluster_counts.csv',index=False)
    pd.DataFrame(train_counts).to_csv(OUT/'cutoff_training_counts.csv',index=False)
    finish('cutoffs',[cache/'host_distances_2048.npz',cache/'phage_distances_2048.npz',cache/'predictions.csv'],
        {'cutoffs':[.025,.05,.1],'seeds':SEEDS,'sketch':'canonical sourmash k21, 2048 hashes',
         'cohort':'same original 398 dual-block test rows; training class counts matched within each cutoff, not across cutoffs',
         'limitation':'threshold changes grouping and eligible training counts; descriptive, not ANI or a deployment-threshold selection'})
    print(pd.DataFrame(results).groupby(['cutoff','regime']).auroc.agg(['min','max']).to_string(),flush=True)


def ablations():
    df,x=design();y=df.label.to_numpy()
    folds={
        'loso_species':list(leave_one_group_out(df,'host_species','loso',3,3)),
        'host_cluster':list(leave_one_group_out(df,'host_cluster','host_cluster',3,3)),
        'phage_cluster':list(leave_one_group_out(df,'phage_cluster','phage_cluster',3,3)),
        'combined_unseen':list(combined_unseen_folds(df,'phage_cluster','host_cluster',5,42,3,3))}
    blocks={'full':slice(0,448),'phage_only':slice(0,212),'host_only':slice(212,424),'pair_only':slice(424,448)}
    predictions=[];fold_metrics=[]
    for regime,items in folds.items():
        for i,fold in enumerate(items):
            for method,cols in blocks.items():
                s=score_model(x[:,cols],y,fold.train_idx,fold.test_idx,42)
                q=df.iloc[fold.test_idx][['phage','host','label']].copy()
                q['row_index']=fold.test_idx;q['fold']=fold.name;q['regime']=regime;q['method']=method;q['score']=s
                predictions.append(q)
                fold_metrics.append({'regime':regime,'fold':fold.name,'method':method,**metrics(q.label.to_numpy(),s)})
            print(f'ablation {regime} {i+1}/{len(items)}',flush=True)
    pred=pd.concat(predictions,ignore_index=True)
    results=[]
    for (regime,method),q in pred.groupby(['regime','method']):
        assert not q.row_index.duplicated().any()
        results.append({'regime':regime,'method':method,**metrics(q.label.to_numpy(),q.score.to_numpy())})
    pred.to_csv(OUT/'ablation_predictions.csv',index=False)
    pd.DataFrame(fold_metrics).to_csv(OUT/'ablation_folds.csv',index=False)
    pd.DataFrame(results).to_csv(OUT/'ablation_summary.csv',index=False)
    finish('ablations',[ROOT/'data/cache_v2/node_features.npz',ROOT/'data/cache_v2/node_features.meta.json',
        ROOT/'data/interim_v2/interactions_modelable.csv',ROOT/'data/interim_v2/seq_pair_features.csv',
        ROOT/'data/cache_v2/clusters_k21_d0.05.json'],
        {'features':{'full':448,'phage_only':212,'host_only':212,'pair_only':24},
         'settings':'fixed XGBoost400,depth5,lr.05,subsample/colsample.8,seed42; same rows/splits/scaler fitting',
         'limitation':'new full control refitted alongside ablations; not a replacement for saved historical model, no causal attribution or independent-row inference'})
    print(pd.DataFrame(results).to_string(index=False),flush=True)


def recalibration():
    path=ROOT/'data/results_v2/external_staph_predictions.csv'
    df=pd.read_csv(path);y=df.label.to_numpy(int)
    assert len(df)==1053 and df.host.nunique()==27
    raw=df.prob_full.to_numpy(float)
    logit=np.log(np.clip(raw,1e-6,1-1e-6)/np.clip(1-raw,1e-6,1-1e-6))[:,None]
    splits=list(StratifiedGroupKFold(n_splits=5,shuffle=True,random_state=20261014).split(logit,y,df.host))
    records=[];count_rows=[]
    for seed in [20261014,20261015,20261016,20261017,20261018]:
        for fraction in [.1,.25,.5,.75,1.]:
            for fi,(development,test) in enumerate(splits):
                train_hosts=sorted(df.iloc[development].host.unique())
                rng=np.random.default_rng(seed+fi*100)
                rng.shuffle(train_hosts)
                chosen=set(train_hosts[:max(1,int(np.ceil(fraction*len(train_hosts))))])
                train=development[df.iloc[development].host.isin(chosen).to_numpy()]
                assert set(df.iloc[train].host).isdisjoint(df.iloc[test].host)
                prior=float(np.mean(y[train]))
                for method in ['uncalibrated','target_prior','platt','isotonic']:
                    if method=='uncalibrated':s=raw[test]
                    elif method=='target_prior' or len(np.unique(y[train]))<2:s=np.full(len(test),prior)
                    elif method=='platt':
                        fit=LogisticRegression(C=1.,solver='lbfgs',random_state=seed).fit(logit[train],y[train])
                        s=fit.predict_proba(logit[test])[:,1]
                    else:s=IsotonicRegression(out_of_bounds='clip').fit(raw[train],y[train]).predict(raw[test])
                    q=df.iloc[test][['phage','host','label']].copy()
                    q['row_index']=test;q['seed']=seed;q['fraction']=fraction;q['fold']=fi;q['method']=method;q['score']=s
                    records.append(q)
                count_rows.append({'seed':seed,'fraction':fraction,'fold':fi,'train_hosts':len(chosen),
                    'train_rows':len(train),'test_hosts':df.iloc[test].host.nunique(),'test_rows':len(test),
                    'train_positive':int(y[train].sum()),'single_class_fallback':len(np.unique(y[train]))<2})
    pred=pd.concat(records,ignore_index=True);summaries=[]
    for (seed,fraction,method),q in pred.groupby(['seed','fraction','method']):
        assert len(q)==1053 and not q.row_index.duplicated().any()
        summaries.append({'seed':int(seed),'fraction':fraction,'method':method,**metrics(q.label.to_numpy(),q.score.to_numpy())})
    pred.to_csv(OUT/'recalibration_predictions.csv',index=False)
    pd.DataFrame(summaries).to_csv(OUT/'recalibration_summary.csv',index=False)
    pd.DataFrame(count_rows).to_csv(OUT/'recalibration_label_budgets.csv',index=False)
    finish('recalibration',[path],{'seeds':[20261014,20261015,20261016,20261017,20261018],
        'fractions':[.1,.25,.5,.75,1.],'outer_splits':5,'grouping':'27 exact target host reference identifiers, not sequence clusters',
        'target':'label-assisted target-source adaptation on hosts disjoint from calibration; phages may recur',
        'single_class_rule':'constant observed calibration prior',
        'limitation':'target labels now train calibrators; not retraining-free transfer, clinical probability validation, or an independent frozen test; no setting chosen from these results'})
    print(pd.DataFrame(summaries).groupby(['fraction','method'])[['auroc','brier','ece']].mean().to_string(),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['cutoffs','ablations','recalibration'])
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    {'cutoffs':cutoffs,'ablations':ablations,'recalibration':recalibration}[args.mode]()
