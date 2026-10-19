"""Independent saved-artifact checks for the fresh model analyses."""
import argparse
import json
import sys
import time
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import joblib
from model_common import ROOT,OUT,sha,assert_split,metric,save_json
from secondsource_matched import fractions
from model_nested_comparison import estimator,folds_for,CANDIDATES
from prediction_portability_policy import assert_saved_predictions

if not __debug__:
    raise RuntimeError('Run without Python optimization; verification assertions are mandatory')

class LeakageGuards(unittest.TestCase):
    def test_entity_control_does_not_use_test_labels(self):
        train=pd.DataFrame({'phage':['p','p','q'],'host':['h','i','h'],'label':[1,0,0]})
        test=pd.DataFrame({'phage':['p','z'],'host':['h','z'],'label':[1,1]})
        a=fractions(train,test);test.label=0;b=fractions(train,test)
        for k in a:np.testing.assert_array_equal(a[k],b[k])
        self.assertEqual(a['phage_positive_fraction'][1],1/3)
    def test_scaler_fit_is_local(self):
        X=np.arange(90,dtype=float).reshape(30,3);y=np.tile([0,1],15)
        m=estimator('GBM',{'max_depth':2,'n_estimators':2},7)
        m.fit(X[:20],y[:20]);np.testing.assert_allclose(m['scale'].mean_,X[:20].mean(0))
        self.assertFalse(np.allclose(m['scale'].mean_,X.mean(0)))
    def test_mlp_has_no_internal_row_holdout(self):
        m=estimator('MLP',{'hidden_layer_sizes':[64],'alpha':.001},7)
        self.assertFalse(m['model'].early_stopping)
        self.assertGreater(m['model'].n_iter_no_change,m['model'].max_iter)
        self.assertEqual(m['model'].tol,0)
    def test_assert_split_rejects_axis_overlap(self):
        df=pd.DataFrame({'host_cluster':[1,1]})
        with self.assertRaises(AssertionError):assert_split(df,[0],[1],['host'])

def main(write=True, include_toy_fit=True, portable_float32=False):
    if portable_float32 and write:raise ValueError('Portable inference must not regenerate scientific reports')
    started=time.monotonic()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(LeakageGuards)
    if not include_toy_fit:
        suite=unittest.TestSuite(test for test in suite if test._testMethodName!='test_scaler_fit_is_local')
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    assert result.wasSuccessful()
    checks={'unit_tests':result.testsRun,'artifact_reload_predictions':0,'split_checks':0,'input_hash_checks':0}
    for analysis in ['secondsource','nested_comparison']:
        dest=OUT/analysis
        completion=json.loads((dest/'completion.json').read_text())
        assert completion['complete']
        plan=json.loads((dest/'retrospective_protocol.json').read_text())
        assert sha(dest/'retrospective_protocol.json')==completion['protocol_sha256']
        for path,digest in plan['inputs'].items():
            assert sha(ROOT/path)==digest;checks['input_hash_checks']+=1
        df=pd.read_csv(dest/'input_rows.csv');X=np.load(dest/'input_features.npz')['X']
        # Round-trip parsing preserves ties exactly; pandas' default fast
        # parser can separate equal fractions by one ULP and change AP/AUROC.
        pred=pd.read_csv(dest/'predictions.csv',float_precision='round_trip');splits=json.loads((dest/'splits.json').read_text())
        pooled=pd.read_csv(dest/'pooled_metrics.csv',float_precision='round_trip')
        if analysis=='secondsource':
            frozen=joblib.load(ROOT/'data/results_v2/external_validation_models/full.joblib')
            old=pd.read_csv(ROOT/'data/results_v2/external_staph_predictions.csv')
            np.testing.assert_array_equal(df[['phage','host','label']],old[['phage','host','label']])
            reproduced=frozen['model'].predict_proba(frozen['scaler'].transform(X))[:,1]
            np.testing.assert_allclose(reproduced,old.prob_full,rtol=1e-6,atol=1e-7)
            checks['staph_assembly_matches_frozen_prediction_max_error']=float(np.max(np.abs(reproduced-old.prob_full)))
            assert len(splits)==48
            counts={};testrows={}
            for split in splits:
                tr,te=np.array(split['train']),np.array(split['test']);regime=split['regime']
                axes=[] if regime=='shared_entities_allowed' else ['host','phage'] if regime=='both_clusters_excluded' else [regime.split('_')[0]]
                assert_split(df,tr,te,axes);checks['split_checks']+=1
                key=(split['fold'],split['seed']);c=tuple(np.bincount(df.iloc[tr].label,minlength=2))
                if key in counts:assert counts[key]==c and testrows[key]==te.tolist()
                else:counts[key]=c;testrows[key]=te.tolist()
                for name,cols in [('phage_only',slice(0,212)),('pair_only',slice(424,448)),('full',slice(0,448))]:
                    model=joblib.load(dest/'fitted_models'/f'{split["fold"]}_{split["seed"]}_{regime}_{name}.joblib')
                    p=model.predict_proba(X[te,cols])[:,1]
                    saved=pred.loc[pred.fold.eq(split['fold'])&pred.seed.eq(split['seed'])&pred.regime.eq(regime)&pred.method.eq(name)]
                    np.testing.assert_array_equal(saved.row_index,te)
                    assert_saved_predictions(saved.score,p,model,portable_float32=portable_float32)
                    checks['artifact_reload_predictions']+=1
                for name,p in fractions(df.iloc[tr],df.iloc[te]).items():
                    saved=pred.loc[pred.fold.eq(split['fold'])&pred.seed.eq(split['seed'])&pred.regime.eq(regime)&pred.method.eq(name)]
                    np.testing.assert_allclose(saved.score,p,rtol=1e-12,atol=1e-12)
            # Report every originally assigned block, including excluded block.
            rng=np.random.default_rng(42)
            pb={c:i%5 for i,c in enumerate(rng.permutation(df.phage_cluster.unique()))}
            hb={c:i%5 for i,c in enumerate(rng.permutation(df.host_cluster.unique()))}
            blocks=[]
            for i in range(5):
                test=df.loc[df.phage_cluster.map(pb).eq(i)&df.host_cluster.map(hb).eq(i)]
                pos=int(test.label.sum());neg=len(test)-pos
                blocks.append({'assigned_block':i,'n':len(test),'positive':pos,'negative':neg,'eligible':pos>=3 and neg>=3})
            if write:save_json(dest/'all_block_eligibility.json',blocks)
            else:assert blocks==json.loads((dest/'all_block_eligibility.json').read_text())
        else:
            assert len(splits)==10
            sels=json.loads((dest/'selections.json').read_text())
            inner_metrics=pd.read_csv(dest/'inner_metrics.csv',float_precision='round_trip')
            assert len(inner_metrics)==180
            for split in splits:
                tr,te=np.array(split['train']),np.array(split['test']);regime=split['regime']
                axes=['host'] if regime=='host_cluster' else ['phage','host']
                assert_split(df,tr,te,axes);checks['split_checks']+=1
                for inner in split['inner']:
                    it,iv=inner['train'],inner['validation'];assert_split(df,it,iv,axes)
                    assert set(it).issubset(tr) and set(iv).issubset(tr)
                    assert not np.intersect1d(te,np.r_[it,iv]).size;checks['split_checks']+=1
                for family in ['GBM','MLP']:
                    sel=next(v for v in sels if v['regime']==regime and v['fold']==split['fold'] and v['family']==family)
                    assert sel['chosen']==int(np.argmax(sel['mean_inner_aurocs']))
                    subset=inner_metrics.loc[inner_metrics.regime.eq(regime)&inner_metrics.outer_fold.eq(split['fold'])&inner_metrics.family.eq(family)]
                    means=subset.groupby('candidate',sort=True).auroc.mean().to_numpy()
                    np.testing.assert_allclose(means,sel['mean_inner_aurocs'],rtol=1e-12,atol=1e-12)
                    assert sel['parameters']==CANDIDATES[family][sel['chosen']]
                    model=joblib.load(dest/'fitted_models'/f'{regime}_{split["fold"]}_{family}.joblib')
                    assert model.n_features_in_==448
                    if family=='MLP':assert model['model'].n_iter_==150
                    np.testing.assert_allclose(model['scale'].mean_,X[tr].mean(0,dtype=np.float64),rtol=1e-10,atol=1e-10)
                    p=model.predict_proba(X[te])[:,1]
                    saved=pred.loc[pred.fold.eq(split['fold'])&pred.regime.eq(regime)&pred.method.eq(family)]
                    np.testing.assert_array_equal(saved.row_index,te)
                    assert_saved_predictions(saved.score,p,model,portable_float32=portable_float32);checks['artifact_reload_predictions']+=1
        group=['regime','method','seed'] if analysis=='secondsource' else ['regime','method']
        for keys,p in pred.groupby(group):
            mask=np.ones(len(pooled),bool)
            for k,v in zip(group,keys):mask &= pooled[k].eq(v)
            row=pooled.loc[mask].iloc[0]
            for key,value in metric(p.label,p.score).items():np.testing.assert_allclose(row[key],value,rtol=1e-6,atol=1e-7)
        checks[analysis]={'prediction_rows':len(pred),'split_records':len(splits),'completion':completion}
    checks['elapsed_seconds']=time.monotonic()-started;checks['passed']=True
    if write:save_json(OUT/'model_verification.json',checks)
    checks['verification_files_written']=2 if write else 0
    checks['toy_scaler_guard_fit_executed']=include_toy_fit
    checks['analyzed_models_refitted']=0
    checks['portable_float32_inference_requested']=portable_float32
    print(json.dumps(checks,indent=2))
    return checks

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--read-only',action='store_true',help='Do not rewrite saved reports or eligibility tables')
    parser.add_argument('--skip-toy-fit',action='store_true',help='Skip only the two-tree synthetic scaler guard; reload and check every analyzed model')
    parser.add_argument('--portable-float32',action='store_true',help='Explicit float32 MLP inference-only comparison; use verify_extension_analysis.py with this flag for the mandatory pooled/decision diagnostic guards')
    args=parser.parse_args()
    sys.dont_write_bytecode=args.read_only
    if args.portable_float32 and not args.read_only:parser.error('Portable inference must not regenerate scientific reports')
    main(write=not args.read_only,include_toy_fit=not args.skip_toy_fit,portable_float32=args.portable_float32)
