"""Retrospective reporting audit and train-only degree controls for Frontiers.

No saved model or outcome is replaced. No parameter is selected on test results.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from precisionphage.splits import leave_one_group_out, combined_unseen_folds
OUT = ROOT/'data/results_v2/frontiers_revision_2026_10_14/reporting_audit'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(z):
    return {'pairs':len(z),'positive':int(z.label.sum()),'prevalence':float(z.label.mean()) if len(z) else None,
            'phages':z.phage.nunique(),'host_identifiers':z.host.nunique(),
            'host_species_identifiers':z.host_species.nunique(),
            **({k+'_count':z[k].nunique() for k in ['phage_cluster','host_cluster']} if 'phage_cluster' in z else {})}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    full_path = ROOT/'data/interim_v2/interactions.csv'
    model_path = ROOT/'data/interim_v2/interactions_modelable.csv'
    cluster_path = ROOT/'data/cache_v2/clusters_k21_d0.05.json'
    full = pd.read_csv(full_path); cov = pd.read_csv(model_path)
    assert len(full)==8849 and len(cov)==1947
    clusters = json.loads(cluster_path.read_text())
    for axis in ['phage','host']:
        cov[axis+'_cluster'] = cov[axis].map(clusters[axis])
    assert not cov[['phage_cluster','host_cluster']].isna().any().any()
    keys = ['phage','host','study']
    coverage = full.merge(cov[keys].assign(covered=True),on=keys,how='left',validate='one_to_one')
    coverage['covered'] = coverage.covered.eq(True)
    coverage.to_csv(OUT/'coverage_records.csv',index=False)
    coverage_rows=[]
    for dimension in ['study','host_genus','host_species','label']:
        for value,z in coverage.groupby(dimension,dropna=False):
            for state,p in z.groupby('covered'):
                coverage_rows.append({'dimension':dimension,'value':str(value),'covered':bool(state),**summarize(p)})
    pd.DataFrame(coverage_rows).to_csv(OUT/'coverage_strata.csv',index=False)
    design = {
        'loso_species':list(leave_one_group_out(cov,'host_species','loso',3,3)),
        'logo_genus':list(leave_one_group_out(cov,'host_genus','logo',3,3)),
        'host_cluster':list(leave_one_group_out(cov,'host_cluster','host_cluster',3,3)),
        'phage_cluster':list(leave_one_group_out(cov,'phage_cluster','phage_cluster',3,3)),
        'combined_unseen':list(combined_unseen_folds(cov,'phage_cluster','host_cluster',5,42,3,3)),
    }
    old = np.load(ROOT/'data/results_v2/significance_pooled_preds.npz',allow_pickle=False)
    folds_report, cohorts, degree_rows, fold_metrics, verified_gbm = [],[],[],[],[]
    for regime,folds in design.items():
        indices=np.concatenate([f.test_idx for f in folds])
        assert len(np.unique(indices))==len(indices)
        if regime+'__y' in old:
            np.testing.assert_array_equal(cov.iloc[indices].label.to_numpy(),old[regime+'__y'])
        cohorts.append({'regime':regime,'folds':len(folds),**summarize(cov.iloc[indices]),
                        'not_tested_pairs':len(cov)-len(indices)})
        offset=0
        for fold in folds:
            tr,te=cov.iloc[fold.train_idx],cov.iloc[fold.test_idx]
            if regime+'__GBM' in old:
                score=old[regime+'__GBM'][offset:offset+len(te)]
                verified_gbm.append({'regime':regime,'fold':fold.name,'n':len(te),
                                     'positive':int(te.label.sum()),'auroc':roc_auc_score(te.label,score),
                                     'average_precision':average_precision_score(te.label,score)})
            offset+=len(te)
            folds_report.append({'regime':regime,'fold':fold.name,
                **{'train_'+k:v for k,v in summarize(tr).items()},
                **{'test_'+k:v for k,v in summarize(te).items()},
                'discarded_in_fold':len(cov)-len(tr)-len(te)})
            pos=tr[tr.label==1]
            pdg=pos.groupby('phage').size(); hdg=pos.groupby('host').size()
            p=te.phage.map(pdg).fillna(0).to_numpy(float)
            h=te.host.map(hdg).fillna(0).to_numpy(float)
            # Fixed normalization by maxima in training, then equal weighting.
            additive=(p/max(1,pdg.max())+h/max(1,hdg.max()))/2
            for method,score in [('phage_positive_degree',p),('host_positive_degree',h),('additive_positive_degree',additive)]:
                rows=te[['phage','host','label']].copy()
                rows['row_index']=fold.test_idx;rows['regime']=regime;rows['fold']=fold.name
                rows['method']=method;rows['score']=score
                degree_rows.append(rows)
                fold_metrics.append({'regime':regime,'fold':fold.name,'method':method,
                    'pairs':len(te),'auroc':roc_auc_score(te.label,score),
                    'average_precision':average_precision_score(te.label,score)})
    pd.DataFrame(cohorts).to_csv(OUT/'evaluation_cohorts.csv',index=False)
    pd.DataFrame(folds_report).to_csv(OUT/'evaluation_fold_counts.csv',index=False)
    verified_gbm=pd.DataFrame(verified_gbm)
    verified_gbm.to_csv(OUT/'verified_gbm_fold_metrics.csv',index=False)
    verified_gbm.groupby('regime').auroc.agg(['mean','std','min','max','size']).reset_index().to_csv(OUT/'verified_gbm_fold_summary.csv',index=False)
    degree=pd.concat(degree_rows,ignore_index=True)
    degree.to_csv(OUT/'degree_predictions.csv',index=False)
    pd.DataFrame(fold_metrics).to_csv(OUT/'degree_fold_metrics.csv',index=False)
    degree_metrics=[]
    for (regime,method),p in degree.groupby(['regime','method']):
        row={'regime':regime,'method':method,'n':len(p),'prevalence':p.label.mean(),
             'auroc':roc_auc_score(p.label,p.score),'average_precision':average_precision_score(p.label,p.score)}
        if regime+'__GBM' in old:
            row['saved_gbm_auroc']=roc_auc_score(old[regime+'__y'],old[regime+'__GBM'])
        degree_metrics.append(row)
    pd.DataFrame(degree_metrics).to_csv(OUT/'degree_summary.csv',index=False)
    feature_rows=[
        ('Canonical nucleotide 4-mers',136,272,'genomic.kmer_spectrum','Per-genome frequencies, reverse-complement equivalent'),
        ('Canonical dinucleotides',10,20,'genomic.kmer_spectrum','Per-genome frequencies'),
        ('Codon frequencies',64,128,'genomic.codon_usage','All three forward frames; no CDS annotation requirement'),
        ('GC fraction and log10 length plus one',2,4,'genomic._node_vector','Per-genome composition and length'),
        ('Local composition comparisons',None,4,'genomic.edge_features_from_spectra','Cosine distance, L1, Pearson, abundance Jaccard'),
        ('Source supplied VHIP features',None,4,'VirusHostInter.csv','k3dist, k6dist, GCdiff, Homology; not recomputed locally'),
        ('Exact nucleotide word proxies',None,8,'seqmatch.homology_features','Four descriptors at k16 and k20'),
        ('CRISPR-like spacer proxies',None,4,'seqmatch','Heuristic repeat/spacer matches; not verified immunity measurements'),
        ('Translated amino-acid word proxies',None,4,'proteins','Six-frame translation and exact amino-acid word sharing; not protein alignment'),
    ]
    pd.DataFrame(feature_rows,columns=['family','per_genome_dimensions','model_input_dimensions','source','calculation']).to_csv(OUT/'feature_schema.csv',index=False)
    assert sum(r[2] for r in feature_rows)==448
    # Preserve historical nested-tuning outcomes exactly as recorded, with schema caveat.
    tuning=[]
    for logfile in ['sig_tuned.log','sig_tuned2.log']:
        text=(ROOT/'data/results_v2'/logfile).read_text(encoding='utf-8')
        import re
        rows=re.findall(r'^\s*(loso_species|host_cluster|phage_cluster|combined_unseen)\s+GBM\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',text,flags=re.M)
        assert len(rows)==4
        for regime,auc,_,__,ece in rows:
            tuning.append({'historical_log':logfile,'regime':regime,'auroc_as_logged':float(auc),
                           'ece_as_logged':float(ece),'pair_descriptors_as_logged':20,
                           'scope':'Historical exploratory run, 424 node +20 pair descriptors; not the final 448-column run. Exact selected code/settings not recoverable from log alone.'})
    pd.DataFrame(tuning).to_csv(OUT/'historical_tuning_results.csv',index=False)
    audit={'created_utc':datetime.now(timezone.utc).isoformat(),'status':'retrospective audit; no new model fitting',
           'ncbi_coverage':{'full':summarize(full[full.study=='NCBI_HR']),
                            'covered':summarize(cov),'uncovered':summarize(coverage[(coverage.study=='NCBI_HR') & ~coverage.covered])},
           'genome_bundle_scope':'Coverage describes the preserved analysis bundle, not genome availability in all public databases.',
           'degree_definition':'Training-only positive row counts for exact phage/host IDs; unknown IDs score zero; additive combines each count divided by its training maximum with equal weights.',
           'degree_limitations':'No sequence similarity or tuning; descriptive control. Failure to match a sequence model does not rule out homology/base-rate shortcuts.',
           'uncertainty':'No independent-row p-values or fold confidence intervals calculated.',
           'fold_metric_correction':'The old coldstart_fold_aucs.json does not reproduce the saved pooled predictions. The revision derives every primary fold AUROC from significance_pooled_preds.npz and the verified fold order.',
           'sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),full_path,model_path,cluster_path,ROOT/'data/results_v2/significance_pooled_preds.npz',ROOT/'data/results_v2/sig_tuned.log',ROOT/'data/results_v2/sig_tuned2.log']}}
    (OUT/'audit.json').write_text(json.dumps(audit,indent=2))
    print(pd.DataFrame(cohorts).to_string(index=False));print(pd.DataFrame(degree_metrics).to_string(index=False))


if __name__=='__main__':main()
