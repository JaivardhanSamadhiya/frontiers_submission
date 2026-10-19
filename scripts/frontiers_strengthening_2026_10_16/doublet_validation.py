"""Retrospective measured same-size menu analysis; see DOUBLET_PROTOCOL.md."""
from __future__ import annotations
import hashlib
import importlib.util
import json
from pathlib import Path
import platform

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'data/results_v2/frontiers_strengthening_2026_10_16/doublets'
UP = ROOT / 'external/gaborieau_coli_phage_interactions_2023'
FROZEN = ROOT / 'data/results_v2/external_cocktail_validation'
OLD = ROOT / 'data/results_v2/frontiers_revision_2026_10_14/cocktail_controls'
AUDIT = ROOT / 'data/results_v2/post_rejection_audit'
PRIMARY = 'repaired_product'
REPS, SEED, TOL = 2000, 20261016, 1e-12
INPUTS = {
    'protocol': Path(__file__).with_name('DOUBLET_PROTOCOL.md'),
    'script': Path(__file__),
    'compositions': UP / 'dev/cocktails/data/cocktails_composition.csv',
    'outcomes': UP / 'dev/cocktails/data/cml_cocktails.csv',
    'pairs': FROZEN / 'frozen_pair_predictions.csv',
    'manifest': FROZEN / 'frozen_prediction_manifest.json',
    'reference': AUDIT / 'development_rank_reference.npz',
    'repaired_saved': AUDIT / 'batch_invariant_cocktail_diagnostic.csv',
    'training_matrix': UP / 'data/interactions/interaction_matrix.csv',
    'groups': OLD / 'host_group_assignments.csv',
    'selected_triplets': OLD / 'selected_triplet_scores.csv',
}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def percentiles(raw, ref):
    if raw.shape[1] != ref.shape[1]:
        raise ValueError('Mismatched phage order')
    result = np.empty_like(raw)
    for j in range(raw.shape[1]):
        v = np.sort(ref[:, j])
        result[:, j] = (np.searchsorted(v, raw[:, j], side='left') +
                         np.searchsorted(v, raw[:, j], side='right') + 1) / (2 * (len(v) + 1))
    return result

def aggregate(x, method):
    if not np.isfinite(x).all() or ((x < 0) | (x > 1)).any():
        raise ValueError('Scores must be finite and in [0,1]')
    return {'product': lambda: 1 - np.prod(1-x, axis=-1),
            'mean': lambda: x.mean(axis=-1), 'max': lambda: x.max(axis=-1)}[method]()

def choose_expected(scores, y):
    tied = np.isclose(scores, scores.max(axis=1, keepdims=True), atol=TOL, rtol=0)
    return (tied * y).sum(axis=1) / tied.sum(axis=1), tied

def pair_concordance(y, score, w=None):
    """Return weighted concordant-pair numerator and total positive-negative mass."""
    y, score = np.asarray(y), np.asarray(score)
    w = np.ones(len(y)) if w is None else np.asarray(w)
    p, n = y == 1, y == 0
    diff = score[p, None] - score[None, n]
    mass = w[p, None] * w[None, n]
    credit = (diff > TOL).astype(float) + .5 * (np.abs(diff) <= TOL)
    return float((mass * credit).sum()), float(mass.sum())

def bootstrap_weights(groups, compositions, crossed=False):
    _, gi = np.unique(groups, return_inverse=True)
    _, ci = np.unique(compositions, return_inverse=True)
    ng, nc = gi.max()+1, ci.max()+1
    rng = np.random.default_rng(SEED + int(crossed))
    weights = []
    for _ in range(REPS):
        wg = np.bincount(rng.integers(ng, size=ng), minlength=ng)[gi].astype(float)
        if crossed:
            wg *= np.bincount(rng.integers(nc, size=nc), minlength=nc)[ci]
        weights.append(wg)
    return np.asarray(weights)

def ci(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    return (*np.quantile(values, [.025, .975]), len(values)) if len(values) else (np.nan, np.nan, 0)

def markdown_table(frame):
    """Small report table without an optional formatting dependency."""
    cols = frame.columns.tolist()
    rows = ['| ' + ' | '.join(cols) + ' |', '| ' + ' | '.join(['---']*len(cols)) + ' |']
    for row in frame.itertuples(index=False, name=None):
        rows.append('| ' + ' | '.join(f'{x:.6f}' if isinstance(x,float) else str(x) for x in row) + ' |')
    return '\n'.join(rows)

def weighted_means(w, values):
    return np.divide(w @ values, w.sum(axis=1), out=np.full(len(w), np.nan), where=w.sum(axis=1)>0)

def conditional_metrics(y, s, compositions, weights):
    micro_n = micro_d = 0.
    macro = []
    for key in np.unique(compositions):
        mask = compositions == key
        numerator, denominator = pair_concordance(y[mask], s[mask], weights[mask])
        if denominator > 0:
            macro.append(numerator/denominator)
            micro_n += numerator
            micro_d += denominator
    return np.array([micro_n/micro_d if micro_d else np.nan,
                     np.mean(macro) if macro else np.nan])

def load():
    composition = pd.read_csv(INPUTS['compositions'], sep=';').replace({'bacteria': {'BCH857': 'BCH857-72'}})
    outcome = pd.read_csv(INPUTS['outcomes'], sep=';').replace({'bacteria': {'BCH857': 'BCH857-72'}})
    triplet = pd.read_csv(INPUTS['selected_triplets'])
    hosts = triplet.bacteria.to_list()
    assert len(hosts) == len(set(hosts)) == 100
    assert composition.bacteria.is_unique and outcome.bacteria.is_unique
    assert set(hosts) == set(composition.bacteria) == set(outcome.bacteria)
    composition = composition.set_index('bacteria').loc[hosts]
    outcome = outcome.set_index('bacteria').loc[hosts]
    assert not outcome[['p1+p2', 'p1+p3', 'p2+p3']].isna().any().any()
    groups = pd.read_csv(INPUTS['groups']).set_index('bacteria').loc[hosts].host_group.to_numpy()
    assert len(np.unique(groups)) == 81
    pair = pd.read_csv(INPUTS['pairs'])
    assert len(pair) == 9600 and not pair.duplicated(['bacteria','phage']).any()
    manifest = json.loads(INPUTS['manifest'].read_text())
    assert sha(INPUTS['pairs']) == manifest['checksums']['frozen_pair_predictions']
    ref = np.load(INPUTS['reference'], allow_pickle=False)
    phages = ref['phages'].astype(str).tolist()
    pivot = lambda col: pair.pivot(index='bacteria', columns='phage', values=col).loc[hosts, phages].to_numpy()
    rawlog, rawknn = pivot('prob_logistic_raw'), pivot('prob_knn_raw')
    matrices = {'historical': pivot('score_ensemble'), 'raw_logistic': rawlog, 'raw_knn': rawknn,
                'repaired': .5*percentiles(rawlog, ref['logistic'])+.5*percentiles(rawknn, ref['knn']),
                'training_prevalence': pivot('training_prevalence')}
    train = pd.read_csv(INPUTS['training_matrix'], sep=';', index_col=0)
    train_frac = (train.gt(0).sum()/train.notna().sum()).loc[phages].to_numpy()
    np.testing.assert_allclose(matrices['training_prevalence'], np.tile(train_frac, (100,1)), atol=1e-15, rtol=0)
    for i in range(100):
        np.testing.assert_array_equal(matrices['repaired'][i],
            (.5*percentiles(rawlog[i:i+1], ref['logistic'])+.5*percentiles(rawknn[i:i+1], ref['knn']))[0])
    selected = np.array([[phages.index(p) for p in composition.loc[h, ['p1','p2','p3']]] for h in hosts])
    assert all(len(set(row)) == 3 for row in selected)
    trip_arrays = {name: np.take_along_axis(m, selected, axis=1) for name,m in matrices.items()}
    np.testing.assert_allclose(aggregate(trip_arrays['historical'], 'product'), triplet.pp_historical_product, atol=2e-15,rtol=0)
    saved = pd.read_csv(INPUTS['repaired_saved']).set_index('bacteria').loc[hosts]
    np.testing.assert_allclose(aggregate(trip_arrays['repaired'], 'product'), saved.development_reference_score, atol=2e-15,rtol=0)
    np.testing.assert_array_equal((outcome['p1+p2+p3']>0).astype(int), triplet.observed_lytic)
    keys = np.array(['|'.join(sorted(composition.loc[h,['p1','p2','p3']])) for h in hosts])
    np.testing.assert_array_equal(keys, triplet.composition_key)
    scores = {}
    for source, aggs in [('repaired',['product','max','mean']), ('historical',['product']),
                         ('raw_logistic',['product']), ('raw_knn',['product']),
                         ('training_prevalence',['mean','product'])]:
        for ag in aggs:
            scores[f'{source}_{ag}'] = aggregate(trip_arrays[source][:, [[0,1],[0,2],[1,2]]], ag)
    return hosts, groups, keys, outcome, composition, triplet, scores

def run():
    OUT.mkdir(parents=True, exist_ok=True)
    hashes = {k: sha(p) for k,p in INPUTS.items()}
    hosts, groups, keys, outcomes, compositions, triplet, scores = load()
    y = (outcomes[['p1+p2','p1+p3','p2+p3']].to_numpy()>0).astype(int)
    random = y.mean(axis=1)
    oracle = y.max(axis=1)
    choices, rows = {}, []
    for method,s in scores.items():
        choices[method], ties = choose_expected(s,y)
        for i,h in enumerate(hosts):
            for j,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
                rows.append({'bacteria':h, 'host_group':int(groups[i]), 'parent_composition':keys[i],
                    'option':['p1+p2','p1+p3','p2+p3'][j], 'phage1':compositions.loc[h,['p1','p2','p3']].iloc[a],
                    'phage2':compositions.loc[h,['p1','p2','p3']].iloc[b], 'assay_value':int(outcomes.iloc[i][['p1+p2','p1+p3','p2+p3'][j]]),
                    'active':int(y[i,j]),'method':method,'score':s[i,j],
                    'selection_probability': float(ties[i,j]/ties[i].sum())})
    pd.DataFrame(rows).to_csv(OUT/'doublet_option_scores.csv',index=False)
    hostframe = pd.DataFrame({'bacteria':hosts,'host_group':groups,'composition':keys,
        'random_choice_activity':random,'oracle_activity':oracle, **choices})
    hostframe.to_csv(OUT/'doublet_host_decisions.csv',index=False)
    boot = {kind:bootstrap_weights(groups,keys,kind=='crossed') for kind in ['host_cluster','crossed']}
    np.savez_compressed(OUT/'bootstrap_weights.npz',**boot)
    summaries=[]
    base = choices['training_prevalence_mean']
    for method, activity in {**choices,'uniform_choice':random,'observed_oracle':oracle}.items():
        for kind,w in boot.items():
            for contrast,values in [('activity',activity),('gain_vs_uniform',activity-random),('gain_vs_training_mean',activity-base)]:
                lo,hi,valid=ci(weighted_means(w,values))
                summaries.append({'method':method,'resampling':kind,'metric':contrast,'estimate':values.mean(),
                    'ci_low':lo,'ci_high':hi,'valid_replicates':valid,'total_replicates':REPS})
    summary=pd.DataFrame(summaries)
    summary.to_csv(OUT/'doublet_decision_summary.csv',index=False)
    ranking=[]
    for method,s in scores.items():
        nd=np.array([pair_concordance(y[i],s[i]) for i in range(100)])
        informative=nd[:,1]>0
        aucs=nd[informative,0]/nd[informative,1]
        for kind,w in boot.items():
            pairden=w@nd[:,1]
            bs=np.divide(w@nd[:,0],pairden,out=np.full(REPS,np.nan),where=pairden>0)
            hm=weighted_means(w[:,informative],aucs)
            for metric,estimate,values in [('pair_weighted_concordance',nd[:,0].sum()/nd[:,1].sum(),bs),
                                            ('host_mean_concordance',aucs.mean(),hm)]:
                lo,hi,valid=ci(values)
                ranking.append({'method':method,'resampling':kind,'metric':metric,'estimate':estimate,
                    'ci_low':lo,'ci_high':hi,'valid_replicates':valid,'informative_hosts':int(informative.sum()),
                    'positive_negative_pairs':int(nd[:,1].sum())})
    pd.DataFrame(ranking).to_csv(OUT/'doublet_conditional_ranking.csv',index=False)
    # Exact history-specific compositions, with no outcome-fitted adjustment.
    methods={'historical_product':'pp_historical_product','repaired_product':'pp_development_reference_product',
             'historical_max':'pp_historical_max','historical_mean':'pp_historical_mean',
             'raw_logistic_product':'pp_logistic_raw_product', 'raw_knn_product':'pp_knn_raw_product',
             'training_prevalence_mean':'bank_phage_prevalence_mean','host_phylogroup':'training_phylogroup_susceptibility'}
    ty=triplet.observed_lytic.to_numpy(int)
    bycomp=[]; conditional=[]; leaveout=[]
    counts=pd.Series(keys).value_counts().rename_axis('composition').reset_index(name='n').sort_values(['n','composition'],ascending=[False,True])
    for method,col in methods.items():
        s=triplet[col].to_numpy(float)
        for key in sorted(set(keys)):
            mask=keys==key
            num,den=pair_concordance(ty[mask],s[mask])
            bycomp.append({'method':method,'composition':key,'n':int(mask.sum()),'positive':int(ty[mask].sum()),
                'negative':int(mask.sum()-ty[mask].sum()),'positive_negative_pairs':int(den),'conditional_auroc':num/den if den else np.nan})
        point=conditional_metrics(ty,s,keys,np.ones(100))
        for kind,w in boot.items():
            values=np.array([conditional_metrics(ty,s,keys,wi) for wi in w])
            for j,name in enumerate(['pair_weighted_conditional_auroc','composition_mean_auroc']):
                lo,hi,valid=ci(values[:,j])
                conditional.append({'method':method,'metric':name,'resampling':kind,'estimate':point[j],
                                    'ci_low':lo,'ci_high':hi,'valid_replicates':valid,'total_replicates':REPS})
        exclusions={f'omit_composition_{i+1}':[k] for i,k in enumerate(counts.composition)}
        exclusions.update({'omit_largest_3':counts.composition.head(3).tolist(),'all_hosts':[]})
        for label,excluded in exclusions.items():
            mask=~np.isin(keys,excluded)
            both=np.unique(ty[mask]).size==2
            leaveout.append({'method':method,'sensitivity':label,'excluded_compositions':';'.join(excluded),
                'n':int(mask.sum()),'positive':int(ty[mask].sum()),'negative':int(mask.sum()-ty[mask].sum()),
                'auroc':roc_auc_score(ty[mask],s[mask]) if both else np.nan,
                'ap':average_precision_score(ty[mask],s[mask]) if both else np.nan,
                'prevalence':ty[mask].mean() if mask.any() else np.nan})
    pd.DataFrame(bycomp).to_csv(OUT/'triplet_composition_metrics.csv',index=False)
    cond=pd.DataFrame(conditional);cond.to_csv(OUT/'triplet_conditional_summary.csv',index=False)
    pd.DataFrame(leaveout).to_csv(OUT/'triplet_composition_sensitivities.csv',index=False)
    sans=np.array(hosts)!='AN03'
    sensitivities=[]
    for method,activity in choices.items():
        sensitivities.append({'analysis':'doublet_decision','method':method,'n':99,
            'metric':'activity','estimate':activity[sans].mean(),'gain_vs_uniform':(activity-random)[sans].mean()})
    for method,col in methods.items():
        s=triplet[col].to_numpy(float)
        vals=conditional_metrics(ty[sans],s[sans],keys[sans],np.ones(99))
        sensitivities.append({'analysis':'triplet_within_composition','method':method,'n':99,
            'metric':'pair_weighted_conditional_auroc','estimate':vals[0],'gain_vs_uniform':np.nan})
    pd.DataFrame(sensitivities).to_csv(OUT/'exclude_AN03_sensitivity.csv',index=False)
    report={'status':'completed retrospective extension; no predictor refitted','hosts':100,'host_clusters':81,
        'doublets':300,'doublet_positive':int(y.sum()),'doublet_prevalence':y.mean(),
        'binary_discordant_hosts':int((y.min(axis=1)!=y.max(axis=1)).sum()),
        'ordinal_discordant_hosts':int((outcomes[['p1+p2','p1+p3','p2+p3']].nunique(axis=1)>1).sum()),
        'parent_compositions':18,'informative_triplet_compositions':int(pd.DataFrame(bycomp).query("method == 'repaired_product'").positive_negative_pairs.gt(0).sum()),
        'uniform_choice_activity':random.mean(),'oracle_activity':oracle.mean(),
        'ties':{m:int((choose_expected(s,y)[1].sum(axis=1)>1).sum()) for m,s in scores.items()},
        'primary_decision':summary.query("method == 'repaired_product'").to_dict(orient='records'),
        'primary_conditional_triplet':cond.query("method == 'repaired_product'").to_dict(orient='records'),
        'upstream_commit':'99656d8c020fcd1dbbc30a02ee10274a07655562',
        'scope':'Retrospective prioritization among three measured two-phage subsets of externally selected triplets. No arbitrary cocktail design, unseen cohort, new organism or therapeutic efficacy validated.',
        'interval_scope':'Fixed models and supplied compositions; host-cluster primary and crossed parent-composition sensitivity; paired contrasts. Historical labels and scores were already accessible.',
        'shared_phage_limitation':'Options share constituents within hosts and phages recur across different compositions. Keeping each host menu intact and crossing parent compositions does not create independent phage banks or resample phage identities. Inference is conditional on this fixed phage bank; the intervals do not establish transfer to new phages or banks.',
        'singleton_batch_exact_hosts':100,'reconstructed_historical_and_repaired_triplet_scores':True,
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'sklearn':sklearn.__version__},
        'inputs':{k:{'path':str(p.relative_to(ROOT)),'sha256':hashes[k]} for k,p in INPUTS.items()},
        'input_hashes_unchanged':hashes=={k:sha(p) for k,p in INPUTS.items()}}
    assert report['input_hashes_unchanged']
    (OUT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    draw(summary,cond)
    primary=summary.query("method=='repaired_product' and resampling=='host_cluster'")
    text=['# Measured doublets and composition-conditioned host discrimination','',report['scope'],'',
          f"There are 300 measured doublets across 100 hosts, {int(y.sum())} active results, and 26 hosts with mixed binary outcomes. The uniform-choice expected activity is {random.mean():.3f}; the observed oracle bound is {oracle.mean():.3f}. Only {report['informative_triplet_compositions']} parent compositions contain both active and inactive triplets. Five compositions contain negatives, but two of those have only one host and no positive observations.",'',
          '## Primary doublet decision','',markdown_table(primary),'',
          'All methods and paired comparisons are in doublet_decision_summary.csv. Ties receive equal probability across tied maxima. A host-only phylogroup control is identical to uniform choice for this menu. Choosing among equal-size subsets avoids the structural superset advantage of product aggregation. Any favorable result remains conditional on the original study having chosen the three parent phages.','',
          '## Conditional host discrimination','',markdown_table(cond.query("method=='repaired_product'")),'',
          'The pair-weighted result compares active and inactive hosts only when they share exactly the same triplet. The composition-mean result weights each unique informative composition equally. Only three compositions identify these results, containing 72 hosts and 301 positive-negative host pairs; these 301 comparisons are not independent observations. Crossed intervals additionally multiply row weights by resampled parent-composition counts; the macro summary averages the unique informative compositions retained in each replicate. The repaired pair-weighted host-cluster interval excludes chance, but its crossed interval includes chance and is wide. Therefore, uncertainty across compositions remains substantial.','',
          report['shared_phage_limitation'],'',
          'The primary doublet gain intervals include zero against both uniform choice and the training susceptibility control. The experiment does not establish better choice of a measured doublet. Repaired product and mean select the same highest-scoring doublet here: with three constituent scores and equal-size pairs, both ordinarily choose the top two constituents. Thus their agreement is algebraic, not independent support for the aggregation formula.','',
          'Leave-one-composition-out and largest-three exclusions, AN03 sensitivity, all raw option scores, observed outcomes, selection probabilities and resampling weights are included. No ordinal efficacy optimization was performed. All input hashes were unchanged.']
    (OUT/'REPORT.md').write_text('\n'.join(text),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['inputs','versions']},indent=2))

def draw(summary,cond):
    plt.rcParams.update({'font.size':9,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(12,5.3),layout='constrained')
    selected=['repaired_product','historical_product','repaired_max','repaired_mean','raw_logistic_product','raw_knn_product','training_prevalence_mean']
    a=summary.query("resampling=='host_cluster' and metric=='gain_vs_uniform'").set_index('method').loc[selected]
    labels=['Repaired product','Historical product','Repaired maximum','Repaired mean','Raw logistic product','Raw neighbor product','Training susceptibility']
    axcross=summary.query("resampling=='crossed' and metric=='gain_vs_uniform'").set_index('method').loc[selected]
    axes[0].errorbar(axcross.estimate*100,np.arange(len(a))+.10,xerr=np.vstack([(axcross.estimate-axcross.ci_low)*100,(axcross.ci_high-axcross.estimate)*100]),fmt='s',capsize=3,color='#999999',label='Crossed host/composition')
    axes[0].errorbar(a.estimate*100,np.arange(len(a))-.10,xerr=np.vstack([(a.estimate-a.ci_low)*100,(a.ci_high-a.estimate)*100]),fmt='o',capsize=3,color='#21618c',label='Host cluster')
    axes[0].set(yticks=np.arange(len(a)),yticklabels=labels,xlabel='Activity gain over uniform choice (percentage points)',title='A. Choose one measured doublet per host')
    axes[0].invert_yaxis();axes[0].axvline(0,color='gray',ls='--')
    methods=['historical_product','repaired_product','historical_max','historical_mean','raw_logistic_product','raw_knn_product','training_prevalence_mean']
    b=cond.query("resampling=='host_cluster' and metric=='pair_weighted_conditional_auroc'").set_index('method').loc[methods]
    bcross=cond.query("resampling=='crossed' and metric=='pair_weighted_conditional_auroc'").set_index('method').loc[methods]
    axes[1].errorbar(bcross.estimate,np.arange(len(b))+.10,xerr=np.vstack([bcross.estimate-bcross.ci_low,bcross.ci_high-bcross.estimate]),fmt='s',capsize=3,color='#999999')
    axes[1].errorbar(b.estimate,np.arange(len(b))-.10,xerr=np.vstack([b.estimate-b.ci_low,b.ci_high-b.estimate]),fmt='o',capsize=3,color='#117864')
    axes[1].set(yticks=np.arange(len(b)),yticklabels=['Historical product','Repaired product','Historical maximum','Historical mean','Raw logistic product','Raw neighbor product','Training susceptibility'],xlabel='Within-composition concordance',title='B. Compare hosts given the same triplet',xlim=(-.02,1.02))
    axes[1].invert_yaxis();axes[1].axvline(.5,color='gray',ls='--')
    fig.legend(*axes[0].get_legend_handles_labels(), loc='outside lower center', ncol=2, fontsize=8)
    fig.suptitle('Retrospective fixed-menu analyses; percentile 95% intervals',fontsize=10)
    for ax in axes: ax.spines[['top','right']].set_visible(False)
    fig.savefig(OUT/'measured_menu_and_conditional_host.png',dpi=180)
    fig.savefig(OUT/'measured_menu_and_conditional_host.pdf')
    plt.close(fig)

if __name__=='__main__': run()
