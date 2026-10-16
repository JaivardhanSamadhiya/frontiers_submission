"""Retrospective replication of a fixed selector family on published banks."""
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import re
import zipfile
import numpy as np
import pandas as pd
from Bio import Phylo
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'data/external_validation/genophi_2026'
OUT=ROOT/'data/results_v2/additional_banks'
spec=importlib.util.spec_from_file_location('nested',ROOT/'experiments/36_nested_corrected_selector.py')
nested=importlib.util.module_from_spec(spec); spec.loader.exec_module(nested)
GRID=[5,10,20,40]
SEED=20261005
BANKS={'Klebsiella-1':'klebsiella1','Klebsiella-2':'klebsiella2',
       'Pseudomonas':'pseudomonas','Vibrionaceae':None}


def accession_alias(name):
    # Only explicit assembly accession suffixes; no arbitrary name shortening.
    return re.sub(r'^(GC[AF]_\d+)\.\d+$',r'\1',str(name))


def tree_distances(tree,names):
    tips=tree.get_terminals()
    mapping={}
    for tip in tips:
        key=accession_alias(tip.name)
        if key in mapping: raise ValueError(f'Ambiguous tree alias {key}')
        mapping[key]=tip
    missing=[n for n in names if accession_alias(n) not in mapping]
    if missing: raise ValueError(f'Missing tree tips: {missing}')
    if any(c.branch_length is not None and c.branch_length<0 for c in tree.find_clades()):
        raise ValueError('Negative branch length')
    dd=np.zeros((len(names),len(names)))
    for i,name in enumerate(names):
        for j in range(i):
            dd[i,j]=dd[j,i]=tree.distance(mapping[accession_alias(name)],mapping[accession_alias(names[j])])
    assert np.isfinite(dd).all()
    return dd,{n:mapping[accession_alias(n)].name for n in names}


def triplet(scores):
    return np.argsort(-np.asarray(scores),kind='stable')[:3]


def greedy_coverage(y,observed):
    covered=np.zeros(len(y),dtype=bool); selected=[]
    positive=(y>0)&observed
    breadth=np.divide(positive.sum(axis=0),observed.sum(axis=0),out=np.zeros(y.shape[1]),where=observed.sum(axis=0)>0)
    for _ in range(3):
        candidates=[j for j in range(y.shape[1]) if j not in selected]
        j=max(candidates,key=lambda c:(int((positive[:,c]&~covered).sum()),breadth[c],-c))
        selected.append(j); covered|=positive[:,j]
    return np.asarray(selected)


def tune(train,y,observed,groups,distance):
    values={k:[] for k in GRID}
    for it,iv in GroupKFold(n_splits=min(3,len(np.unique(groups[train])))).split(train,groups=groups[train]):
        tr,te=train[it],train[iv]
        for k in GRID:
            scores=nested.nearest_scores(tr,te,y,observed,distance,k)
            for i,h in enumerate(te):
                lo,up,frac,complete=nested.outcome(h,triplet(scores[i]),y,observed)
                values[k].append((lo,frac))
    metrics={k:np.mean(v,axis=0).tolist() for k,v in values.items()}
    return max(GRID,key=lambda k:(*metrics[k],-k)),metrics


def bootstrap(frame,baseline):
    unique,idx=np.unique(frame.group,return_inverse=True)
    values=frame[['design_lower','design_upper',f'{baseline}_lower',f'{baseline}_upper']].to_numpy()
    rng=np.random.default_rng(SEED); draws=[]
    for _ in range(2000):
        count=np.bincount(rng.integers(len(unique),size=len(unique)),minlength=len(unique))
        m=np.average(values,axis=0,weights=count[idx])
        draws.append([m[0]-m[2],m[0]-m[3],m[1]-m[2]])
    a=np.asarray(draws)
    return {'difference':float(values[:,0].mean()-values[:,2].mean()),
        'ci':np.quantile(a[:,0],[.025,.975]).tolist(),
        'missingness_robust_interval':[float(np.quantile(a[:,1],.025)),float(np.quantile(a[:,2],.975))]}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    z=zipfile.ZipFile(RAW/'processed_data.zip')
    data=pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'),dtype={'strain':str,'phage':str})
    assert not data.duplicated(['dataset','strain','phage']).any()
    assert set(data.interaction.dropna().unique()) <= {0,1}
    tz=zipfile.ZipFile(RAW/'published_phylogenies.zip')
    records,pair_records,summaries,inventory,settings=[],[],[],[],[]
    for bank,tree_name in BANKS.items():
        sub=data[data.dataset==bank]
        matrix=sub.pivot(index='strain',columns='phage',values='interaction').sort_index().sort_index(axis=1)
        names=matrix.index.tolist(); phages=matrix.columns.tolist()
        inv={'bank':bank,'hosts':len(names),'phages':len(phages),'assayed_pairs':int(matrix.notna().sum().sum()),
             'positive_pairs':int(matrix.fillna(0).sum().sum()),'missing_pairs':int(matrix.isna().sum().sum())}
        if tree_name is None:
            inv['status']='not evaluated: no published host phylogeny supplied in source archive'
            inventory.append(inv); continue
        tree=Phylo.read(io.StringIO(tz.read(f'SuppData1_newick_strain_trees/{tree_name}.newick').decode()),'newick')
        try: distance,aliases=tree_distances(tree,names)
        except ValueError as e:
            inv['status']=str(e); inventory.append(inv); continue
        inv['status']='resolved'; inv['tree_tip_mapping']=aliases; inventory.append(inv)
        np.savez_compressed(OUT/f'{tree_name}_patristic_distances.npz',names=names,distance=distance)
        y=matrix.fillna(0).to_numpy(int); observed=matrix.notna().to_numpy()
        for threshold in [None,1e-4,1e-3,1e-2]:
            group_label='identity' if threshold is None else str(threshold)
            groups=np.arange(len(names)) if threshold is None else connected_components(csr_matrix(distance<=threshold),directed=False)[1]
            n_groups=len(np.unique(groups))
            if n_groups<5:
                summaries.append({'bank':bank,'grouping':group_label,'groups':n_groups,'status':'skipped: fewer than five groups'});continue
            local_rows,local_pairs=[],[]
            for fold,(tr,te) in enumerate(GroupKFold(5).split(y,groups=groups)):
                assert set(groups[tr]).isdisjoint(groups[te])
                if len(np.unique(groups[tr]))<2: raise RuntimeError('Insufficient inner training groups')
                k,inner=tune(tr,y,observed,groups,distance)
                scores=nested.nearest_scores(tr,te,y,observed,distance,k)
                prevalence=np.divide((y[tr]*observed[tr]).sum(axis=0),observed[tr].sum(axis=0),
                    out=np.zeros(y.shape[1]),where=observed[tr].sum(axis=0)>0)
                breadth,greedy=triplet(prevalence),greedy_coverage(y[tr],observed[tr])
                settings.append({'bank':bank,'grouping':group_label,'fold':fold,'k':k,'inner':inner})
                for i,h in enumerate(te):
                    row={'bank':bank,'grouping':group_label,'fold':fold,'group':int(groups[h]),'host':names[h],
                         'k':k,'nearest_training_distance':float(distance[h,tr].min())}
                    for method,chosen in [('design',triplet(scores[i])),('breadth',breadth),('greedy',greedy)]:
                        lo,up,frac,complete=nested.outcome(h,chosen,y,observed)
                        row.update({f'{method}_phages':'|'.join(phages[j] for j in chosen),
                            f'{method}_lower':lo,f'{method}_upper':up,f'{method}_complete':complete})
                    local_rows.append(row)
                    for j in np.flatnonzero(observed[h]):
                        local_pairs.append({'bank':bank,'grouping':group_label,'fold':fold,'host':names[h],
                            'phage':phages[j],'label':int(y[h,j]),'score':float(scores[i,j]),'breadth_score':float(prevalence[j])})
            frame=pd.DataFrame(local_rows); pairs=pd.DataFrame(local_pairs)
            assert len(frame)==len(names) and not frame.host.duplicated().any()
            summary={'bank':bank,'grouping':group_label,'groups':n_groups,'hosts':len(frame),'status':'completed',
                **{f'{method}_coverage':float(frame[f'{method}_lower'].mean()) for method in ['design','breadth','greedy']},
                'design_upper':float(frame.design_upper.mean()),'vs_breadth':bootstrap(frame,'breadth'),
                'vs_greedy':bootstrap(frame,'greedy'),'pair_auroc':float(roc_auc_score(pairs.label,pairs.score)),
                'pair_ap':float(average_precision_score(pairs.label,pairs.score)),
                'pair_positive_prevalence':float(pairs.label.mean()),
                'breadth_pair_auroc':float(roc_auc_score(pairs.label,pairs.breadth_score)),
                'breadth_pair_ap':float(average_precision_score(pairs.label,pairs.breadth_score))}
            summaries.append(summary); records.extend(local_rows); pair_records.extend(local_pairs)
            print(json.dumps(summary),flush=True)
    pd.DataFrame(records).to_csv(OUT/'host_selections.csv',index=False)
    pd.DataFrame(pair_records).to_csv(OUT/'pair_predictions.csv',index=False)
    report={'status':'retrospective replication with bank-specific training; not frozen cross-source transport',
        'outcome':'observed individual-constituent union coverage; no mixtures tested',
        'selection':'unconstrained three-phage score ranking; source phage diversity metadata not assumed',
        'inventory':inventory,'summaries':summaries,'settings':settings,
        'sources':{'matrices':'https://doi.org/10.5281/zenodo.21696332','trees':'https://doi.org/10.1038/s41564-026-02482-5',
                   'license':'CC BY 4.0'},
        'checksums':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [RAW/'processed_data.zip',RAW/'published_phylogenies.zip',ROOT/'SECOND_BANK_AND_CANONICAL_PROTOCOL.md']}}
    (OUT/'summary.json').write_text(json.dumps(report,indent=2))


if __name__=='__main__': main()
