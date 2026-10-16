"""Canonical sourmash grouping and fixed-test matched training sensitivity."""
import hashlib
import importlib.util
from importlib.metadata import version
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score, roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from precisionphage.data.genomes import GenomeIndex
from precisionphage.splits.canonical import sketch_fasta, distance_matrix
from precisionphage.splits.cluster import combined_unseen_folds

spec = importlib.util.spec_from_file_location('matched', ROOT / 'experiments/35_matched_cohort_leakage_audit.py')
matched = importlib.util.module_from_spec(spec)
spec.loader.exec_module(matched)
OUT = ROOT / 'data/results_v2/canonical_sensitivity'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df, xx = matched.load_design()
    historical_folds = list(combined_unseen_folds(df, 'phage_cluster', 'host_cluster', 5, 42))
    old = df.copy()
    summaries, hashes, cluster_rows = [], {}, []
    for axis, directory in [('phage', 'phages'), ('host', 'hosts')]:
        names = sorted(df[axis].unique())
        index = GenomeIndex([ROOT / 'external/phist_run' / directory])
        paths = [index.resolve(name) for name in names]
        assert all(path is not None for path in paths)
        for name, path in zip(names, paths):
            hashes[f'{axis}:{name}'] = hashlib.sha256(path.read_bytes()).hexdigest()
        cache = OUT / f'{axis}_canonical_2048.npz'
        hashfile = OUT / f'{axis}_fasta_checksums.json'
        axis_hashes = {name: hashes[f'{axis}:{name}'] for name in names}
        if cache.exists() and hashfile.exists() and json.loads(hashfile.read_text()) == axis_hashes:
            import sourmash
            z = np.load(cache)
            sketches = [sourmash.MinHash(n=2048, ksize=21, seed=42, mins=z[f's{i}']) for i in range(len(names))]
        else:
            with ThreadPoolExecutor(max_workers=2) as pool:
                sketches = []
                for i, mh in enumerate(pool.map(sketch_fasta, paths)):
                    sketches.append(mh)
                    if i % 50 == 0:
                        print(f'{axis}: sketched {i+1}/{len(names)}', flush=True)
            np.savez_compressed(cache, **{f's{i}': np.array(sorted(mh.hashes), dtype=np.uint64) for i,mh in enumerate(sketches)})
            hashfile.write_text(json.dumps(axis_hashes, indent=2))
        for size in [256, 2048]:
            small = [mh.downsample(num=size) for mh in sketches]
            dd = distance_matrix(small)
            np.savez_compressed(OUT / f'{axis}_distances_{size}.npz', names=names, distance=dd)
            count, groups = connected_components(csr_matrix(dd <= 0.05), directed=False)
            old_groups = old.drop_duplicates(axis).set_index(axis).loc[names, f'{axis}_cluster'].to_numpy()
            summaries.append({'axis': axis, 'sketch_size': size, 'entities': len(names),
                'historical_clusters': len(np.unique(old_groups)), 'canonical_clusters': count,
                'adjusted_rand_vs_historical': float(adjusted_rand_score(old_groups, groups)),
                'newly_joined_entity_pairs': int(np.triu((groups[:,None] == groups) & (old_groups[:,None] != old_groups),1).sum()),
                'newly_separated_entity_pairs': int(np.triu((groups[:,None] != groups) & (old_groups[:,None] == old_groups),1).sum())})
            cluster_rows.extend({'axis':axis,'entity':name,'sketch_size':size,'cluster':int(group)} for name,group in zip(names,groups))
            if size == 2048:
                df[f'{axis}_cluster'] = df[axis].map(dict(zip(names,groups)))
            print(summaries[-1],flush=True)
    pd.DataFrame(cluster_rows).to_csv(OUT / 'canonical_clusters.csv', index=False)
    pd.DataFrame(summaries).to_csv(OUT / 'cluster_comparison.csv', index=False)
    yy = df.label.to_numpy()
    preds, fold_rows = [], []
    regimes = ['shared_entities_allowed','host_cluster_excluded','phage_cluster_excluded','both_clusters_excluded']
    for fi, fold in enumerate(historical_folds):
        te = fold.test_idx
        ho = df.host_cluster.isin(df.iloc[te].host_cluster).to_numpy()
        po = df.phage_cluster.isin(df.iloc[te].phage_cluster).to_numpy()
        nt = ~np.isin(np.arange(len(df)),te)
        masks = [nt, nt & ~ho, nt & ~po, nt & ~ho & ~po]
        counts = {c:min(int((mask & (yy==c)).sum()) for mask in masks) for c in [0,1]}
        if not all(counts.values()):
            fold_rows.append({'fold':fi,'status':'not estimable: missing training class','counts':str(counts)})
            continue
        for seed in matched.SEEDS:
            for ri,(regime,mask) in enumerate(zip(regimes,masks)):
                rng=np.random.default_rng(seed+fi*100+ri)
                tr=np.concatenate([rng.choice(np.flatnonzero(mask & (yy==c)),counts[c],replace=False) for c in [0,1]])
                assert not np.intersect1d(tr,te).size
                if ri in [1,3]: assert set(df.iloc[tr].host_cluster).isdisjoint(df.iloc[te].host_cluster)
                if ri in [2,3]: assert set(df.iloc[tr].phage_cluster).isdisjoint(df.iloc[te].phage_cluster)
                scaler=StandardScaler().fit(xx[tr])
                model=XGBClassifier(n_estimators=400,max_depth=5,learning_rate=.05,subsample=.8,
                    colsample_bytree=.8,random_state=seed,tree_method='hist',n_jobs=2,eval_metric='logloss',verbosity=0)
                model.fit(scaler.transform(xx[tr]),yy[tr])
                score=model.predict_proba(scaler.transform(xx[te]))[:,1]
                part=df.iloc[te][['phage','host','label','phage_cluster','host_cluster']].copy()
                part['row_index'],part['fold'],part['seed'],part['regime'],part['score']=te,fi,seed,regime,score
                preds.append(part)
                fold_rows.append({'fold':fi,'seed':seed,'regime':regime,'status':'completed',
                    'train_n':len(tr),'train_positive':counts[1],'test_n':len(te),
                    'auroc':float(roc_auc_score(yy[te],score)),
                    'average_precision':float(average_precision_score(yy[te],score))})
            print(f'canonical block {fi}, seed {seed} completed',flush=True)
    pred=pd.concat(preds,ignore_index=True)
    pooled=[]
    for (regime,seed),p in pred.groupby(['regime','seed']):
        assert not p.row_index.duplicated().any()
        pooled.append({'regime':regime,'seed':int(seed),'n':len(p),'positive':int(p.label.sum()),
            'auroc':float(roc_auc_score(p.label,p.score)),'average_precision':float(average_precision_score(p.label,p.score))})
    pred.to_csv(OUT/'predictions.csv',index=False)
    pd.DataFrame(fold_rows).to_csv(OUT/'folds.csv',index=False)
    pd.DataFrame(pooled).to_csv(OUT/'summary.csv',index=False)
    (OUT/'manifest.json').write_text(json.dumps({'sourmash_version':version('sourmash'),
        'k':21,'hash_seed':42,'main_sketch_size':2048,'distance_threshold':.05,
        'fasta_checksums':hashes,'protocol_sha256':hashlib.sha256((ROOT/'SECOND_BANK_AND_CANONICAL_PROTOCOL.md').read_bytes()).hexdigest(),
        'test_cohort':'original five dual blocks; same test rows for all regimes',
        'interpretation':'retrospective canonical sketch sensitivity; not ANI or independent validation'},indent=2))
    print(pd.DataFrame(pooled).to_string(index=False),flush=True)


if __name__=='__main__': main()
