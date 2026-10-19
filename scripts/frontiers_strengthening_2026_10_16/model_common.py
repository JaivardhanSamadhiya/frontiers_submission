"""Shared, read-only input assembly for retrospective model strengthening."""
from pathlib import Path
import importlib.util
import hashlib
import json
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
OUT = ROOT / 'data/results_v2/frontiers_strengthening_2026_10_16/models'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

def load_ncbi():
    return module('historical_matched', 'experiments/35_matched_cohort_leakage_audit.py').load_design()

def load_staph():
    """Rebuild deterministic genome descriptors; never modify old input caches."""
    from precisionphage.data import GenomeIndex
    from precisionphage.features.genomic import build_node_features, edge_features_from_spectra, kmer_spectrum
    from precisionphage.splits.canonical import sketch_fasta, distance_matrix
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    df = pd.read_csv(ROOT / 'data/interim_v2/interactions.csv')
    df = df.loc[df.study.eq('StaphStudy')].reset_index(drop=True)
    assert len(df) == 1053 and df.label.sum() == 333
    nodes, mapped, genome_hashes, groups = {}, {}, {}, []
    for axis, directory in [('phage', 'data/external_validation/staphstudy/phages_by_study_id'),
                            ('host', 'external/upstream_vhip_tool/example/host_genomes')]:
        names = sorted(df[axis].unique())
        index = GenomeIndex([ROOT / directory])
        paths = [index.resolve(name) for name in names]
        assert all(path is not None for path in paths)
        genome_hashes.update({axis + ':' + name: sha(path) for name, path in zip(names, paths)})
        nodes[axis], _ = build_node_features(names, index, k=4, use_codon=True, use_dinuc=True, n_workers=2)
        mapped[axis] = df[axis].map(dict(zip(names, range(len(names))))).to_numpy()
        sketches = [sketch_fasta(path) for path in paths]
        distances = distance_matrix(sketches)
        _, labels = connected_components(csr_matrix(distances <= .05), directed=False)
        df[axis + '_cluster'] = df[axis].map(dict(zip(names, labels)))
        groups.extend({'axis': axis, 'entity': name, 'cluster': int(label)} for name, label in zip(names, labels))
    P, H = nodes['phage'][mapped['phage']], nodes['host'][mapped['host']]
    kdim = len(kmer_spectrum('ACGT' * 4, 4))
    local = np.asarray([edge_features_from_spectra(p[:kdim], h[:kdim]) for p, h in zip(P,H)])
    cache = pd.read_csv(ROOT / 'data/interim_v2/external_staph_seq_pair_features.csv')
    extra = df[['phage','host']].merge(cache, on=['phage','host'], how='left', validate='one_to_one')
    extra = extra.drop(columns=['phage','host']).to_numpy(np.float32)
    X = np.hstack([P,H,local,df[['k3dist','k6dist','GCdiff','Homology']].to_numpy(),extra]).astype(np.float32)
    assert X.shape == (1053,448) and np.isfinite(X).all()
    assert df[['phage','host']].drop_duplicates().shape[0] == 39*27
    return df, X, genome_hashes, groups

def metric(y, p):
    from sklearn.metrics import roc_auc_score, average_precision_score
    y, p = np.asarray(y), np.asarray(p)
    return {'n': len(y), 'positive': int(y.sum()), 'prevalence': float(y.mean()),
            'auroc': float(roc_auc_score(y,p)) if len(np.unique(y)) == 2 else None,
            'average_precision': float(average_precision_score(y,p))}

def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2), encoding='utf-8')

def assert_split(df, tr, te, axes):
    assert not np.intersect1d(tr, te).size
    for axis in axes:
        assert set(df.iloc[tr][axis + '_cluster']).isdisjoint(df.iloc[te][axis + '_cluster'])
