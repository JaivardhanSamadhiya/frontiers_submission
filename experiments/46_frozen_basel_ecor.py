"""Frozen sequence-only NCBI_HR -> BASEL-ECOR evaluation, staged prediction.

The predict command cannot read the outcome file. A separate evaluate command
requires the saved prediction/overlap lock and refuses altered artifacts.
This is a newly performed public-data test, not prospectively collected data.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
import joblib
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import sourmash

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.data import GenomeIndex
from precisionphage.features.genomic import build_node_features, edge_features_from_spectra, kmer_spectrum
from precisionphage.features.seqmatch import compute_pair_features, pair_feature_cols
from precisionphage.splits.canonical import sketch_fasta
from precisionphage.eval import binary_metrics, calibration_curve_ece
from precisionphage.utils import load_config

RAW = ROOT/'data/external_validation/basel_ecor_2026'
OUT = ROOT/'data/results_v2/basel_ecor_validation'
MODEL = ROOT/'data/results_v2/external_validation_models/sequence_only.joblib'
CANONICAL = ROOT/'data/results_v2/canonical_sensitivity'
SEED = 20261005


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024*1024), b''): h.update(chunk)
    return h.hexdigest()


def sequence_digest(path):
    # Invariant to FASTA line wrapping, headers, contig order/orientation.
    # Different assembly fragmentation is not called an exact genome match.
    hashes = []
    for r in SeqIO.parse(path, 'fasta'):
        seq = str(r.seq).upper().replace('U', 'T')
        rc = str(r.seq.reverse_complement()).upper().replace('U', 'T')
        hashes.append(hashlib.sha256(min(seq, rc).encode()).hexdigest())
    return hashlib.sha256('\n'.join(sorted(hashes)).encode()).hexdigest()


def feature_names(cfg):
    return ([f'phage_node_{i}' for i in range(212)] + [f'host_node_{i}' for i in range(212)] +
            ['cos_dist', 'l1', 'pearson', 'jaccard'] + pair_feature_cols(cfg))


def assemble(pairs, cfg, pair_features=None):
    phages, hosts = sorted(pairs.phage.unique()), sorted(pairs.host.unique())
    pidx, hidx = GenomeIndex([cfg['paths']['phage_fasta_dir']]), GenomeIndex([cfg['paths']['host_fasta_dir']])
    for names, idx in [(phages, pidx), (hosts, hidx)]:
        if idx.coverage(names)['resolved'] != len(names): raise ValueError('Incomplete exact genome resolution')
    f = cfg['features']
    P, pm = build_node_features(phages, pidx, k=f['kmer_k'], use_codon=f['use_codon'], use_dinuc=f['use_dinuc'], n_workers=1)
    H, hm = build_node_features(hosts, hidx, k=f['kmer_k'], use_codon=f['use_codon'], use_dinuc=f['use_dinuc'], n_workers=1)
    if not pm.all() or not hm.all() or P.shape[1] != 212 or H.shape[1] != 212: raise ValueError('Invalid node schema')
    pi = pairs.phage.map(dict(zip(phages, range(len(phages))))).to_numpy()
    hi = pairs.host.map(dict(zip(hosts, range(len(hosts))))).to_numpy()
    kdim = len(kmer_spectrum('ACGT'*f['kmer_k'], f['kmer_k']))
    edge = np.stack([edge_features_from_spectra(P[i, :kdim], H[j, :kdim]) for i, j in zip(pi, hi)])
    if pair_features is None:
        pair_features = compute_pair_features(pairs[['phage', 'host']], cfg, n_jobs=1)
    if pair_features.duplicated(['phage', 'host']).any(): raise ValueError('Duplicate feature rows')
    extra = pairs[['phage', 'host']].merge(pair_features, on=['phage', 'host'], how='left', validate='one_to_one')
    X = np.hstack([P[pi], H[hi], edge, extra[pair_feature_cols(cfg)].to_numpy(np.float32)]).astype(np.float32)
    if X.shape != (len(pairs), 444) or not np.isfinite(X).all(): raise ValueError('Invalid frozen feature matrix')
    return X, pair_features


def regression_check(cfg, frozen):
    # Verify feature order against saved predictions on a previously evaluated
    # cohort, never the new outcomes. Covers every phage against its first host.
    old = pd.read_csv(ROOT/'data/results_v2/external_staph_predictions.csv').drop_duplicates('phage')
    probe = old.iloc[:3].reset_index(drop=True)
    check_cfg = load_config(ROOT/'configs/default.yaml')
    check_cfg['paths']['phage_fasta_dir'] = ROOT/'data/external_validation/staphstudy/phages_by_study_id'
    check_cfg['paths']['host_fasta_dir'] = ROOT/'external/upstream_vhip_tool/example/host_genomes'
    cache = pd.read_csv(ROOT/'data/interim_v2/external_staph_seq_pair_features.csv')
    X, _ = assemble(probe, check_cfg, cache)
    p = frozen['model'].predict_proba(frozen['scaler'].transform(X))[:, 1]
    delta = float(np.max(np.abs(p-probe.prob_sequence_only.to_numpy())))
    if delta > 1e-6: raise ValueError(f'Saved-model regression mismatch: {delta}')
    return {'previously_evaluated_pairs': len(probe), 'maximum_absolute_difference': delta,
            'note': 'Regression check used the existing pair-feature cache and recomputed node/composition features.'}


def overlaps(cfg, pairs):
    result, detail = {}, []
    for axis, rawdir in [('phage', 'phage_genomes'), ('host', 'host_genomes')]:
        names = sorted(pairs[axis].unique())
        train_names = np.load(CANONICAL/f'{axis}_distances_2048.npz')['names'].tolist()
        source_hashes = json.loads((CANONICAL/f'{axis}_fasta_checksums.json').read_text())
        if sorted(source_hashes) != train_names: raise ValueError('Canonical training cache order mismatch')
        archive = np.load(CANONICAL/f'{axis}_canonical_2048.npz')
        train_idx = GenomeIndex([cfg['paths'][f'{axis}_fasta_dir']])
        train_sketches, exact = [], {}
        for i, name in enumerate(train_names):
            path = train_idx.resolve(name)
            if path is None or sha(path) != source_hashes[name]: raise ValueError('Canonical reference genome changed')
            train_sketches.append(sourmash.MinHash(n=2048, ksize=21, seed=42, mins=archive[f's{i}'].tolist()))
            exact.setdefault(sequence_digest(path), []).append(name)
        sketches = [sketch_fasta(RAW/rawdir/f'{name}.fasta') for name in names]
        dd = np.eye(len(names))
        for i, a in enumerate(sketches):
            distances = []
            for b in train_sketches:
                j = a.jaccard(b)
                distances.append(float(-np.log(2*j/(1+j))/21) if j else 1.)
            nearest = int(np.argmin(distances))
            match = exact.get(sequence_digest(RAW/rawdir/f'{names[i]}.fasta'), [])
            detail.append({'axis': axis, 'entity': names[i], 'nearest_training_entity': train_names[nearest],
                           'nearest_training_distance': distances[nearest],
                           'exact_contig_set_matches': '|'.join(match), 'identifier_seen': names[i] in train_names})
            for j in range(i):
                similarity = a.jaccard(sketches[j])
                distance = float(-np.log(2*similarity/(1+similarity))/21) if similarity else 1.
                dd[i, j] = dd[j, i] = distance
        np.fill_diagonal(dd, 0.)
        n, group = connected_components(csr_matrix(dd <= .05), directed=False)
        for row, g in zip(detail[-len(names):], group): row['cluster'] = int(g)
        result[axis] = {'entities': len(names), 'training_entities': len(train_names), 'sequence_clusters': int(n),
                        'at_or_below_0_05': sum(r['nearest_training_distance'] <= .05 for r in detail[-len(names):]),
                        'exact_contig_set_matches': sum(bool(r['exact_contig_set_matches']) for r in detail[-len(names):])}
        print(f'Overlap audit finished for {axis}: {result[axis]}', flush=True)
    table = pd.DataFrame(detail)
    table.to_csv(OUT/'sequence_overlap.csv', index=False)
    return result


def predict():
    if (OUT/'evaluation.json').exists(): raise ValueError('Outcomes already evaluated; do not replace frozen predictions')
    if (OUT/'prediction_lock.json').exists(): raise ValueError('Predictions already locked; verify/evaluate instead of replacing')
    cfg = load_config(ROOT/'configs/default.yaml')
    frozen = joblib.load(MODEL)
    if frozen['scaler'].n_features_in_ != 444 or frozen['model'].n_features_in_ != 444: raise ValueError('Wrong frozen model')
    mapping = json.loads((OUT/'genome_download_manifest.json').read_text())
    for r in mapping['genomes']:
        if sha(RAW/f'{r["kind"]}_genomes'/f'{r["identifier"]}.fasta') != r['sha256']: raise ValueError('Genome altered')
    pairs = pd.read_csv(OUT/'eligible_input_identifiers.csv').rename(columns={'strain': 'host'})
    if list(pairs.columns) != ['host', 'phage'] or len(pairs) != 1240: raise ValueError('Unexpected input schema')
    if pairs.duplicated(['phage', 'host']).any(): raise ValueError('Duplicate candidates')
    old_check = regression_check(cfg, frozen)
    print(f'Saved-model regression check: {old_check}', flush=True)
    ov = overlaps(cfg, pairs)
    cfg['paths']['phage_fasta_dir'] = RAW/'phage_genomes'
    cfg['paths']['host_fasta_dir'] = RAW/'host_genomes'
    X, pair_features = assemble(pairs, cfg)
    pair_features.to_csv(OUT/'pair_features.csv', index=False)
    np.savez_compressed(OUT/'frozen_input_matrix.npz', X=X, feature_names=feature_names(cfg))
    pairs['score'] = frozen['model'].predict_proba(frozen['scaler'].transform(X))[:, 1]
    if not np.isfinite(pairs.score).all(): raise ValueError('Nonfinite predictions')
    pairs.to_csv(OUT/'frozen_predictions.csv', index=False)
    files = [MODEL, OUT/'frozen_predictions.csv', OUT/'frozen_input_matrix.npz', OUT/'pair_features.csv',
             OUT/'sequence_overlap.csv', OUT/'genome_download_manifest.json', OUT/'identifier_manifest.json',
             OUT/'eligible_input_identifiers.csv', ROOT/'configs/default.yaml', Path(__file__),
             ROOT/'NEXT_STAGE_VALIDATION_PROTOCOL_2026_10_05.md']
    files += [ROOT/'src/precisionphage'/p for p in ['features/genomic.py', 'features/seqmatch.py',
               'features/proteins.py', 'data/genomes.py', 'splits/canonical.py']]
    lock = {'locked_utc': datetime.now(timezone.utc).isoformat(), 'model': 'previously saved sequence_only.joblib',
            'fit_operations': 0, 'outcomes_read': False, 'pairs': len(pairs), 'overlap': ov,
            'regression_check': old_check,
            'bootstrap': '2000 two-way canonical sequence-cluster replicates, k21/2048 hashes, threshold 0.05',
            'limitations': ['The full 448-input model cannot be evaluated without four unavailable VHIP columns.',
                           'The 444-input sequence-only model was selected by input compatibility before this evaluation.',
                           'Historic frozen features concatenate contigs; no feature definition was changed for this test.',
                           'Public data; staged prediction lock is not independent certification of historical blindness.',
                           'Binary individual phage-host assay labels, not cocktail efficacy or clinical outcomes.'],
            'sha256': {str(p.relative_to(ROOT)): sha(p) for p in files}}
    (OUT/'prediction_lock.json').write_text(json.dumps(lock, indent=2), encoding='utf-8')
    print('PREDICTIONS LOCKED; outcomes not evaluated', flush=True)


def evaluate():
    lock = json.loads((OUT/'prediction_lock.json').read_text())
    for path, checksum in lock['sha256'].items():
        if sha(ROOT/path) != checksum: raise ValueError(f'Locked artifact changed: {path}')
    manifest = json.loads((OUT/'identifier_manifest.json').read_text())
    sealed = ROOT/manifest['sealed_outcomes']
    if sha(sealed) != manifest['sha256']: raise ValueError('Source outcome checksum changed')
    pred = pd.read_csv(OUT/'frozen_predictions.csv')
    outcome = pd.read_csv(sealed).rename(columns={'strain': 'host', 'interaction': 'label'})
    if outcome.duplicated(['host', 'phage']).any() or not np.isin(outcome.label, [0, 1]).all(): raise ValueError('Invalid outcomes')
    table = pred.merge(outcome, on=['host', 'phage'], how='outer', validate='one_to_one', indicator=True)
    if (table._merge != 'both').any(): raise ValueError('Prediction/outcome cohort differs')
    table.drop(columns='_merge', inplace=True)
    ov = pd.read_csv(OUT/'sequence_overlap.csv')
    for axis in ['host', 'phage']:
        a = ov[ov.axis == axis].set_index('entity')
        table[f'{axis}_cluster'] = table[axis].map(a.cluster)
        table[f'{axis}_nearest_distance'] = table[axis].map(a.nearest_training_distance)
    spec = importlib.util.spec_from_file_location('frozen_staph_helpers', ROOT/'experiments/26_frozen_cross_source_validation.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    y, score = table.label.to_numpy(int), table.score.to_numpy()
    metrics = binary_metrics(y, score)
    clusters = [table.phage_cluster.nunique(), table.host_cluster.nunique()]
    ci = (module._two_way_cluster_ci(y, score, table.phage_cluster.to_numpy(), table.host_cluster.to_numpy(), SEED)
          if min(clusters) >= 2 else {'unavailable': 'Fewer than two sequence clusters on an entity axis'})
    calibration = calibration_curve_ece(y, score, n_bins=10)
    subsets = {}
    for name, mask in [('host_beyond_0_05', table.host_nearest_distance > .05),
                       ('phage_beyond_0_05', table.phage_nearest_distance > .05),
                       ('both_axes_beyond_0_05', (table.host_nearest_distance > .05) & (table.phage_nearest_distance > .05))]:
        z = table[mask]
        subsets[name] = {'pairs': len(z), 'positive': int(z.label.sum()),
                         'metrics': binary_metrics(z.label.to_numpy(), z.score.to_numpy()) if len(z) else None}
    result = {'evaluated_utc': datetime.now(timezone.utc).isoformat(), 'prediction_locked_utc': lock['locked_utc'],
              'prediction_lock_sha256': sha(OUT/'prediction_lock.json'), 'metrics': metrics,
              'ci95': ci, 'calibration': calibration, 'overlap': lock['overlap'],
              'distance_subsets': subsets, 'limitations': lock['limitations'],
              'uncertainty_note': 'Cluster bootstrap conditional on saved model, with only 25 candidate host entities; intervals may be unstable.'}
    table.to_csv(OUT/'evaluated_predictions.csv', index=False)
    (OUT/'evaluation.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('stage', choices=['predict', 'evaluate'])
    args = parser.parse_args()
    predict() if args.stage == 'predict' else evaluate()
