"""Pre-outcome cluster definitions and secondary BASEL uncertainty diagnostics."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from precisionphage.splits.canonical import sketch_fasta, distance_matrix
OUT = ROOT/'data/results_v2/basel_ecor_validation'
RAW = ROOT/'data/external_validation/basel_ecor_2026'


def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare():
    if (OUT/'evaluation.json').exists(): raise ValueError('Cannot claim pre-outcome definitions after evaluation')
    if (OUT/'uncertainty_definition_lock.json').exists(): raise ValueError('Definitions already frozen')
    targets = json.loads((OUT/'identifier_manifest.json').read_text())
    names = sorted(targets['hosts'])
    dd = distance_matrix([sketch_fasta(RAW/'host_genomes'/f'{n}.fasta') for n in names])
    frame = pd.DataFrame({'host': names, 'identity_cluster': np.arange(len(names))})
    counts = {}
    for threshold in [.001, .01, .05]:
        count, group = connected_components(csr_matrix(dd <= threshold), directed=False)
        frame[f'cluster_{threshold}'] = group; counts[str(threshold)] = int(count)
    path = OUT/'secondary_host_clusters.csv'; frame.to_csv(path, index=False)
    protocol = ROOT/'BASEL_ECOR_UNCERTAINTY_AMENDMENT_2026_10_05.md'
    lock = {'recorded_utc': datetime.now(timezone.utc).isoformat(), 'host_cluster_counts': counts,
            'protocol_sha256': digest(protocol), 'assignments_sha256': digest(path),
            'script_sha256': digest(__file__), 'outcomes_read': False}
    (OUT/'uncertainty_definition_lock.json').write_text(json.dumps(lock, indent=2))
    print(json.dumps(lock, indent=2))


def evaluate():
    lock = json.loads((OUT/'uncertainty_definition_lock.json').read_text())
    if digest(__file__) != lock['script_sha256']: raise ValueError('Script changed')
    if digest(ROOT/'BASEL_ECOR_UNCERTAINTY_AMENDMENT_2026_10_05.md') != lock['protocol_sha256']: raise ValueError('Protocol changed')
    if digest(OUT/'secondary_host_clusters.csv') != lock['assignments_sha256']: raise ValueError('Assignments changed')
    predictions = pd.read_csv(OUT/'evaluated_predictions.csv')
    frame = predictions.merge(pd.read_csv(OUT/'secondary_host_clusters.csv'), on='host', validate='many_to_one')
    spec = importlib.util.spec_from_file_location('frozen_helpers', ROOT/'experiments/26_frozen_cross_source_validation.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    results = []
    for name, hc, pc in [('host_0.001_phage_0.05', frame['cluster_0.001'], frame.phage_cluster),
                         ('host_0.01_phage_0.05', frame['cluster_0.01'], frame.phage_cluster),
                         ('host_identity_phage_identity', frame.host, frame.phage)]:
        h, p = hc.nunique(), pc.nunique()
        ci = module._two_way_cluster_ci(frame.label.to_numpy(int), frame.score.to_numpy(),
                                       pc.to_numpy(), hc.to_numpy(), 20261005) if min(h, p) >= 2 else {'unavailable': 'Fewer than two clusters on an axis'}
        results.append({'definition': name, 'host_groups': h, 'phage_groups': p, 'ci95': ci})
    result = {'status': 'secondary conditional uncertainty; not a replacement for primary unavailable cluster interval',
              'definition_lock_sha256': digest(OUT/'uncertainty_definition_lock.json'),
              'evaluation_sha256': digest(OUT/'evaluation.json'), 'results': results,
              'caveat': 'These resamplings do not establish transfer beyond the sampled Escherichia source or account for all training uncertainty.'}
    (OUT/'secondary_uncertainty.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('stage', choices=['prepare', 'evaluate'])
    a = p.parse_args(); prepare() if a.stage == 'prepare' else evaluate()
