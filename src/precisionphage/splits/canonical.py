"""Canonical DNA sketches for sensitivity analyses; historical caches unchanged."""
from pathlib import Path
import gzip
import numpy as np
from Bio import SeqIO
import sourmash


def sketch_fasta(path: Path, num: int = 2048, k: int = 21):
    """Sketch records separately so no artificial k-mers cross contig boundaries."""
    mh = sourmash.MinHash(n=num, ksize=k, seed=42)
    opener = gzip.open if Path(path).suffix == '.gz' else open
    with opener(path, 'rt') as stream:
        for record in SeqIO.parse(stream, 'fasta'):
            mh.add_sequence(str(record.seq).upper().replace('U', 'T'), force=True)
    if not len(mh):
        raise ValueError(f'No valid DNA k-mers in {path}')
    return mh


def distance_matrix(sketches, k=21):
    result = np.zeros((len(sketches), len(sketches)), dtype=np.float64)
    for i in range(len(sketches)):
        for j in range(i):
            similarity = sketches[i].jaccard(sketches[j])
            distance = -np.log(2 * similarity / (1 + similarity)) / k if similarity > 0 else 1.0
            result[i, j] = result[j, i] = distance
        if i and i % 250 == 0:
            print(f'compared {i}/{len(sketches)} sketches', flush=True)
    return result
