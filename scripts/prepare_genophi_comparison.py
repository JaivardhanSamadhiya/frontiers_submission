"""Resolve exact public proteomes for a retrospective Pseudomonas comparison.

No fitting or outcome-dependent genome selection. Cached downloads are never
silently replaced; a changed checksum requires explicit investigation.
"""
from pathlib import Path
import hashlib
import io
import json
import sys
import time
import zipfile
from datetime import datetime, timezone

import pandas as pd
import requests
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/external_validation/genophi_comparison'
OUT = ROOT / 'data/results_v2/genophi_comparison'
MATRIX = ROOT / 'data/external_validation/genophi_2026/processed_data.zip'
EXPECTED = '1d98c77fc4068315246a4eb56aa185fc769c06bf70fd20e5090856b25f94a00d'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def download(url, path):
    """Rate-limited, capped download with immutable cache provenance."""
    sidecar = path.with_suffix(path.suffix + '.source.json')
    if path.exists():
        data = path.read_bytes()
        meta = json.loads(sidecar.read_text())
        if meta['url'] != url or meta['sha256'] != sha(data):
            raise ValueError(f'Cache mismatch: {path}')
        return data
    for attempt in range(4):
        time.sleep(0.4 if attempt == 0 else 2 ** attempt)
        try:
            response = requests.get(url, timeout=90)
            response.raise_for_status()
            data = response.content
            if not data or len(data) > 30 * 1024 * 1024:
                raise ValueError(f'Unexpected response size: {len(data)}')
            break
        except requests.RequestException:
            if attempt == 3:
                raise
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    sidecar.write_text(json.dumps({'url': url, 'sha256': sha(data),
        'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'bytes': len(data)}, indent=2))
    return data


def protein_records(data):
    text = data.decode('utf-8')
    if not text.lstrip().startswith('>'):
        raise ValueError('Response is not protein FASTA')
    records = list(SeqIO.parse(io.StringIO(text), 'fasta'))
    ids = [r.id for r in records]
    if not records or len(ids) != len(set(ids)):
        raise ValueError('Empty proteome or duplicate protein identifiers')
    alphabet = set('ACDEFGHIKLMNPQRSTVWYBXZJUO*')
    for record in records:
        if not record.seq or not set(str(record.seq).upper()) <= alphabet:
            raise ValueError(f'Invalid amino acid sequence: {record.id}')
    return records


def immutable_write(path, data):
    if path.exists() and path.read_bytes() != data:
        raise ValueError(f'Refusing to replace changed input: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        path.write_bytes(data)


def host_proteins(entity):
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{entity}/dataset_report'
    report = json.loads(download(url, RAW / 'reports' / f'{entity}.json'))
    rows = report.get('reports', [])
    if len(rows) != 1 or rows[0]['accession'].split('.')[0] != entity:
        raise ValueError(f'Non-unique or mismatched host resolution: {entity}')
    row = rows[0]
    source = row['accession']
    paired = row.get('paired_accession')
    attempts = [source] + ([paired] if paired else [])
    for accession in attempts:
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download?include_annotation_type=PROT_FASTA'
        archive = download(url, RAW / 'archives' / f'{accession}.zip')
        with zipfile.ZipFile(io.BytesIO(archive)) as z:
            member = f'ncbi_dataset/data/{accession}/protein.faa'
            if member not in z.namelist():
                continue
            proteins = z.read(member)
            protein_records(proteins)
            return proteins, {'requested': entity, 'resolved_genbank': source,
                'protein_assembly': accession, 'used_paired_refseq': accession != source,
                'organism': row.get('organism'), 'archive_member': member, 'url': url}
    raise ValueError(f'No supplied protein annotation: {entity}')


def main():
    if sha(MATRIX.read_bytes()) != EXPECTED:
        raise ValueError('Published matrix archive changed')
    with zipfile.ZipFile(MATRIX) as z:
        data = pd.read_csv(z.open('GenoPHI_processed_data/interaction_matrices_all_datasets.csv'))
    bank = data[data.dataset == 'Pseudomonas'].copy()
    if len(bank) != 437 or bank.duplicated(['strain', 'phage']).any():
        raise ValueError('Unexpected Pseudomonas matrix')
    OUT.mkdir(parents=True, exist_ok=True)
    successes, failures = [], []
    for kind, entities in [('host', sorted(bank.strain.unique())),
                           ('phage', sorted(bank.phage.unique()))]:
        for entity in entities:
            try:
                if kind == 'host':
                    proteins, meta = host_proteins(entity)
                else:
                    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
                           f'?db=nuccore&id={entity}&rettype=fasta_cds_aa&retmode=text')
                    proteins = download(url, RAW / 'phage_raw' / f'{entity}.faa')
                    records = protein_records(proteins)
                    prefixes = {r.id.split('|')[1].split('_prot_')[0] for r in records}
                    if len(prefixes) != 1 or next(iter(prefixes)).split('.')[0] != entity:
                        raise ValueError('Phage accession mismatch')
                    meta = {'requested': entity, 'resolved_nucleotide': next(iter(prefixes)), 'url': url}
                records = protein_records(proteins)
                path = RAW / f'{kind}_proteins' / f'{entity}.faa'
                immutable_write(path, proteins)
                successes.append({'kind': kind, 'entity': entity, **meta,
                    'proteins': len(records), 'amino_acids': sum(len(r.seq) for r in records),
                    'sha256': sha(proteins), 'path': path.relative_to(ROOT).as_posix()})
                print(f'{kind} {entity}: {len(records)} proteins', flush=True)
            except Exception as exc:
                failures.append({'kind': kind, 'entity': entity, 'error': repr(exc)})
                print(f'FAILED {kind} {entity}: {exc}', flush=True)
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(),
        'matrix_archive_sha256': EXPECTED,
        'protocol_sha256': sha((ROOT/'GENOPHI_COMPARISON_PROTOCOL.md').read_bytes()),
        'expected_hosts': 23, 'expected_phages': 19,
        'complete': len(successes) == 42 and not failures,
        'proteomes': successes, 'failures': failures,
        'stage': 'input preparation only; no predictor comparison has been fitted'}
    (OUT/'input_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps({'complete': manifest['complete'], 'resolved': len(successes), 'failures': failures}))
    if failures:
        sys.exit(1)


if __name__ == '__main__':
    main()
