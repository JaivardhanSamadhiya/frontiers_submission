"""Resolve BASEL and ECOR genomes from primary accession tables, without labels.

Run with the bundled document Python (openpyxl is used only to read the public
supplement). No phenotype table is read. Raw assemblies retain their contigs.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT/'data/external_validation/basel_ecor_2026'
OUT = ROOT/'data/results_v2/basel_ecor_validation'
last_ncbi = 0.


def retrieve(url):
    global last_ncbi
    for attempt in range(4):
        if 'eutils.ncbi.nlm.nih.gov' in url:
            time.sleep(max(0, .4-(time.monotonic()-last_ncbi)))
            last_ncbi = time.monotonic()
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'PrecisionPhage-reproducibility/1.0'})
            with urllib.request.urlopen(request, timeout=90) as r:
                content = r.read(30_000_001)
            if len(content) > 30_000_000: raise ValueError('Unexpectedly large genome response')
            return content
        except Exception:
            if attempt == 3: raise
            time.sleep(2**attempt)


def cached(url, path):
    if not path.exists():
        content = retrieve(url)
        path.write_bytes(content)
    return path.read_bytes()


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    targets = json.loads((OUT/'identifier_manifest.json').read_text())
    supplement = RAW/'pbio.3003063.s002.xlsx'
    cached('https://ndownloader.figshare.com/files/53473839', supplement)
    wb = openpyxl.load_workbook(supplement, read_only=True, data_only=True)
    sheet = wb['list of all phages']
    records = list(sheet.values)
    assert records[6][0] == 'Bas##' and records[6][14] == 'NCBI GenBank accession'
    phages = {}
    for row_number, row in enumerate(records[7:], start=8):
        if row[0] not in targets['phages']: continue
        if row[0] in phages: raise ValueError('Duplicate phage accession mapping')
        accession = str(row[14]).strip()
        if not re.fullmatch(r'[A-Z]{1,2}\d+(?:\.\d+)?', accession): raise ValueError(accession)
        phages[row[0]] = {'identifier': row[0], 'accession': accession,
            'expected_length': int(row[15]), 'genus': row[8], 'published_name': row[10],
            'source': 'https://doi.org/10.1371/journal.pbio.3003063.s002',
            'sheet': sheet.title, 'row': row_number}
    if set(phages) != set(targets['phages']): raise ValueError('Incomplete BASEL accession coverage')
    article_url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6256646/fullTextXML'
    xml = cached(article_url, RAW/'ECOR_Patel_2018.xml')
    tree = ET.fromstring(xml)
    hosts = {}
    for tr in tree.findall('.//table-wrap/table/tbody/tr'):
        values = [''.join(c.itertext()).strip() for c in tr.findall('td')]
        if len(values) < 7: continue
        name = 'ECOR'+values[0]
        if name not in targets['hosts']: continue
        if not re.fullmatch(r'QO[A-Z]{2}0{8}', values[3]): raise ValueError(values)
        hosts[name] = {'identifier': name, 'wgs_master': values[3],
                       'expected_contigs': int(values[5]),
                       'expected_length': int(values[6].replace(',', '')),
                       'source': 'https://doi.org/10.1128/MRA.01133-18'}
    if set(hosts) != set(targets['hosts']): raise ValueError('Incomplete ECOR accession coverage')
    # GenoPHI data availability independently identifies this same WGS series
    # and BioProject, avoiding substitution of the older contaminated series.
    reference = 'https://www.nature.com/articles/s41564-026-02482-5'
    mapping = {'recorded_utc': datetime.now(timezone.utc).isoformat(),
               'genophi_host_provenance': reference,
               'genophi_host_provenance_detail': 'PRJNA230969, QOWX00000000-QOXX00000000; Supplementary Table 9',
               'phages': phages, 'hosts': hosts,
               'sources_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [supplement, RAW/'ECOR_Patel_2018.xml']}}
    (OUT/'genome_accession_mapping.json').write_text(json.dumps(mapping, indent=2), encoding='utf-8')
    downloaded = []
    for kind, entries in [('phage', phages), ('host', hosts)]:
        folder = RAW/f'{kind}_genomes'; folder.mkdir(exist_ok=True)
        for name, meta in entries.items():
            path = folder/f'{name}.fasta'
            info = folder/f'{name}.provenance.json'
            if path.exists() and info.exists():
                row = json.loads(info.read_text())
                if hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']: raise ValueError('Genome changed')
                downloaded.append(row); continue
            if kind == 'phage':
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(
                    {'db': 'nuccore', 'id': meta['accession'], 'rettype': 'fasta', 'retmode': 'text'})
                content = retrieve(url)
                extra = {}
            else:
                query = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'+urllib.parse.urlencode(
                    {'db': 'assembly', 'term': meta['wgs_master'], 'retmode': 'json'})
                search = json.loads(cached(query, folder/f'{name}.assembly_search.json'))
                ids = search['esearchresult']['idlist']
                if len(ids) != 1: raise ValueError(f'Ambiguous assembly for {name}: {ids}')
                summary_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urllib.parse.urlencode(
                    {'db': 'assembly', 'id': ids[0], 'retmode': 'json'})
                result = json.loads(cached(summary_url, folder/f'{name}.assembly_summary.json'))['result'][ids[0]]
                expected_wgs = meta['wgs_master'][:4]+'01'
                if result.get('wgs') != expected_wgs: raise ValueError('Assembly WGS series/version mismatch')
                if not any(p['bioprojectaccn'] == 'PRJNA230969' for p in result.get('gb_bioprojects', [])):
                    raise ValueError('Unexpected ECOR BioProject')
                ftp = result['ftppath_genbank'].replace('ftp://', 'https://')
                url = ftp+'/'+ftp.rsplit('/', 1)[1]+'_genomic.fna.gz'
                content = gzip.decompress(retrieve(url))
                extra = {'assembly_accession': result['synonym']['genbank'],
                         'assembly_name': result['assemblyname'], 'wgs_series_version': result['wgs']}
            text = content.decode('utf-8').strip()
            if not text.startswith('>'): raise ValueError(f'Not FASTA: {name}')
            headers = [line[1:] for line in text.splitlines() if line.startswith('>')]
            sequence = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))
            if not sequence or set(sequence.upper())-set('ACGTRYSWKMBDHVN'): raise ValueError('Invalid DNA')
            if len(sequence) != meta['expected_length']: raise ValueError(f'Published length mismatch {name}: {len(sequence)} vs {meta["expected_length"]}')
            if kind == 'phage' and (len(headers) != 1 or headers[0].split()[0].split('.')[0] != meta['accession'].split('.')[0]):
                raise ValueError('Phage accession mismatch')
            if kind == 'host' and len(headers) != meta['expected_contigs']: raise ValueError('Contig count mismatch')
            path.write_bytes(content)
            row = {'identifier': name, 'kind': kind, 'source_url': url, **extra,
                   'sha256': hashlib.sha256(content).hexdigest(), 'bases': len(sequence),
                   'contigs': len(headers), 'sequence_accessions': [h.split()[0] for h in headers],
                   'retrieved_utc': datetime.now(timezone.utc).isoformat()}
            info.write_text(json.dumps(row, indent=2)); downloaded.append(row)
            print(f'{kind} {name}: {len(sequence)} bases, {len(headers)} contigs', flush=True)
    manifest = {'status': 'all published lengths and accession/contig counts verified',
                'genomes': downloaded, 'outcomes_read': False,
                'contig_handling_note': 'Raw contigs retained. Frozen historical feature extraction concatenates contigs; this limitation must be documented, not silently changed.'}
    (OUT/'genome_download_manifest.json').write_text(json.dumps(manifest, indent=2))
    print(f'Complete: {len(phages)} phages and {len(hosts)} hosts', flush=True)


if __name__ == '__main__': main()
