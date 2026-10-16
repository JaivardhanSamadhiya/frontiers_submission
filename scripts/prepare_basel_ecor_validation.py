"""Quarantine the published BASEL–ECOR labels; emit identifier metadata only."""
from pathlib import Path
import csv
import hashlib
import io
import json
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/external_validation/basel_ecor_2026'
OUT = ROOT / 'data/results_v2/basel_ecor_validation'
COMMIT = 'f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f'
URL = f'https://raw.githubusercontent.com/Noonanav/GenoPHI/{COMMIT}/data/experimental_validation/BASEL_ECOR_interaction_matrix.csv'


def main():
    RAW.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
    path = RAW / 'SEALED_BASEL_ECOR_interaction_matrix.csv'
    if path.exists():
        content = path.read_bytes()
    else:
        with urllib.request.urlopen(URL, timeout=60) as response:
            content = response.read(2_000_001)
        if len(content) > 2_000_000: raise ValueError('Unexpected matrix size')
        path.write_bytes(content)
    reader = csv.DictReader(io.StringIO(content.decode('utf-8-sig')))
    headers = reader.fieldnames
    rows = list(reader)
    columns = {v.lower():v for v in headers}
    if 'strain' not in columns or 'phage' not in columns:
        raise ValueError(f'Unrecognized schema; no outcomes displayed: {headers}')
    pairs = [(row[columns['strain']], row[columns['phage']]) for row in rows]
    if len(set(pairs)) != len(pairs): raise ValueError('Duplicate identifier pairs')
    hosts = sorted({h for h,p in pairs}); phages = sorted({p for h,p in pairs})
    result = {'stage':'identifier-only preparation; no outcomes evaluated',
              'source_url':URL, 'commit':COMMIT, 'retrieved_or_checked_utc':datetime.now(timezone.utc).isoformat(),
              'sha256':hashlib.sha256(content).hexdigest(), 'bytes':len(content),
              'columns':headers, 'pairs':len(pairs), 'hosts':hosts, 'phages':phages,
              'sealed_outcomes':str(path.relative_to(ROOT)),
              'limitations':'A preserved staged workflow does not independently prove historical blindness.'}
    old = OUT/'identifier_manifest.json'
    if old.exists() and json.loads(old.read_text())['sha256'] != result['sha256']:
        raise ValueError('External matrix changed')
    old.write_text(json.dumps(result, indent=2), encoding='utf-8')
    with (OUT/'eligible_input_identifiers.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['host','phage']);writer.writerows(pairs)
    print(json.dumps(result,indent=2))


if __name__ == '__main__': main()
