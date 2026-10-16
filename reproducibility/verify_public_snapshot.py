"""Check the explicit public scientific payload; no model fits or network calls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT/'reproducibility/public_payload_sha256.json'
TEXT = {'.py','.md','.txt','.json','.csv','.tsv','.yaml','.yml','.toml','.cff','.log','.xml'}


def digest(data):
    return hashlib.sha256(data).hexdigest()


def lf_bytes(data):
    return data.replace(b'\r\n',b'\n').replace(b'\r',b'\n')


def safe_path(name):
    p=PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or '\\' in name or '.git' in p.parts:
        raise ValueError('Unsafe public inventory path')
    result=(ROOT/name).resolve()
    if not result.is_relative_to(ROOT.resolve()):
        raise ValueError('Inventory path escapes checkout')
    return result


def build():
    raw=subprocess.check_output(['git','-c',f'safe.directory={ROOT.as_posix()}',
        '-C',str(ROOT),'ls-files','-z'])
    names=[x.decode('utf-8') for x in raw.split(b'\0') if x]
    rows=[]
    for name in sorted(set(names)):
        path=safe_path(name)
        if path==MANIFEST or not path.is_file():
            continue
        content=path.read_bytes()
        row={'path':name,'bytes':len(content),'sha256_exact':digest(content)}
        if path.suffix.lower() in TEXT or path.name=='.gitignore':
            row['sha256_lf']=digest(lf_bytes(content))
        rows.append(row)
    manifest={'schema':1,'snapshot_version':'2026.10.15',
        'scope':'Explicit tracked scientific payload; excludes this self-referential inventory and submodule contents.',
        'text_policy':'Exact bytes preferred; explicitly inventoried text may also match after CRLF/CR-to-LF normalization. No other normalization is accepted.',
        'files':rows}
    MANIFEST.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'Built public inventory for {len(rows)} tracked files')


def verify():
    manifest=json.loads(MANIFEST.read_text(encoding='utf-8'))
    exact=normalized=0
    for row in manifest['files']:
        path=safe_path(row['path'])
        if not path.is_file():
            raise AssertionError('Missing public payload file: '+row['path'])
        content=path.read_bytes()
        if digest(content)==row['sha256_exact']:
            exact+=1
        elif row.get('sha256_lf') and digest(lf_bytes(content))==row['sha256_lf']:
            normalized+=1
        else:
            raise AssertionError('Public payload checksum mismatch: '+row['path'])
    print(json.dumps({'status':'passed','files':len(manifest['files']),
                      'exact_byte_matches':exact,'documented_text_eol_matches':normalized},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true',help='Maintainer-only: inventory currently tracked files before a release commit')
    args=parser.parse_args()
    if args.build: build()
    else: verify()
