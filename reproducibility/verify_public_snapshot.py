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
TEXT_NAMES = {'.gitignore', '.gitattributes', '.gitmodules'}


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
        # Git archives contain staged blob bytes, which can differ from a
        # Windows working tree under autocrlf even for extensionless files.
        # Record that exact representation rather than broadening text rules.
        git_blob=subprocess.check_output(['git','-c',f'safe.directory={ROOT.as_posix()}',
            '-C',str(ROOT),'cat-file','blob',':'+name])
        row={'path':name,'bytes':len(content),'sha256_exact':digest(content),
             'git_blob_bytes':len(git_blob),'sha256_git_blob':digest(git_blob)}
        if path.suffix.lower() in TEXT or path.name in TEXT_NAMES:
            row['sha256_lf']=digest(lf_bytes(content))
        rows.append(row)
    manifest={'schema':2,'snapshot_version':'2026.10.15.1',
        'scope':'Explicit tracked scientific payload; excludes this self-referential inventory and submodule contents.',
        'text_policy':'Accept the exact recorded working-tree or staged Git-blob SHA-256. Explicitly inventoried text may also match its recorded CRLF/CR-to-LF-normalized SHA-256. No other normalization or unrecorded representation is accepted.',
        'files':rows}
    MANIFEST.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(f'Built public inventory for {len(rows)} tracked files')


def verify():
    manifest=json.loads(MANIFEST.read_text(encoding='utf-8'))
    exact=git_blob=normalized=0
    for row in manifest['files']:
        path=safe_path(row['path'])
        if not path.is_file():
            raise AssertionError('Missing public payload file: '+row['path'])
        content=path.read_bytes()
        observed=digest(content)
        if observed==row['sha256_exact']:
            exact+=1
        elif row.get('sha256_git_blob') and observed==row['sha256_git_blob']:
            git_blob+=1
        elif row.get('sha256_lf') and digest(lf_bytes(content))==row['sha256_lf']:
            normalized+=1
        else:
            raise AssertionError('Public payload checksum mismatch: '+row['path'])
    print(json.dumps({'status':'passed','files':len(manifest['files']),
                      'exact_byte_matches':exact,'recorded_git_blob_matches':git_blob,
                      'documented_text_eol_matches':normalized},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true',help='Maintainer-only: inventory currently tracked files before a release commit')
    args=parser.parse_args()
    if args.build: build()
    else: verify()
