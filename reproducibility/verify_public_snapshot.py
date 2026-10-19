"""Check the explicit public scientific payload; no model fits or network calls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
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


def index_entries(raw):
    """Read exact stage-zero index object IDs, preserving only the old gitlink."""
    entries={}
    case_names=set()
    preserved_link=('external/PHIST','849efbe1894e321c1ff06b13f7f1aca377bfddb5')
    links=[]
    if raw and not raw.endswith(b'\0'):
        raise ValueError('Truncated NUL-delimited Git index')
    for record in raw.split(b'\0'):
        if not record:
            continue
        header, separator, encoded_name=record.partition(b'\t')
        fields=header.decode('ascii').split()
        if not separator or len(fields)!=3:
            raise ValueError('Malformed Git index entry')
        mode,identifier,stage=fields
        name=encoded_name.decode('utf-8')
        safe_path(name)
        if (stage!='0' or not re.fullmatch('[0-9a-f]{40}',identifier)
                or name.casefold() in case_names):
            raise ValueError('Unmerged, malformed or case-colliding Git index')
        case_names.add(name.casefold())
        if mode=='160000' and (name,identifier)==preserved_link:
            links.append((name,identifier))
            continue
        if mode not in {'100644','100755'}:
            raise ValueError('Symlink or unreviewed Git index mode/gitlink')
        entries[name]=identifier
    if links!=[preserved_link]:
        raise ValueError('The exact preserved PHIST gitlink must remain in the index')
    return entries


def batch_blobs(raw,identifiers):
    """Parse Git binary batch output without interpreting paths or text EOLs."""
    blobs={}
    offset=0
    for identifier in identifiers:
        end=raw.find(b'\n',offset)
        if end<0:
            raise ValueError('Truncated Git batch header')
        fields=raw[offset:end].decode('ascii').split()
        if len(fields)!=3 or fields[0]!=identifier or fields[1]!='blob' or not fields[2].isdigit():
            raise ValueError('Unexpected Git batch object identity/type')
        size=int(fields[2])
        start,stop=end+1,end+1+size
        if stop>=len(raw) or raw[stop:stop+1]!=b'\n':
            raise ValueError('Truncated Git batch content')
        content=raw[start:stop]
        actual=hashlib.sha1(b'blob '+str(size).encode('ascii')+b'\0'+content).hexdigest()
        if actual!=identifier:
            raise ValueError('Git batch object content does not match its exact blob ID')
        blobs[identifier]=content
        offset=stop+1
    if offset!=len(raw):
        raise ValueError('Unrequested trailing Git batch content')
    return blobs


def build():
    command=['git','-c',f'safe.directory={ROOT.as_posix()}','-C',str(ROOT)]
    raw=subprocess.check_output(command+['ls-files','--stage','-z'])
    entries=index_entries(raw)
    identifiers=sorted(set(entries.values()))
    payload=('\n'.join(identifiers)+'\n').encode('ascii')
    # One binary Git subprocess supplies every exact index blob. Object IDs,
    # types, lengths, binary boundaries and content hashes are all checked.
    batch=subprocess.check_output(command+['cat-file','--batch'],input=payload)
    blobs=batch_blobs(batch,identifiers)
    rows=[]
    for name,identifier in sorted(entries.items()):
        path=safe_path(name)
        if path==MANIFEST:
            continue
        if not path.is_file():
            raise ValueError('Missing tracked regular scientific file: '+name)
        content=path.read_bytes()
        # Git archives contain staged blob bytes, which can differ from a
        # Windows working tree under autocrlf even for extensionless files.
        # Record that exact representation rather than broadening text rules.
        git_blob=blobs[identifier]
        row={'path':name,'bytes':len(content),'sha256_exact':digest(content),
             'git_blob_bytes':len(git_blob),'sha256_git_blob':digest(git_blob)}
        if path.suffix.lower() in TEXT or path.name in TEXT_NAMES:
            row['sha256_lf']=digest(lf_bytes(content))
        rows.append(row)
    manifest={'schema':2,'snapshot_version':'2026.10.18.3',
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
