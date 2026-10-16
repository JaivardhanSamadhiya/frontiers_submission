"""Fetch pinned public data archives used by experiment 41, verifying SHA-256."""
from pathlib import Path
import hashlib
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/external_validation/genophi_2026'
FILES={
    'processed_data.zip':(
        'https://zenodo.org/api/records/21696332/files/GenoPHI_processed_data.zip/content',
        '1d98c77fc4068315246a4eb56aa185fc769c06bf70fd20e5090856b25f94a00d'),
    'published_phylogenies.zip':(
        'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-026-02482-5/MediaObjects/41564_2026_2482_MOESM4_ESM.zip',
        '166864bcf2f2931f8becb7286ef1db1f3efc7a8ebe23480b1a624b9f795f9cd2')}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for name,(url,expected) in FILES.items():
        path=OUT/name
        if path.exists():
            content=path.read_bytes()
        else:
            with urllib.request.urlopen(url,timeout=60) as response:
                content=response.read(20*1024*1024+1)
            if len(content)>20*1024*1024: raise RuntimeError('Unexpected archive size')
        if hashlib.sha256(content).hexdigest()!=expected:
            raise RuntimeError(f'Input checksum mismatch: {name}; preserve and investigate')
        if not path.exists(): path.write_bytes(content)
        print(f'verified {name}')


if __name__=='__main__': main()
