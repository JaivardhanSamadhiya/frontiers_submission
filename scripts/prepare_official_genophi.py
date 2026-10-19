"""Fetch pinned official software only; exclude phenotype tables and notebooks."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.request
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
COMMIT='f7e30abc1b0bbba1ae76b8f58a3363bbbc5cd54f'
DEST=ROOT/'data/external_validation/genophi_comparison/official_source'
OUT=ROOT/'data/results_v2/genophi_comparison'


def get(url):
    req=urllib.request.Request(url,headers={'User-Agent':'PrecisionPhage-reproducibility-audit'})
    with urllib.request.urlopen(req,timeout=60) as response:
        data=response.read(2_000_001)
    if len(data)>2_000_000:raise ValueError('Unexpected individual source size')
    return data


def main():
    tree=json.loads(get(f'https://api.github.com/repos/Noonanav/GenoPHI/git/trees/{COMMIT}?recursive=1'))
    members=[r['path'] for r in tree['tree'] if r['type']=='blob' and
             (r['path'].startswith('genophi/') and r['path'].endswith('.py') or
              r['path'] in ['README.md','LICENSE','setup.py','pyproject.toml'])]
    def fetch(member):
        url=f'https://raw.githubusercontent.com/Noonanav/GenoPHI/{COMMIT}/{member}'
        path=DEST/member
        if path.exists(): content=path.read_bytes()
        else:
            content=get(url);path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(content)
        return {'path':member,'url':url,'sha256':hashlib.sha256(content).hexdigest(),'bytes':len(content)}
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:records=list(pool.map(fetch,members))
    report={'commit':COMMIT,'recorded_utc':datetime.now(timezone.utc).isoformat(),
            'scope':'official software only; no external outcome files or result notebooks', 'files':records}
    path=OUT/'official_source_manifest.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old['files']!=records:raise ValueError('Previously pinned software differs')
    path.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({'files':len(records),'bytes':sum(r['bytes'] for r in records),'commit':COMMIT,'source':str(DEST)}))


if __name__=='__main__':main()
