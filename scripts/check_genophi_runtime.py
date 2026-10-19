"""Isolated, CPU-only MMseqs2 smoke check; does not evaluate a predictor."""
from pathlib import Path
import json
import os
import subprocess
import sys
import zipfile

from Bio import SeqIO
from prepare_genophi_comparison import RAW, OUT, ROOT, download, immutable_write, sha


def main():
    source = 'https://github.com/soedinglab/MMseqs2/releases/download/15-6f452/mmseqs-win64.zip'
    archive = RAW/'software/mmseqs-win64-15-6f452.zip'
    download(source, archive)
    destination = RAW/'software/runtime'
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for member in z.infolist():
            if member.is_dir():
                continue
            target = (destination/member.filename).resolve()
            if not target.is_relative_to(destination.resolve()):
                raise ValueError('Unsafe archive member')
            immutable_write(target, z.read(member))
    bindir = destination/'mmseqs/bin'
    # Direct invocation avoids the upstream batch wrapper's optional UAC path.
    # BusyBox's own link installation fails on this filesystem. Multicall binary
    # copies provide equivalent applet dispatch without privilege or PATH edits.
    applets = ['sh','bash','awk','cat','cut','dirname','basename','expr','grep',
               'head','tail','ln','mkdir','mv','rm','sed','sort','touch','tr',
               'uniq','wc','readlink','realpath','cp','sleep','date','xargs']
    for applet in applets:
        immutable_write(bindir/f'{applet}.exe', (bindir/'busybox.exe').read_bytes())
    env = os.environ.copy()
    env['PATH'] = str(bindir) + os.pathsep + env['PATH']
    work = RAW/'software/smoke'
    work.mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(RAW/'phage_proteins/ON857925.faa', 'fasta'))[:10]
    sample = ''.join(f'>p{i}\n{r.seq}\n' for i,r in enumerate(records))
    immutable_write(work/'sample.faa', sample.encode())
    commands = [
        ['version'],
        ['createdb','sample.faa','sample_db','-v','1'],
        ['cluster','sample_db','clusters','tmp','--min-seq-id','0.4','-c','0.8',
         '-s','7.5','--threads','2','-v','1'],
        ['createtsv','sample_db','sample_db','clusters','clusters.tsv','-v','1'],
    ]
    runs = []
    success = True
    for args in commands:
        try:
            result = subprocess.run([str(bindir/'mmseqs.exe'), *args], cwd=work,
                env=env, capture_output=True, text=True, timeout=180)
            item = {'arguments': args, 'returncode': result.returncode,
                'stdout': result.stdout, 'stderr': result.stderr}
            runs.append(item)
            print(json.dumps(item), flush=True)
            if result.returncode:
                success = False
                break
        except subprocess.TimeoutExpired:
            runs.append({'arguments': args, 'error': 'timeout after 180 seconds'})
            success = False
            break
    if success:
        lines = (work/'clusters.tsv').read_text().splitlines()
        success = {line.split('\t')[1] for line in lines} == {f'p{i}' for i in range(10)}
    report = {'complete': success, 'scope': 'ten-protein software smoke test only',
        'archive_sha256': sha(archive.read_bytes()),
        'executable_sha256': sha((bindir/'mmseqs.exe').read_bytes()),
        'source': source, 'runs': runs}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'runtime_check.json').write_text(json.dumps(report, indent=2))
    if not success:
        sys.exit(1)


if __name__ == '__main__':
    main()
