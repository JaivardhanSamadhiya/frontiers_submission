"""Training-only protein-family representation for a GenoPHI-derived benchmark.

This is a Windows-compatible adapter, not the complete upstream workflow.
See GENOPHI_TRAINING_AMENDMENT_2026_10_05.md for the scientific differences.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import time
import shutil
import zipfile

import numpy as np
import pandas as pd
from Bio import SeqIO


def digest(data):
    return hashlib.sha256(data).hexdigest()


def merged_features(raw, fitting_rows, groups=None):
    """Filter and merge using fitting rows only; OR merged members at inference."""
    a=np.asarray(raw,dtype=np.uint8)
    fit=np.asarray(fitting_rows,dtype=int)
    if a.ndim!=2 or not len(fit):
        raise ValueError('Expected a nonempty training partition')
    if groups is None:
        support=a[fit].sum(axis=0)
    else:
        g=np.asarray(groups)
        support=sum(a[fit[g[fit]==v]].max(axis=0) for v in np.unique(g[fit]))
    keep=np.flatnonzero(support>=2)
    buckets={}
    for c in keep:
        key=a[fit,c].tobytes()  # exact pattern, not only a hash
        buckets.setdefault(key,[]).append(int(c))
    members=list(buckets.values())
    if not members:
        raise ValueError('No shared protein-family features remain')
    return np.column_stack([a[:,cs].max(axis=1) for cs in members]),members


def top_three(scores):
    a=np.asarray(scores)
    if len(a)<3 or not np.isfinite(a).all():
        raise ValueError('Need at least three finite scores')
    return np.argsort(-a,kind='stable')[:3]


def unique_queries(records):
    """Deduplicate identical query strings, retaining every original identifier."""
    sequences={};aliases={};unique=[];seen=set()
    for identifier,sequence in records:
        if identifier in seen:raise ValueError('Duplicate query identifier')
        seen.add(identifier)
        if sequence not in sequences:
            sequences[sequence]=identifier;unique.append((identifier,sequence));aliases[identifier]=[]
        aliases[sequences[sequence]].append(identifier)
    return unique,aliases


class ProteinReferences:
    def __init__(self,root):
        self.root=Path(root)
        self.raw=self.root/'data/external_validation/genophi_comparison'
        self.cache=self.raw/'feature_cache_v2'
        self.cache.mkdir(parents=True,exist_ok=True)
        # Native Cygwin/MMseqs temporary files live on the internal NTFS volume,
        # not the exFAT project drive. Only verified compact caches persist there.
        # Avoid Microsoft Store Python AppData redirection: native child
        # processes do not necessarily share the redirected filesystem view.
        self.scratch=Path(os.environ['USERPROFILE'])/'.cache/precisionphage/genophi_native_v2'
        self.scratch.mkdir(parents=True,exist_ok=True)
        archive=self.raw/'software/mmseqs-win64-15-6f452.zip'
        recorded=json.loads((self.root/'data/results_v2/genophi_comparison/runtime_check.json').read_text())
        if digest(archive.read_bytes())!=recorded['archive_sha256']:
            raise ValueError('Native runtime archive changed')
        with zipfile.ZipFile(archive) as z:
            for member in z.infolist():
                if member.is_dir():continue
                target=(self.scratch/member.filename).resolve()
                if not target.is_relative_to(self.scratch.resolve()):raise ValueError('Unsafe runtime member')
                content=z.read(member)
                target.parent.mkdir(parents=True,exist_ok=True)
                if target.exists() and target.read_bytes()!=content:raise ValueError('Runtime file changed')
                if not target.exists():target.write_bytes(content)
        self.bin=self.scratch/'mmseqs/bin'
        for applet in ['sh','bash','awk','cat','cut','dirname','basename','expr','grep','head',
                'tail','ln','mkdir','mv','rm','sed','sort','touch','tr','uniq','wc','readlink',
                'realpath','cp','sleep','date','xargs']:
            p=self.bin/f'{applet}.exe';content=(self.bin/'busybox.exe').read_bytes()
            if p.exists() and p.read_bytes()!=content:raise ValueError('Runtime applet changed')
            if not p.exists():p.write_bytes(content)
        self.exe=self.bin/'mmseqs.exe'
        self.env=os.environ.copy()
        self.env['PATH']=str(self.bin)+os.pathsep+self.env['PATH']
        manifest=json.loads((self.root/'data/results_v2/genophi_comparison/input_manifest.json').read_text())
        self.entries={e['entity']:e for e in manifest['proteomes']}
        self.sequences={}
        for name,e in self.entries.items():
            path=self.root/e['path']
            if digest(path.read_bytes())!=e['sha256']:
                raise ValueError(f'Proteome checksum changed: {name}')
            self.sequences[name]=[(f'{name}|p{i:06d}',str(r.seq))
                                  for i,r in enumerate(SeqIO.parse(path,'fasta'))]

    def command(self,args,work):
        command=[str(self.exe),*map(str,args)]
        started=time.time()
        with (work/'commands.log').open('a',encoding='utf-8') as log:
            log.write(json.dumps(command)+'\n')
        # Do not share a Python-owned append handle with the Cygwin subprocess.
        result=subprocess.run(command,cwd=work,env=self.env,capture_output=True)
        with (work/'commands.log').open('a',encoding='utf-8') as log:
            log.write(result.stdout.decode('utf-8',errors='replace'))
            log.write(result.stderr.decode('utf-8',errors='replace'))
            log.write(f'EXIT {result.returncode}, seconds {time.time()-started:.3f}\n')
        if result.returncode:
            raise RuntimeError(f'MMseqs2 {args[0]} failed; see {work / "commands.log"}')

    def write_fasta(self,names,path):
        with path.open('w',encoding='ascii',newline='\n') as out:
            for name in names:
                for identifier,sequence in self.sequences[name]:
                    out.write(f'>{identifier}\n{sequence}\n')

    def build(self,training,all_names,kind):
        """Cluster only training entities; all other entities are search queries."""
        training=sorted(training); all_names=list(all_names)
        if not set(training)<=set(all_names): raise ValueError('Unknown training entity')
        spec={'schema':2,'deduplicate_identical_queries':True,'training':training,'all_names':all_names,'kind':kind,
            'sequences':{n:self.entries[n]['sha256'] for n in all_names},
            'mmseqs_sha256':digest(self.exe.read_bytes()),'identity':0.4,'coverage':0.8,'sensitivity':7.5}
        key=digest(json.dumps(spec,sort_keys=True).encode())[:20]
        saved=self.cache/f'{kind}_{key}'; saved.mkdir(parents=True,exist_ok=True)
        if (saved/'complete.json').exists():
            complete=json.loads((saved/'complete.json').read_text())
            if complete['spec']!=spec or digest((saved/'raw_features.npz').read_bytes())!=complete['features_sha256']:
                raise ValueError('Feature cache provenance mismatch')
            with np.load(saved/'raw_features.npz',allow_pickle=False) as z:
                return z['presence'],z['clusters'].tolist(),str(saved.relative_to(self.root))
        if shutil.disk_usage(self.scratch).free<2*1024**3:
            raise RuntimeError('Less than 2 GiB free on internal scratch volume; stop safely')
        work=self.scratch/'work'/f'{kind}_{key}';work.mkdir(parents=True,exist_ok=True)
        print(f'Building {kind} references: {len(training)} training entities ({key})',flush=True)
        self.write_fasta(training,work/'reference.faa')
        self.command(['createdb','reference.faa','reference','-v','1'],work)
        self.command(['cluster','reference','clusters','tmp_cluster','--min-seq-id','0.4',
            '-c','0.8','-s','7.5','--threads','4','--split-memory-limit','4G','-v','1'],work)
        self.command(['createtsv','reference','reference','clusters','clusters.tsv','-v','1'],work)
        table=pd.read_csv(work/'clusters.tsv',sep='\t',header=None,names=['cluster','protein'],dtype=str)
        expected={i for n in training for i,s in self.sequences[n]}
        if set(table.protein)!=expected or table.protein.duplicated().any():
            raise ValueError('Clustering lost or duplicated training proteins')
        protein_to_cluster=dict(zip(table.protein,table.cluster))
        clusters=sorted(table.cluster.unique()); cidx={c:i for i,c in enumerate(clusters)}
        nidx={n:i for i,n in enumerate(all_names)}
        a=np.zeros((len(all_names),len(clusters)),dtype=np.uint8)
        for row in table.itertuples(index=False):
            a[nidx[row.protein.split('|')[0]],cidx[row.cluster]]=1
        query=[n for n in all_names if n not in training]
        matches=0;unique_count=0
        if query:
            unique,aliases=unique_queries([record for n in query for record in self.sequences[n]])
            unique_count=len(unique)
            with (work/'query.faa').open('w',encoding='ascii',newline='\n') as out:
                for identifier,sequence in unique:out.write(f'>{identifier}\n{sequence}\n')
            (saved/'query_aliases.json').write_text(json.dumps(aliases))
            self.command(['createdb','query.faa','query','-v','1'],work)
            self.command(['search','query','reference','matches','tmp_search','--min-seq-id','0.4',
                '-c','0.8','-s','7.5','--threads','4','--split-memory-limit','4G','-v','1'],work)
            self.command(['convertalis','query','reference','matches','matches.tsv',
                '--format-output','query,target,bits,evalue','-v','1'],work)
            if (work/'matches.tsv').stat().st_size:
                hits=pd.read_csv(work/'matches.tsv',sep='\t',header=None,names=['query','target','bits','evalue'])
                hits=hits.sort_values(['query','bits','evalue','target'],ascending=[True,False,True,True]).drop_duplicates('query')
                for row in hits.itertuples(index=False):
                    if row.query not in aliases or row.target not in protein_to_cluster:
                        raise ValueError('Query/reference membership mismatch')
                    for original in aliases[row.query]:
                        name=original.split('|')[0]
                        if name not in query:raise ValueError('Query assigned outside query hosts')
                        a[nidx[name],cidx[protein_to_cluster[row.target]]]=1;matches+=1
        np.savez_compressed(work/'raw_features.npz',presence=a,clusters=np.asarray(clusters),names=np.asarray(all_names))
        complete={'spec':spec,'features_sha256':digest((work/'raw_features.npz').read_bytes()),
            'training_proteins':len(expected),'query_proteins':sum(len(self.sequences[n]) for n in query),
            'unique_query_sequences':unique_count,'query_proteins_matched':matches,'clusters':len(clusters)}
        for name in ['raw_features.npz','commands.log','clusters.tsv','matches.tsv']:
            if (work/name).exists():shutil.copyfile(work/name,saved/name)
        if digest((saved/'raw_features.npz').read_bytes())!=complete['features_sha256']:
            raise ValueError('Cache copy checksum mismatch; scratch retained')
        (saved/'complete.json').write_text(json.dumps(complete,indent=2))
        # Delete only this validated, task-owned, reproducible native scratch dir.
        resolved=work.resolve();allowed=(self.scratch/'work').resolve()
        if resolved.parent!=allowed or not resolved.name.startswith(f'{kind}_'):
            raise ValueError('Unsafe scratch cleanup path')
        shutil.rmtree(resolved)
        return a,clusters,str(saved.relative_to(self.root))
