#!/usr/bin/env python3
"""Preserve compact complete null evidence and hashes of detailed native models."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def sha(b):return hashlib.sha256(b).hexdigest()
results=json.loads((a.results/'results.json').read_text());assert not results['pending']
complete=json.loads((a.root/'native-complete.json').read_text());assert complete['status']=='finished-all-declared-attempts'
files=[];external=[];remote=Path('/Users/n/numivivo-null-benchmark-20260910')
inventory=json.loads((a.root/'native-report-inventory.json').read_text())
for name in ['preparation-check.json','native-runs.json','native-complete.json','reference-runs.json','native-run.log','reference-run.log','prepare-kang.log','prepare-hagai.log','duplicate-cleanup.json','cow-consolidation.json']:
 files.append(a.root/name)
for name in ['native-report-inventory.json','local-model-mirror-cleanup.json','local-annotation-mirror-cleanup.json','restoration.log','restoration-verify.json','restoration-verify.log','restoration-integrity.log','final-check.log','null-hit-audit.log','attempts.json']:
 files.append(a.root/name)
files.append(a.root/'restoration-integrity/checks.json')
files+=list((a.root/'restoration-integrity').glob('*/rejection.log'))
for study in ['kang','hagai']:
 d=a.root/study
 for name in ['preparation.json','native-input.json','annotation-plan.json','annotation.log','annotation-check.json','annotation-check.log']:
  files.append(d/name)
 for seed in range(1,11):
  s=d/str(seed)
  for name in ['reference.npz','memberships.json','default-plan.json','active-plan.json','default-publish.log','active-publish.log','default-verify.log','active-verify.log','r-run.log']:
   if (s/name).exists():files.append(s/name)
  for sub in ['r-input','r-reference']:files+=sorted(f for f in (s/sub).iterdir() if f.is_file())
  for mode in ['default','active']:
   m=s/mode
   if not (m/'archive.json').exists():continue
   files+=[m/'archive.json',m/'plan.json.gz',m/'receipt.json.gz']
   ar=json.loads((m/'archive.json').read_text())
   for e in ar['entries']:
    if e['logicalPath']!='report.json':continue
    info=inventory[str((m/e['path']).relative_to(a.root))];assert info['sha256']==e['sha256']
    external.append(dict(path=str(remote/(m.relative_to(a.root))/e['path']),bytes=info['bytes'],sha256=e['sha256'],logicalSHA256=e['logicalSHA256'],kind='complete-native-model-report'))
 for name in ['original.h5ad','annotated.h5ad']:
  if name=='original.h5ad':
   manifest=json.loads((d/'preparation.json').read_text());digest=manifest['sourceSHA256']
  else:
   manifest=json.loads((d/'annotation-check.json').read_text());digest=manifest['annotatedSHA256']
  external.append(dict(path=str(remote/study/name),sha256=digest,kind='retained-complete-H5AD'))
entries=[]
for f in sorted(set(files)):
 raw=f.read_bytes();rel=f.relative_to(a.root).as_posix();stored=raw if f.suffix in ['.gz','.npz'] else gzip.compress(raw,mtime=0)
 path=rel if f.suffix in ['.gz','.npz'] else rel+'.gz';out=a.out/path;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(stored)
 entries.append(dict(path=path,bytes=len(stored),sha256=sha(stored),logicalPath=rel,logicalBytes=len(raw),logicalSHA256=sha(raw)))
for f in sorted(a.results.iterdir()):
 if not f.is_file():continue
 raw=f.read_bytes();rel='results/'+f.name;stored=raw if f.suffix=='.gz' else gzip.compress(raw,mtime=0)
 path=rel if f.suffix=='.gz' else rel+'.gz';out=a.out/path;out.parent.mkdir(parents=True,exist_ok=True);out.write_bytes(stored)
 entries.append(dict(path=path,bytes=len(stored),sha256=sha(stored),logicalPath=rel,logicalBytes=len(raw),logicalSHA256=sha(raw)))
manifest=dict(schemaVersion=1,status=results['status'],entries=entries,external=external,
 qualification='All declared sham split attempts preserved; overlapping untreated-cell resamples do not establish general FDR calibration or power. Detailed native model parameters remain at hashed remote paths.')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
for e in entries:
 b=(a.out/e['path']).read_bytes();assert sha(b)==e['sha256'] and len(b)==e['bytes']
 decoded=gzip.decompress(b) if e['path']!=e['logicalPath'] else b
 assert sha(decoded)==e['logicalSHA256'] and len(decoded)==e['logicalBytes']
print(json.dumps(dict(status='archived-and-verified',files=len(entries),storedBytes=sum(e['bytes'] for e in entries),external=len(external))))
