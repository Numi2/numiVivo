#!/usr/bin/env python3
"""Native full-cohort integration, immutable replay, frozen trajectory and confounding rejection."""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
import struct
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','oracle','prepared','hagai','baron','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[];checks=[]
protocol=json.loads((a.prepared/'protocol.json').read_text())
def save(path,v):path.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def run(label,binary,args,ok=True):
 command=['/usr/bin/time','-l',str(binary),*map(str,args)];r=subprocess.run(command,capture_output=True,text=True)
 (a.out/(label+'.log')).write_text(r.stdout+r.stderr);rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
 commands.append(dict(label=label,command=command,exitCode=r.returncode,expectedSuccess=ok,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None));save(a.out/'commands.json',commands)
 assert r.returncode==(0 if ok else 65),(label,r.stderr)
 print(label+' passed',flush=True);return r
for name,source in [('kang',a.prepared/'kang.h5ad'),('hagai',a.hagai),('baron',a.baron)]:
 root=a.out/name;root.mkdir()
 run(name+'-fit',a.binary,['singlecell-h5ad-pca',source,'--plan',a.prepared/(name+'-fit.json'),'--output',root/'pca'])
 dimensions=json.loads((root/'pca/model.json').read_text())['options']['components']
 for seed in protocol['seeds'] if name!='baron' else [7]:
  options=dict(protocol['integration'],seed=seed);plan=root/f'plan-{seed}.json';save(plan,dict(schemaVersion=1,inputKind='fitted',integration=options));opts=root/f'options-{seed}.json';save(opts,options)
  target=root/f'seed-{seed}'
  result=run(name+f'-{seed}-integrate',a.binary,['singlecell-pca-integrate',root/'pca','--plan',plan,'--output',target],name!='baron')
  if name=='baron':
   assert 'confounded with condition' in result.stderr and not target.exists();checks.append(dict(cohort=name,fullCohortConfoundingRejected=True));break
  oracle=root/f'legacy-{seed}.json';run(name+f'-{seed}-legacy',a.oracle,[root/'pca',dimensions,opts,oracle]);legacy=json.loads(oracle.read_text());report=json.loads((target/'report.json').read_text());metadata=json.loads((target/'metadata.json').read_text());assert legacy['cells']==[{k:c[k] for k in ['sampleID','barcode']} for c in metadata['cells']]
  for file,key,width in [('scores.bin','scores',dimensions),('memberships.bin','memberships',options['clusters']),('assignment-scores.bin','assignmentScores',dimensions)]:
   assert (target/file).stat().st_size==len(metadata['cells'])*width*16
   with (target/file).open('rb') as stream:
    for i,row in enumerate(legacy[key]):
     actual=stream.read(width*16)
     expected=b''.join(struct.pack('<IId',i,j,v) for j,v in enumerate(row))
     assert actual==expected,(name,seed,key,i)
    assert stream.read(1)==b''
  for key in ['method','levels','cellLevels','assignmentCenters','objectives','relativeImprovements','stoppingReason','maximumRidgeResidual']:assert report[key]==legacy[key],key
  checks.append(dict(cohort=name,seed=seed,cells=len(metadata['cells']),allLegacyMatrixBitsExact=True,allLegacyDiagnosticsExact=True));save(a.out/'completed.json',checks)
  run(name+f'-{seed}-verify',a.binary,['singlecell-pca-integrate-verify',target])
  # Large legacy JSON is oracle evidence; keep a lossless deterministic archive.
  import gzip
  with oracle.with_suffix('.json.gz').open('wb') as out:
   with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as g:g.write(oracle.read_bytes())
  oracle.unlink()
save(a.out/'checks.json',dict(status='passed',commands=len(commands),checks=checks,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),oracleSHA256=hashlib.sha256(a.oracle.read_bytes()).hexdigest(),qualification='Full-cohort source reconstruction and numerical trajectory; biological preservation requires independent evaluation.'))
