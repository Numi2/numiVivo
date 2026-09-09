#!/usr/bin/env python3
"""Rebuild source PCA, retain fixed-penalty trajectories, and run adaptive candidates."""
import argparse,hashlib,json,re,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','prepared','hagai','baron','frozen','protocol','out']:
 p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
protocol=json.loads(a.protocol.read_text());commands=[];checks=[]
def save(name,v):(a.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def run(label,args,ok=True):
 command=['/usr/bin/time','-l',str(a.binary),*map(str,args)]
 r=subprocess.run(command,capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
 commands.append(dict(label=label,command=command,exitCode=r.returncode,expectedSuccess=ok,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None));save('commands.json',commands)
 assert r.returncode==(0 if ok else 65),(label,r.stderr)
 print(label+' passed',flush=True);return r
for name,source in [('kang',a.prepared/'kang.h5ad'),('hagai',a.hagai),('baron',a.baron)]:
 root=a.out/name;root.mkdir()
 run(name+'-fit',['singlecell-h5ad-pca',source,'--plan',a.prepared/(name+'-fit.json'),'--output',root/'pca'])
 assert sha(root/'pca/scores.bin')==sha(a.frozen/name/'pca/scores.bin')
 assert json.loads((root/'pca/metadata.json').read_text())==json.loads((a.frozen/name/'pca/metadata.json').read_text())
 for mode in ['fixed','adaptive']:
  for seed in protocol['seeds'] if name!='baron' else [7]:
   if mode=='fixed':
    options=json.loads((a.prepared/'protocol.json').read_text())['integration']
   else:options=dict(protocol['nativeOptions'])
   options['seed']=seed;plan=root/f'{mode}-plan-{seed}.json'
   plan.write_text(json.dumps(dict(schemaVersion=1,inputKind='fitted',integration=options),indent=2)+'\n')
   target=root/(f'fixed-{seed}' if mode=='fixed' else f'seed-{seed}')
   r=run(f'{name}-{mode}-{seed}',['singlecell-pca-integrate',root/'pca','--plan',plan,'--output',target],name!='baron')
   if name=='baron':
    assert 'confounded with condition' in r.stderr and not target.exists()
    checks.append(dict(cohort=name,mode=mode,fullCohortConfoundingRejected=True));continue
   if mode=='fixed':
    prior=a.frozen/name/f'seed-{seed}'
    for file in ['scores.bin','memberships.bin','assignment-scores.bin','report.json']:
     assert sha(target/file)==sha(prior/file),(name,seed,file)
    checks.append(dict(cohort=name,seed=seed,allFixedMatrixAndReportBytesExact=True))
   run(f'{name}-{mode}-{seed}-verify',['singlecell-pca-integrate-verify',target])
   save('completed.json',checks)
save('checks.json',dict(status='passed',checks=checks,commands=len(commands),binarySHA256=sha(a.binary),protocolSHA256=sha(a.protocol),qualification='Fresh binary-bound inputs. Frozen outputs are byte comparison oracles only. Adaptive biological preservation requires separate evaluation.'))
