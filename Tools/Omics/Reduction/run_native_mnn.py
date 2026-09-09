#!/usr/bin/env python3
"""Fresh binary-bound source PCA and full-cohort MNN runs on the owning Mac mini."""
import argparse,hashlib,json,re,shlex,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--host',required=True);p.add_argument('--resume',action='store_true',help='Verify completed transfers and continue an interrupted run; preserve every command attempt.')
for name in ['binary','hdf5','root','spec','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=a.resume);spec=json.loads(a.spec.read_text());commands=json.loads((a.out/'commands.json').read_text()) if a.resume else [];artifacts=json.loads((a.out/'artifacts.json').read_text()) if a.resume else []
def write(name,v):(a.out/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def remote(args):return subprocess.run(['ssh',a.host,shlex.join(list(map(str,args)))],capture_output=True,text=True)
def run(label,args,ok=True):
 attempt=sum(c.get('baseLabel',c['label'])==label for c in commands)+1;base=label;label=label if attempt==1 else label+'-attempt-'+str(attempt)
 command=['env','NUMIVIVO_HDF5_LIBRARY='+str(a.hdf5),'/usr/bin/time','-l',a.binary,*args];r=remote(command);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
 rss=re.search(r'(\d+)\s+maximum resident set size',r.stderr);seconds=re.search(r'([\d.]+) real',r.stderr)
 commands.append(dict(label=label,baseLabel=base,command=list(map(str,command)),exitCode=r.returncode,expectedSuccess=ok,seconds=float(seconds[1]) if seconds else None,maximumResidentBytes=int(rss[1]) if rss else None));write('commands.json',commands);assert r.returncode==(0 if ok else 65),(label,r.stderr);print(label+' passed',flush=True);return r
result=remote(['shasum','-a','256',a.binary]);assert result.returncode==0 and result.stdout.split()[0]==spec['binarySHA256']
result=remote(['test','-d',a.root] if a.resume else ['mkdir',a.root]);assert result.returncode==0,result.stderr
for cohort in spec['cohorts']:
 name=cohort['name'];root=a.root/name
 completed=next((x for x in artifacts if x['cohort']==name),None)
 if completed:
  assert json.loads((a.out/(name+'-fit.json')).read_text())==cohort['fit'],'Completed PCA plan changed'
  assert json.loads((a.out/(name+'-mnn.json')).read_text())==cohort['plan'],'Completed MNN plan changed'
  for e in completed['files']:
   if e['path'].endswith('.h5ad'):assert e['sha256']==cohort['sourceSHA256'],'Completed source changed'
   result=remote(['shasum','-a','256',root/e['path']]);assert result.returncode==0 and result.stdout.split()[0]==e['sha256']
   if not e['path'].endswith('.h5ad'):
    path=a.out/name/e['path'];assert path.stat().st_size==e['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==e['sha256']
  print(name+' completed artifacts rechecked',flush=True);continue
 exists=remote(['test','-d',root]).returncode==0
 assert a.resume or not exists
 if not exists:
  result=remote(['mkdir',root]);assert result.returncode==0
 assert remote(['test','!','-e',root/'mnn']).returncode==0,'Untransferred published output requires manual inspection'
 for filename,key in [('fit.json','fit'),('mnn.json','plan')]:
  path=a.out/(name+'-'+filename)
  if path.exists():
   prior=json.loads(path.read_text())
   if key=='fit':assert prior==cohort[key],'Cannot resume with changed PCA plan'
   elif prior!=cohort[key]:
    before=json.loads(json.dumps(prior));after=json.loads(json.dumps(cohort[key]));before['mnn'].pop('maximumWork',None);after['mnn'].pop('maximumWork',None);assert before==after,'Only work admission can change in a resumed plan'
    old=a.out/(name+'-mnn-before-resume-'+str(len(commands))+'.json');assert not old.exists();old.write_bytes(path.read_bytes())
  path.write_text(json.dumps(cohort[key],indent=2)+'\n');subprocess.run(['rsync',str(path),a.host+':'+str(root/filename)],check=True)
 check=remote(['shasum','-a','256',cohort['source']]);assert check.returncode==0 and check.stdout.split()[0]==cohort['sourceSHA256']
 if not any(c.get('baseLabel',c['label'])==name+'-fit' and c['exitCode']==0 for c in commands):run(name+'-fit',['singlecell-h5ad-pca',cohort['source'],'--plan',root/'fit.json','--output',root/'pca'])
 r=run(name+'-mnn',['singlecell-pca-integrate',root/'pca','--plan',root/'mnn.json','--output',root/'mnn'],cohort.get('expectedSuccess',True))
 if cohort.get('expectedSuccess',True):run(name+'-verify',['singlecell-pca-integrate-verify',root/'mnn'])
 else:
  assert 'confounded with condition' in r.stderr
  check=remote(['test','!','-e',root/'mnn']);assert check.returncode==0
 local=a.out/name;local.mkdir();subprocess.run(['rsync','-az','--exclude','*.h5ad',a.host+':'+str(root)+'/',str(local)+'/'],check=True)
 script='''import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);entries=[]
for p in sorted(root.rglob('*')):
 if not p.is_file():continue
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 entries.append(dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=h.hexdigest()))
print(json.dumps(entries))
'''
 result=remote(['python3','-c',script,root]);assert result.returncode==0;entries=json.loads(result.stdout)
 for e in entries:
  if e['path'].endswith('.h5ad'):assert e['sha256']==cohort['sourceSHA256'];continue
  path=local/e['path'];assert path.stat().st_size==e['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==e['sha256']
 artifacts.append(dict(cohort=name,files=entries,allTransferredBytesExact=True));write('artifacts.json',artifacts)
write('checks.json',dict(status='passed',binarySHA256=spec['binarySHA256'],commands=len(commands),priorFailedAttempts=[c['label'] for c in commands if c['exitCode']!=(0 if c['expectedSuccess'] else 65)],allPublishedNativeRunsReplayed=True,allTransferredBytesExact=True,qualification='Fresh native PCA and MNN artifacts, full source cohorts. All failed attempts remain recorded. Confounding rejection and numerical replay are separate from biological evaluation. Canonical H5AD inputs remain external, with every omitted alias hash checked.'))
