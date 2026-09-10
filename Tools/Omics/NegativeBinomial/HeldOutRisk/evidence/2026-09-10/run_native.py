#!/usr/bin/env python3
"""Run disjoint fold processes, freezing models before opening test counts for scoring."""
import argparse,concurrent.futures,gzip,hashlib,json,subprocess,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--jobs',type=int,default=4)
a=p.parse_args();assert 1<=a.jobs<=4
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((a.root/'protocol.json').read_text());binary_sha=sha(a.binary)
assert sha(a.root/'metadata.json.gz')==protocol['metadataSHA256']
def execute(rec):
 d=a.root/rec['id'];assert sha(d/'training.json.gz')==rec['trainingSHA256'];records=[]
 for stage,inputs in [('fit',[a.root/'metadata.json.gz',d/'training.json.gz']),('score',[d/'model.json.gz',d/'test.json.gz'])]:
  if stage=='score':assert sha(d/'test.json.gz')==rec['testSHA256']
  assert sha(a.binary)==binary_sha
  target=d/('model.json.gz' if stage=='fit' else 'scores.json.gz');assert not target.exists()
  start=time.monotonic();command=['/usr/bin/time','-l',str(a.binary),stage]+list(map(str,inputs))
  with (d/(stage+'.log')).open('wb') as f:result=subprocess.run(command,stdout=subprocess.PIPE,stderr=f)
  target.write_bytes(gzip.compress(result.stdout,mtime=0))
  record=dict(case=rec['id'],stage=stage,exitCode=result.returncode,seconds=time.monotonic()-start,command=command,
   outputSHA256=sha(target),logicalSHA256=hashlib.sha256(result.stdout).hexdigest(),outputBytes=len(result.stdout))
  records.append(record)
  (d/'runs.json').write_text(json.dumps(records,sort_keys=True,indent=2)+'\n')
  if result.returncode:break
 return records
records=[];start=time.monotonic()
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
 futures={pool.submit(execute,r):r['id'] for r in protocol['cases']}
 for future in concurrent.futures.as_completed(futures):
  result=future.result();records.extend(result)
  summary=dict(binarySHA256=binary_sha,protocolSHA256=sha(a.root/'protocol.json'),jobs=a.jobs,runs=sorted(records,key=lambda r:(r['case'],r['stage'])))
  (a.root/'native-runs.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
  print(json.dumps(dict(case=futures[future],stages=len(result),exitCodes=[r['exitCode'] for r in result])),flush=True)
(a.root/'native-complete.json').write_text(json.dumps(dict(cases=len(protocol['cases']),seconds=time.monotonic()-start,
 fitSuccesses=sum(r['stage']=='fit' and r['exitCode']==0 for r in records),scoreSuccesses=sum(r['stage']=='score' and r['exitCode']==0 for r in records),
 failedCommands=[r for r in records if r['exitCode']],binarySHA256=binary_sha),sort_keys=True,indent=2)+'\n')
if any(r['exitCode'] for r in records):raise SystemExit(1)
