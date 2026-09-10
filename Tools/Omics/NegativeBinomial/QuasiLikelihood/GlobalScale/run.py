#!/usr/bin/env python3
"""Checkpoint native global scale/refit families and preserve all failures."""
import argparse,concurrent.futures,gzip,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--reference-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',required=True);p.add_argument('--jobs',type=int,default=1);p.add_argument('--case');a=p.parse_args();assert 1<=a.jobs<=3
sha=lambda b:hashlib.sha256(b).hexdigest();protocol=json.loads((a.ql_root/'protocol.json').read_text());binarySHA=subprocess.check_output(['ssh',protocol['host'],'shasum','-a','256',a.binary]).decode().split()[0]
manifest=dict(cases=protocol['cases'],host=protocol['host'],binary=a.binary,binarySHA256=binarySHA,protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md').read_bytes()),priorProtocolSHA256=sha((a.ql_root/'protocol.json').read_bytes()))
a.out.mkdir(parents=True,exist_ok=True)
if (a.out/'protocol.json').exists():assert json.loads((a.out/'protocol.json').read_text())==manifest
else:(a.out/'protocol.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
def execute(task):
 c,method=task;dest=a.out/c['id']/method;dest.mkdir(parents=True,exist_ok=True)
 source=a.ql_root/c['id'];rr=a.reference_root/c['id'];ir=json.loads((rr/'reference-receipt.json').read_text());assert ir['exitCode']==0
 raw=(source/'input.json.gz').read_bytes();assert sha(raw)==ir['inputSHA256'];i=json.loads(gzip.decompress(raw))
 raw=(rr/'reference.json.gz').read_bytes();assert sha(raw)==ir['outputSHA256'];reference=json.loads(gzip.decompress(raw))['results'][method]
 payload=json.dumps(dict(counts=i['counts'],design=i['design'],offsets=i['offsets'],contrast=i['contrast'],trendDispersions=reference['trendDispersions'],abundanceCovariates=reference['abundanceCovariates']),separators=(',',':'),allow_nan=False).encode()
 if (dest/'run.json').exists():
  r=json.loads((dest/'run.json').read_text());assert r['binarySHA256']==binarySHA and r['requestSHA256']==sha(payload) and r['outputSHA256']==sha((dest/'native.json.gz').read_bytes());return r
 print('start',c['id'],method,len(i['counts']),flush=True);t0=time.monotonic()
 with (dest/'native.log').open('wb') as log:
  run=subprocess.run(['ssh',protocol['host'],'/usr/bin/time','-l',a.binary],input=payload,stdout=subprocess.PIPE,stderr=log)
 elapsed=time.monotonic()-t0
 receipt=dict(case=c['id'],method=method,exitCode=run.returncode,seconds=elapsed,binarySHA256=binarySHA,inputSHA256=ir['inputSHA256'],referenceSHA256=ir['outputSHA256'],requestSHA256=sha(payload))
 packed=gzip.compress(run.stdout,mtime=0);(dest/'native.json.gz').write_bytes(packed);receipt.update(outputSHA256=sha(packed),logicalSHA256=sha(run.stdout))
 if run.returncode!=0:receipt.update(status='failed-process')
 else:
  output=json.loads(run.stdout)
  if 'error' in output:receipt.update(status='failed-input',error=output['error'])
  else:
   result=output['fit'];receipt.update(status='passed-native-stage' if result['completed'] else 'incomplete-native-stage',genes=len(i['counts']),failures=result['failures'])
   if result['completed']:
    relative=lambda x,y:float(np.max(np.abs(np.asarray(x)-np.asarray(y))/np.maximum(1,np.abs(y))))
    receipt.update(averageQuasiDispersion=result['averageQuasiDispersion'],referenceScale=reference['averageQuasiDispersion'],scaleRelativeDifference=abs(result['averageQuasiDispersion']/reference['averageQuasiDispersion']-1),initialMeanRelativeDifference=relative([f['means'] for f in result['initialFits']],reference['initialMeans']),finalMeanRelativeDifference=relative([f['means'] for f in result['refittedFits']],reference['finalMeans']),maximumScaledScore=max(f['maximumScaledScore'] for f in result['initialFits']+result['refittedFits']),updateScales=[u['outputScale'] for u in result['updates']],omittedIndices=[u['omittedIndices'] for u in result['updates']],momentEvaluations=sum(u['momentEvaluations'] for u in result['updates'])+sum(m['evaluatedCounts'] for residual in result['adjustedResiduals'] for m in residual['moments']))
 (dest/'run.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt),flush=True);return receipt
jobs=[(c,m) for c in protocol['cases'] if a.case is None or c['id']==a.case for m in ['nativeTrend-adjusted','edgeRTrend-adjusted']];assert jobs
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as executor:results=list(executor.map(execute,jobs))
(a.out/('pilot-complete.json' if a.case else 'complete.json')).write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
