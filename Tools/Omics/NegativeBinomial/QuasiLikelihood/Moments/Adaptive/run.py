#!/usr/bin/env python3
"""Checkpoint complete native residual stages by frozen case and method."""
import argparse,concurrent.futures,gzip,hashlib,json,subprocess,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--binary',required=True);p.add_argument('--jobs',type=int,default=1);p.add_argument('--case');a=p.parse_args();assert 1<=a.jobs<=4
sha=lambda b:hashlib.sha256(b).hexdigest();protocol=json.loads((a.ql_root/'protocol.json').read_text());host=protocol['host'];binaryHash=subprocess.check_output(['ssh',host,'shasum','-a','256',a.binary]).decode().split()[0]
manifest=dict(sourceProtocolSHA256=sha((a.ql_root/'protocol.json').read_bytes()),protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md').read_bytes()),binarySHA256=binaryHash,host=host,binary=a.binary,cases=protocol['cases'])
a.out.mkdir(parents=True,exist_ok=True);mp=a.out/'protocol.json'
if mp.exists():assert json.loads(mp.read_text())==manifest
else:mp.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
tasks=[(c,m) for c in protocol['cases'] if a.case is None or c['id']==a.case for m in ['nativeTrend-adjusted','edgeRTrend-adjusted']];assert tasks

def execute(task):
 c,method=task;dest=a.out/c['id']/method;dest.mkdir(parents=True,exist_ok=True)
 if (dest/'run.json').exists():
  receipt=json.loads((dest/'run.json').read_text());assert receipt['binarySHA256']==binaryHash and sha((dest/'native.json.gz').read_bytes())==receipt['outputSHA256'];return receipt
 source=a.ql_root/c['id'];raw=(source/'input.json.gz').read_bytes();ir=json.loads((source/'input-receipt.json').read_text());assert sha(raw)==ir['inputSHA256'];i=json.loads(gzip.decompress(raw))
 rr=json.loads((source/'run.json').read_text());raw=(source/'reference.json.gz').read_bytes();assert sha(raw)==rr['outputSHA256'];reference=json.loads(gzip.decompress(raw))['results'][method]
 rows=[dict(counts=y,means=mu,design=i['design'],dispersion=t['actualDispersion'],averageQuasiDispersion=t['averageQLScale']) for y,mu,t in zip(i['counts'],reference['means'],reference['table'],strict=True)]
 payload=json.dumps(dict(residuals=rows),allow_nan=False,separators=(',',':')).encode();t0=time.monotonic()
 print('start',c['id'],method,len(rows),flush=True)
 with (dest/'stderr.log').open('wb') as log:
  process=subprocess.run(['ssh',host,'/usr/bin/time','-l',a.binary],input=payload,stdout=subprocess.PIPE,stderr=log)
 elapsed=time.monotonic()-t0
 if process.returncode!=0:
  (dest/'failed-stdout.txt').write_bytes(process.stdout);(dest/'failure.json').write_text(json.dumps(dict(exitCode=process.returncode,seconds=elapsed,binarySHA256=binaryHash),indent=2)+'\n');raise RuntimeError((c['id'],method,process.returncode))
 raw=process.stdout;packed=gzip.compress(raw,mtime=0);(dest/'native.json.gz').write_bytes(packed);output=json.loads(raw)['residuals'];assert len(output)==len(rows)
 details=[];failures=[];terms=0;moments=0;adaptiveMoments=0;maxLeverage=0.;maximumBound=0.
 for j,(out,ref) in enumerate(zip(output,reference['table'],strict=True)):
  if 'error' in out:failures.append(dict(featureIndex=i['featureIndices'][j],error=out['error']));continue
  v=out['value'];mm=v['moments'];terms+=sum(x['evaluatedCounts'] for x in mm);moments+=len(mm)
  adaptiveMoments+=sum('summation' in x for x in mm)
  bound=max(max((x['meanTruncationBound']+x.get('summation',{}).get('meanErrorBound',0))/x['mean'],(x['varianceTruncationBound']+x.get('summation',{}).get('varianceErrorBound',0))/x['variance']) for x in mm);maximumBound=max(maximumBound,bound)
  leverage=float(np.max(np.abs(np.asarray(v['leverage'])-reference['leverage'][j])));maxLeverage=max(maxLeverage,leverage)
  if leverage>1e-9 or bound>1e-10:failures.append(dict(featureIndex=i['featureIndices'][j],leverageError=leverage,maximumRelativeBound=bound))
  relative=lambda x,y:abs(x-y)/max(abs(y),1e-30)
  details.append(dict(featureIndex=i['featureIndices'][j],relativeDevianceDifference=relative(v['deviance'],ref['residualDeviance']),relativeDFDifference=relative(v['degreesOfFreedom'],ref['residualDF']),maximumLeverageError=leverage,evaluatedCounts=sum(x['evaluatedCounts'] for x in mm)))
 receipt=dict(case=c['id'],method=method,status='passed-conditional-stage' if not failures else 'completed-with-failures',genes=len(rows),moments=moments,adaptiveMoments=adaptiveMoments,evaluatedCounts=terms,seconds=elapsed,maximumLeverageError=maxLeverage,maximumRelativeCombinedBound=maximumBound,failures=failures,relativeDevianceDifferenceQuantiles=np.quantile([v['relativeDevianceDifference'] for v in details],[0,.5,.9,.99,1]).tolist(),relativeDFDifferenceQuantiles=np.quantile([v['relativeDFDifference'] for v in details],[0,.5,.9,.99,1]).tolist(),inputSHA256=ir['inputSHA256'],referenceSHA256=rr['outputSHA256'],requestSHA256=sha(payload),outputSHA256=sha(packed),logicalSHA256=sha(raw),binarySHA256=binaryHash)
 (dest/'gene-checks.json.gz').write_bytes(gzip.compress(json.dumps(details,separators=(',',':')).encode(),mtime=0));(dest/'run.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='failures'}),flush=True);return receipt
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as ex:results=list(ex.map(execute,tasks))
(a.out/('selected-complete.json' if a.case else 'complete.json')).write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
