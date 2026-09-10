#!/usr/bin/env python3
"""Check terminal native families with native-input LOWESS and fixed-scale tight GLMs."""
import argparse,concurrent.futures,gzip,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--ql-root',type=Path,required=True);p.add_argument('--reference-root',type=Path,required=True);p.add_argument('--root',type=Path,required=True);p.add_argument('--pilot',action='store_true');p.add_argument('--jobs',type=int,default=1);a=p.parse_args();assert 1<=a.jobs<=2
sha=lambda b:hashlib.sha256(b).hexdigest();script=Path(__file__).with_suffix('.R');runs=json.loads((a.root/('pilot-complete.json' if a.pilot else 'complete.json')).read_text())
def check(run):
 dest=a.root/run['case']/run['method'];path=dest/'independent.json.gz'
 if (dest/'independent-receipt.json').exists():
  r=json.loads((dest/'independent-receipt.json').read_text());assert r['nativeSHA256']==run['outputSHA256'] and r['outputSHA256']==sha(path.read_bytes()) and r['scriptSHA256']==sha(script.read_bytes());return r
 if run['status']!='passed-native-stage':return dict(case=run['case'],method=run['method'],status='native-stage-unavailable')
 raw=(dest/'native.json.gz').read_bytes();assert sha(raw)==run['outputSHA256'];n=json.loads(gzip.decompress(raw))['fit']
 raw=(a.reference_root/run['case']/'reference.json.gz').read_bytes();assert sha(raw)==run['referenceSHA256'];ref=json.loads(gzip.decompress(raw))['results'][run['method']]
 request=dict(trendDispersions=ref['trendDispersions'],abundanceCovariates=ref['abundanceCovariates'],averageQuasiDispersion=n['averageQuasiDispersion'],initialMeans=[f['means'] for f in n['initialFits']],finalMeans=[f['means'] for f in n['refittedFits']],initialEffects=[f['effect'] for f in n['initialFits']],finalEffects=[f['effect'] for f in n['refittedFits']],updates=[dict(eligibleIndices=u['eligibleIndices'],quasiDispersions=u['quasiDispersions'],smoother=u['smoother']['fitted']) for u in n['updates']])
 data=gzip.compress(json.dumps(request,separators=(',',':')).encode(),mtime=0);(dest/'independent-input.json.gz').write_bytes(data)
 source=a.ql_root/run['case']/'input.json.gz';assert sha(source.read_bytes())==run['inputSHA256']
 with (dest/'independent.log').open('wb') as log:
  process=subprocess.run(['/opt/homebrew/bin/Rscript',str(script),str(source),str(dest/'independent-input.json.gz'),str(dest/'independent.json')],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 r=dict(case=run['case'],method=run['method'],exitCode=process.returncode,nativeSHA256=run['outputSHA256'],requestSHA256=sha(data),scriptSHA256=sha(script.read_bytes()))
 if (dest/'independent.json').exists():
  raw=(dest/'independent.json').read_bytes();result=json.loads(raw);path.write_bytes(gzip.compress(raw,mtime=0));assert gzip.decompress(path.read_bytes())==raw;(dest/'independent.json').unlink();r.update(outputSHA256=sha(path.read_bytes()),logicalSHA256=sha(raw),status=result['status'],maximumFitRelativeError=result['maximumFitRelativeError'],maximumSmootherRelativeError=result['maximumSmootherRelativeError'])
 else:r['status']='reference-process-failed'
 (dest/'independent-receipt.json').write_text(json.dumps(r,sort_keys=True,indent=2)+'\n');print(json.dumps(r),flush=True);return r
with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(check,runs))
(a.root/('pilot-independent.json' if a.pilot else 'independent-complete.json')).write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
