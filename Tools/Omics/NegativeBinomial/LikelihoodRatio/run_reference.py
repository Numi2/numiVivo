#!/usr/bin/env python3
"""Run pinned edgeR full/null sensitivity and retain every numerical discrepancy."""
import argparse,gzip,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True);p.add_argument('--tighter-full-fit',action='store_true');a=p.parse_args()
suffix='-tight' if a.tighter_full_fit else ''
protocol=json.loads((a.root/'protocol.json').read_text());results=[]
for rec in protocol['cases']:
 if rec['family']!='treatment':continue
 d=a.root/rec['id'];output=d/('r-output'+suffix+'.json');assert not output.exists()
 with (d/('reference'+suffix+'.log')).open('w') as log:
  run=subprocess.run(['/opt/homebrew/bin/Rscript',str(Path(__file__).with_name('reference.R')),str(d/'r-input.json.gz'),str(output)]+(['--tighter-full-fit'] if a.tighter_full_fit else []),stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':str(a.r_library)})
 result=dict(case=rec['id'],exitCode=run.returncode,referenceScriptSHA256=hashlib.sha256(Path(__file__).with_name('reference.R').read_bytes()).hexdigest())
 if run.returncode==0:
  raw=output.read_bytes();reference=json.loads(raw);native=json.loads(gzip.decompress((d/'native.json.gz').read_bytes()))
  fitted={g['feature']['featureIndex']:g['likelihoodRatio'] for g in native['genes'] if g['feature']['status']=='tested'}
  assert len(reference['results'])==len(fitted);errors=[];stat=0.0;prob=0.0
  for row in reference['results']:
   fit=fitted[row['featureIndex']];sd=abs(row['statistic']-fit['statistic']);pd=abs(row['pValue']-fit['pValue']);stat=max(stat,sd);prob=max(prob,pd)
   if row['failed'] or sd>2e-5 or pd>2e-6:errors.append(dict(featureIndex=row['featureIndex'],statisticDifference=sd,pValueDifference=pd,failed=row['failed']))
  packed=gzip.compress(raw,mtime=0);(d/('r-output'+suffix+'.json.gz')).write_bytes(packed);assert gzip.decompress(packed)==raw;output.unlink()
  result.update(genes=len(fitted),status='passed' if not errors else 'numerical-disagreements',maximumStatisticDifference=stat,maximumPValueDifference=prob,errors=errors,warnings=reference['warnings'],messages=reference['messages'],outputSHA256=hashlib.sha256(packed).hexdigest())
 (d/('reference'+suffix+'-check.json')).write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');results.append(result);print(json.dumps({k:v for k,v in result.items() if k!='errors'}),flush=True)
(a.root/('reference'+suffix+'-complete.json')).write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
