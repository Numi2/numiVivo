#!/usr/bin/env python3
"""Reproduce all reference fits and compare retained native dispersion stages."""
import argparse,gzip,hashlib,json,os,shlex,subprocess,time
from pathlib import Path
import numpy as np
import pandas as pd

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--r-library',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
owner=Path(__file__).resolve().parent
def sha(b):return hashlib.sha256(b).hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
write(a.out/'declaration.json',dict(protocolSHA256=sha((owner/'PROTOCOL.md').read_bytes()),
    runnerSHA256=sha(Path(__file__).read_bytes()),referenceSHA256=sha((owner/'stages.R').read_bytes()),
    originalResultsSHA256=sha((a.root/'final-results/results.json').read_bytes()),postScoreDiagnostic=True))
hits=json.loads((a.root/'final-results/native-default-null-hits.json').read_text())
# The previously published diagnostics provide the exact declared call identities.
hit_ids={(h['study'],h['seed'],h['featureID']) for h in json.loads((a.root/'final-results/null-hit-diagnostics.json').read_text())['hits']}
records=[];call_records=[]
for study in ['kang','hagai']:
 for seed in range(1,11):
  d=a.root/study/str(seed);out=a.out/study/str(seed);out.parent.mkdir(parents=True,exist_ok=True)
  started=time.monotonic()
  with (out.parent/(str(seed)+'.log')).open('w') as log:
   process=subprocess.run(['/opt/homebrew/bin/Rscript',str(owner/'stages.R'),str(d/'r-input'),str(out)],
    stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'R_LIBS_USER':str(a.r_library),'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
  record=dict(study=study,seed=seed,exitCode=process.returncode,seconds=time.monotonic()-started)
  records.append(record);write(a.out/'runs.json',records)
  if process.returncode:continue
  stage=json.loads((out/'stages.json').read_text())
  r=pd.read_csv(out/'stages.tsv.gz',sep='\t',index_col='featureID',float_precision='round_trip')
  old=pd.read_csv(d/'r-reference/native_size_factors-DESeq2.tsv.gz',sep='\t',index_col='featureID',float_precision='round_trip')
  assert list(r.index)==list(old.index)
  errors={}
  for newcol,oldcol in [('finalDispersion','dispersion'),('log2FoldChange','log2FoldChange'),('standardError','standardError'),('pValue','pValue'),('adjustedPValue','adjustedPValue')]:
   v=r[newcol].to_numpy(float);w=old[oldcol].to_numpy(float)
   assert np.allclose(v,w,rtol=1e-9,atol=1e-12,equal_nan=True),(study,seed,newcol)
   errors[newcol]=float(np.nanmax(np.abs(v-w)))
  assert np.array_equal(r.dispersionOutlier,old.dispersionOutlier)
  assert np.array_equal(r.betaConverged,old.betaConverged)
  ar=json.loads((d/'default/archive.json').read_text())
  report=next(e for e in ar['entries'] if e['logicalPath']=='report.json')
  remote='/Users/n/numivivo-null-benchmark-20260910/'+study+'/'+str(seed)+'/default/report.json.gz'
  # Hash the full model on its owning host; return only small stage diagnostics.
  code="""import gzip,hashlib,json,sys
from pathlib import Path
b=Path(sys.argv[1]).read_bytes();assert hashlib.sha256(b).hexdigest()==sys.argv[2]
raw=gzip.decompress(b);assert hashlib.sha256(raw).hexdigest()==sys.argv[3]
c=json.loads(raw)['contrasts'][0];diagnostics={v['featureIndex']:v for v in c['negativeBinomial']['features']}
rows=[]
for f in c['features']:
 v=diagnostics[f['featureIndex']]
 rows.append(dict(featureID=f['featureID'],status=f['status'],meanNormalizedCount=f['meanNormalizedCount'],**{k:v.get(k) for k in ['geneWiseDispersion','geneWiseLowerBoundary','geneWiseUpperBoundary','trendDispersion','dispersionOutlier','finalDispersion']}))
sys.stdout.buffer.write(gzip.compress(json.dumps(dict(trend=c['negativeBinomial']['trend'],residualDF=c['design']['residualDegreesOfFreedom'],features=rows)).encode(),mtime=0))
"""
  cmd='python3 -c '+shlex.quote(code)+' '+shlex.join([remote,report['sha256'],report['logicalSHA256']])
  compressed=subprocess.check_output(['ssh','macmini',cmd],timeout=60)
  native=json.loads(gzip.decompress(compressed));(out/'native-stages.json.gz').write_bytes(compressed)
  n=pd.DataFrame(native['features']).set_index('featureID');joint=n.loc[r.index]
  assert np.allclose(joint.meanNormalizedCount,r.meanNormalizedCount,rtol=1e-12,atol=1e-10)
  def ratios(col):
   valid=joint[col].notna() & (joint[col]>0) & r[col].notna()
   v=(r.loc[valid,col]/joint.loc[valid,col]).to_numpy(float)
   return dict(genes=len(v),minimum=float(v.min()),median=float(np.median(v)),maximum=float(v.max()))
  record.update(status='reproduced-and-audited',reference=stage,
    nativeTrend={k:v for k,v in native['trend'].items() if not isinstance(v,list)},
    nativeReferenceGenes=len(native['trend']['referenceFeatureIndices']),
    nativeResidualDF=native['residualDF'],referenceReproductionErrors=errors,
    originalReportSHA256=report['logicalSHA256'],nativeStagesSHA256=sha(compressed),
    dispersionRatios={col:ratios(col) for col in ['geneWiseDispersion','trendDispersion','finalDispersion']})
  for gene in r.index:
   if (study,seed,gene) not in hit_ids:continue
   row=dict(study=study,seed=seed,featureID=gene)
   for col in ['geneWiseDispersion','trendDispersion','finalDispersion']:
    row[col]=dict(native=float(joint.loc[gene,col]),reference=float(r.loc[gene,col]))
   row['nativeOutlier']=bool(joint.loc[gene,'dispersionOutlier']);row['referenceOutlier']=bool(r.loc[gene,'dispersionOutlier'])
   row['referenceAdjustedPValue']=float(r.loc[gene,'adjustedPValue']);call_records.append(row)
  write(a.out/'runs.json',records);write(a.out/'null-calls.json',call_records)
  print(json.dumps(dict(study=study,seed=seed,status=record['status'],nativePrior=native['trend']['priorLogVariance'],referencePrior=stage['priorLogVariance'])),flush=True)
assert len(records)==20 and all(r.get('status')=='reproduced-and-audited' for r in records)
assert len(call_records)==20
write(a.out/'complete.json',dict(status='completed-all-twenty-stage-audits',runs=20,calls=20,
    qualification='Post-score stage diagnosis; no native inference changes or new calibration/power qualification.'))
