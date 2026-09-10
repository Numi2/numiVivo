#!/usr/bin/env python3
"""Independently check every reference QL stage and summarize frozen families."""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import f as f_distribution
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
sha=lambda b:hashlib.sha256(b).hexdigest()
def bh(p):
 o=np.argsort(p,kind='stable');q=np.empty(len(p));q[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
protocol=json.loads((root/'protocol.json').read_text());complete=json.loads((root/'complete.json').read_text());assert len(complete)==len(protocol['cases'])==29
summaries=[]
for c in protocol['cases']:
 d=root/c['id'];receipt=json.loads((d/'input-receipt.json').read_text());raw=(d/'input.json.gz').read_bytes();assert sha(raw)==receipt['inputSHA256'];assert sha(gzip.decompress(raw))==receipt['logicalSHA256'];i=json.loads(gzip.decompress(raw));run=json.loads((d/'run.json').read_text());assert run['exitCode']==0
 raw=(d/'reference.json.gz').read_bytes();assert sha(raw)==run['outputSHA256'];assert sha(gzip.decompress(raw))==run['logicalSHA256'];r=json.loads(gzip.decompress(raw));assert set(r['results'])=={o+'-'+v for o in ['nativeTrend','edgeRTrend'] for v in ['adjusted','legacy']}
 details=[]
 for name,result in r['results'].items():
  rows=result['table'];assert [v['featureIndex'] for v in rows]==i['featureIndices']
  col=lambda name:np.asarray([x[name] for x in rows],dtype=float)
  residual=col('residualDF');prior=col('priorDF');scale=col('posteriorScale');deviance=col('residualDeviance');priorScale=col('priorScale');lr=col('LR')
  statistic=lr/scale;df=np.minimum(prior+residual,np.sum(result['ordinaryResidualDF']));probability=np.maximum(f_distribution.sf(statistic,1,df),col('poissonBound'));q=bh(probability)
  finite=np.isfinite(prior);posterior=priorScale.copy();posterior[finite]=(deviance[finite]+prior[finite]*priorScale[finite])/(residual[finite]+prior[finite])
  rel=lambda actual,expected:np.abs(actual-expected)/np.maximum(1,np.abs(expected))
  metrics=dict(F=rel(col('F'),statistic),DF=rel(col('denominatorDF'),df),pValue=np.abs(col('pValue')-probability),BH=np.abs(col('BH')-q),posterior=rel(scale,posterior))
  if name.endswith('adjusted'):
   metrics['unitDF']=rel(np.sum(result['unitDF'],axis=1),residual);metrics['unitDeviance']=rel(np.sum(result['unitDeviance'],axis=1),deviance)
  limits=dict(F=1e-10,DF=1e-10,pValue=2e-10,BH=2e-10,posterior=1e-10,unitDF=1e-9,unitDeviance=1e-9)
  bad=np.zeros(len(rows),dtype=bool)
  for k,v in metrics.items():bad|=~np.isfinite(v)|(v>limits[k])
  maxima={k:float(np.max(v)) for k,v in metrics.items()}
  failed=[dict(featureIndex=rows[j]['featureIndex'],metrics={k:float(v[j]) for k,v in metrics.items()}) for j in np.flatnonzero(bad)]
  stage=dict(case=c['id'],family=c['family'],method=name,genes=len(rows),calls=int(np.sum(col('BH')<.05)),methodConstraint=result['methodConstraint'],arithmeticStatus='passed' if not failed else 'failed',failures=failed,maximumErrors=maxima,maximumScaledScore=float(np.max(col('scaledScore'))),optimizerFailureFlags=int(sum(v['failed'] for v in rows)),residualDFQuantiles=np.quantile(residual,[0,.25,.5,.75,1]).tolist(),posteriorScaleQuantiles=np.quantile(scale,[0,.25,.5,.75,1]).tolist(),stages=r['stages'][name])
  details.append(stage);summaries.append(stage)
 (d/'check.json').write_text(json.dumps(dict(case=c['id'],checkerSHA256=sha(Path(__file__).read_bytes()),inputSHA256=receipt['inputSHA256'],outputSHA256=run['outputSHA256'],stages=details),sort_keys=True,indent=2,allow_nan=False)+'\n')
null=[]
for study in ['kang','hagai']:
 for method in ['nativeTrend-adjusted','nativeTrend-legacy','edgeRTrend-adjusted','edgeRTrend-legacy']:
  cases=[s for s in summaries if s['family']=='null' and s['case'].startswith(study) and s['method']==method];assert len(cases)==10
  null.append(dict(study=study,method=method,calls=sum(s['calls'] for s in cases),splitsWithCalls=sum(s['calls']>0 for s in cases),arithmeticFailures=sum(s['arithmeticStatus']!='passed' for s in cases),methodConstraintFailures=sum(s['methodConstraint']!='passed' for s in cases)))
result=dict(cases=29,stages=len(summaries),geneAnalyses=sum(s['genes'] for s in summaries),status='completed' if all(s['arithmeticStatus']=='passed' and s['methodConstraint']=='passed' for s in summaries) else 'completed-with-retained-failures',arithmeticFailures=[s for s in summaries if s['arithmeticStatus']!='passed'],methodConstraintFailures=[dict(case=s['case'],method=s['method']) for s in summaries if s['methodConstraint']!='passed'],null=null,treatment=[{k:s[k] for k in ['case','method','genes','calls','arithmeticStatus','methodConstraint']} for s in summaries if s['family']=='treatment'],maximumErrors={k:max(s['maximumErrors'].get(k,0) for s in summaries) for k in ['F','DF','pValue','BH','posterior','unitDF','unitDeviance']},maximumScaledScore=max(s['maximumScaledScore'] for s in summaries),optimizerFailureFlags=sum(s['optimizerFailureFlags'] for s in summaries),qualification='Previously inspected full-support families; reference QL stage measurements, not a native QL method, FDR calibration, effect truth or power')
(root/'summary.json').write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['treatment','arithmeticFailures']},indent=2))
