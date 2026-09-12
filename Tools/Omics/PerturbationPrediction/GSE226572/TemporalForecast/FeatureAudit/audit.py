from pathlib import Path
import json,hashlib,gzip,csv,io
import numpy as np
r=Path(__file__).parent;p=Path('/Users/n/numivivo-duration-forecast-20260912');source=Path('/Users/n/numivivo-duration-source-20260912.npz');h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();assert h(source)=='534a0f5bcd676af37cda6af9b7085da012dc68edc90ac99c3fdffc06b971e341'
a=np.load(source);freeze=json.loads((p/'freeze.json').read_text());g=len(a['featureIDs']);errors={k:np.zeros(g) for k in ['candidate','trainingMean','noChange']};donor_counts={d:sum(x['donor']==d for x in freeze['records']) for d in {x['donor'] for x in freeze['records']}}
for x in freeze['records']:
 f=p/'outputs'/f"{x['fold']}.json";assert h(f)==x['outputSHA256'];b=json.loads(f.read_text());prefix='query:GSE226572:'+x['donor'];truth=a[prefix+'|truth|IFNB:'+str(int(x['time']))+'h'];control=a[prefix+'|noChange'];w=1/(len(donor_counts)*donor_counts[x['donor']])
 for k,v in [('candidate',b['predicted']),('trainingMean',b['baseline']),('noChange',control)]:errors[k]+=w*(np.asarray(v)-truth)**2
summary={}
for baseline in ['trainingMean','noChange']:
 delta=errors['candidate']-errors[baseline];order=np.argsort(delta);summary[baseline]={'lowerMSEGenes':int(np.sum(delta < -1e-12)),'higherMSEGenes':int(np.sum(delta > 1e-12)),'tiedWithin1e12':int(np.sum(np.abs(delta)<=1e-12)),'meanGeneMSECandidate':float(errors['candidate'].mean()),'meanGeneMSEBaseline':float(errors[baseline].mean()),'largestRegressions':[{'featureID':str(a['featureIDs'][i]),'excessMSE':float(delta[i])} for i in order[-10:][::-1]]}
buf=io.StringIO();writer=csv.writer(buf);writer.writerow(['featureID','candidateMSE','trainingMeanMSE','noChangeMSE'])
for i,name in enumerate(a['featureIDs']):writer.writerow([name]+[repr(float(errors[k][i])) for k in errors])
(r/'all-features.csv.gz').write_bytes(gzip.compress(buf.getvalue().encode(),mtime=0));out={'scope':'Retrospective all-feature error audit, no refit or feature selection; equal donor then equal within-donor time weighting','features':g,'absoluteMSETieTolerance':1e-12,'sourceSHA256':h(source),'predictionFreezeSHA256':h(p/'freeze.json'),'tableSHA256':h(r/'all-features.csv.gz'),'comparisons':summary};(r/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:{n:v[n] for n in ['lowerMSEGenes','higherMSEGenes','tiedWithin1e12']} for k,v in summary.items()}))
