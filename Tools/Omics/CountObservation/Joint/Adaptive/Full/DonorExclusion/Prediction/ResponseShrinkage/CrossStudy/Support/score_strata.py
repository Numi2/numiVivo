from pathlib import Path
import json,gzip,hashlib,math
import numpy as np
s=Path(__file__).parent;source=Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs');freeze=json.loads((source/'input-freeze.json').read_text());diagnosis=json.loads((s/'control-only-diagnostic.json').read_text());panel=json.loads((source/'panel.json').read_text());folds=json.loads((source/'folds.json').read_text());logs={};out={'scope':'Post-hoc paired-error decomposition of all frozen predictions, by observed training-range membership. No causal or prospective applicability gate is claimed. All genes/donors retained.','controlDiagnosticSHA256':hashlib.sha256((s/'control-only-diagnostic.json').read_bytes()).hexdigest(),'folds':{}}
for name in ['panel.json','folds.json']:assert hashlib.sha256((source/name).read_bytes()).hexdigest()==freeze['files'][name]
for study in ['Kang','HIRISA']:
 p=source/(study+'-source-counts.npz');assert hashlib.sha256(p.read_bytes()).hexdigest()==freeze['files'][p.name]
 with np.load(p) as a:
  counts=a['counts'];ids=a['featureIDs'].tolist();idx=[ids.index(v) for v in panel];logs[study]=np.log1p(counts.astype('f8')/counts.sum(axis=1,dtype=np.uint64)[:,None]*1e6)[:,idx]
for fold in folds:
 tag=fold['id'];p=s/(tag+'-support.jsonl.gz');assert hashlib.sha256(p.read_bytes()).hexdigest()==diagnosis['folds'][tag]['supportSHA256'];rows=[json.loads(l) for l in gzip.open(p,'rt')];assert [r['featureID'] for r in rows]==panel
 truth=logs[fold['queryStudy']][fold['scoringIndex']];base=np.array([r['meanBaseline'] for r in rows]);prediction=np.array([r['frozenPrediction'] for r in rows]);outside=np.array([r['outsideObservedTrainingRange'] for r in rows]);delta=prediction-base;residual=base-truth;error=prediction-truth;decomposition=delta**2+2*delta*residual
 np.testing.assert_allclose(error**2-residual**2,decomposition,rtol=1e-9,atol=1e-12)
 result={'strata':{}}
 for name,mask in [('all',np.ones(len(rows),dtype=bool)),('inside',~outside),('outside',outside)]:
  n=int(mask.sum());result['strata'][name]={'features':n,'candidateRMSE':float(np.sqrt(np.mean(error[mask]**2))) if n else None,'meanBaselineRMSE':float(np.sqrt(np.mean(residual[mask]**2))) if n else None,'squaredPredictionShiftSum':float(np.sum(delta[mask]**2)),'twiceShiftTimesBaselineResidualSum':float(np.sum(2*delta[mask]*residual[mask])),'candidateMinusBaselineSSE':float(np.sum(error[mask]**2-residual[mask]**2))}
 out['folds'][tag]=result
out['status']='complete-all-26-folds';(s/'error-decomposition.json').write_text(json.dumps(out,indent=2)+'\n')
for k,v in out['folds'].items():
 if k.startswith('cross'):
  print(k,[(n,round(z['candidateRMSE']/z['meanBaselineRMSE'],3)) for n,z in v['strata'].items()],flush=True)
