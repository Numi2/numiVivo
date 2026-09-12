from pathlib import Path
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import json,hashlib,itertools,numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score
r=Path(__file__).parent;p=json.loads((r/'protocol.json').read_text());terminal=json.loads((r/'terminal.json').read_text());assert len(terminal['results'])==12 and all(x['exitCode']==0 for x in terminal['results'])
parent=Path('/Users/n/numivivo-native-metal-knn-20260912');expected=json.loads((parent/'pca/metadata.json').read_text())['cells']
record=np.dtype([('row','<u4'),('col','<u4'),('weight','<f8')]);graphs={mode:np.fromfile(parent/mode/'edges.bin',record) for mode in ['cpu','metal']}
labels={};checks=[]
for mode in graphs:
 edges=graphs[mode];degrees=np.bincount(edges['row'],weights=edges['weight'],minlength=p['cells']);total=float(edges['weight'].sum())
 for seed in p['seeds']:
  result=json.loads((r/(mode+'-'+str(seed))/'result.json').read_text());a=np.asarray(result['labels']);assert len(a)==p['cells'] and np.all(a>=0)
  assert result['cells']==[{k:x[k] for k in ['sampleID','barcode']} for x in expected]
  assert result['options']['seed']==seed and result['options']['resolution']==1
  volumes=np.bincount(a,weights=degrees);internal=edges['weight'][a[edges['row']]==a[edges['col']]].sum();q=float(internal/total-np.sum((volumes/total)**2));assert abs(q-result['modularity'])<1e-10
  assert sorted(np.bincount(a).tolist())==sorted(result['clusterSizes'])
  labels[mode,seed]=a;checks.append(dict(mode=mode,seed=seed,clusters=len(np.unique(a)),modularity=result['modularity'],modularityOracleError=abs(q-result['modularity']),termination=result['termination'],disconnectedCommunities=result['disconnectedCommunities']))
comparisons=[]
for seed in p['seeds']:
 a=labels['cpu',seed];b=labels['metal',seed];table=np.zeros((int(a.max())+1,int(b.max())+1),dtype=np.int64);np.add.at(table,(a,b),1);rows,cols=linear_sum_assignment(table,maximize=True);agreement=float(table[rows,cols].sum()/len(a));ari=float(adjusted_rand_score(a,b))
 comparisons.append(dict(seed=seed,ARI=ari,matchedAssignmentAgreement=agreement,matchedChangedCells=int(len(a)-table[rows,cols].sum()),contingency=table.tolist(),passes=ari>=p['acceptance']['minimumARI'] and agreement>=p['acceptance']['minimumMatchedAssignmentAgreement']))
within=[dict(mode=mode,seeds=[a,b],ARI=float(adjusted_rand_score(labels[mode,a],labels[mode,b]))) for mode in graphs for a,b in itertools.combinations(p['seeds'],2)]
result=dict(status='PASS' if all(x['passes'] for x in comparisons) else 'FAIL',scope=p['scope'],protocolSHA256=hashlib.sha256((r/'protocol.json').read_bytes()).hexdigest(),cells=p['cells'],comparisons=comparisons,withinBackendSeedComparisons=within,nativeChecks=checks,biologicalAccuracyQualified=False)
(r/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['comparisons','nativeChecks']}));print(json.dumps([{k:v for k,v in x.items() if k!='contingency'} for x in comparisons]))
