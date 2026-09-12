from pathlib import Path
import json,hashlib,math
import numpy as np
r=Path('/Users/n/numivivo-context-weight-audit-20260912');p=Path('/Users/n/numivivo-parse-context-evaluation-20260912')
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
assert h(p/'inputs/Parse.json')=='f68d55774b70e62d931876552756543467a68b7566271245df894d204f694169'
assert h(p/'native/Parse.json')=='d9fe231a74e301f82f6793544f588cd08796f24664d818218a7a33bc4f6883a6'
a=json.loads((p/'inputs/Parse.json').read_text());b=json.loads((p/'native/Parse.json').read_text());studies=np.array(a['trainingStudies']);w=np.array(b['model']['donorWeights']);beta=np.array(b['model']['queryResponseWeights']);x=np.array(a['controls']);q=np.array(a['queryControls'])
summary={s:{'donors':int(np.sum(studies==s)),'baseWeight':float(w[studies==s].sum())} for s in sorted(set(studies))}
assert all(abs(v['baseWeight']-1/len(summary))<1e-12 for v in summary.values())
rows=[]
for i,donor in enumerate(a['queryDonorIDs']):
 distances=np.sqrt(np.mean((x-q[i])**2,axis=1));nearest=int(np.argmin(distances));row={'donor':donor,'weightSum':float(beta[i].sum()),'minimumDonorWeight':float(beta[i].min()),'maximumDonorWeight':float(beta[i].max()),'negativeWeightMass':float(-beta[i][beta[i]<0].sum()),'totalAbsoluteWeight':float(np.abs(beta[i]).sum()),'responseStudyWeights':{s:float(beta[i][studies==s].sum()) for s in summary},'nearestControlDonor':a['trainingDonorIDs'][nearest],'nearestControlRMSDistance':float(distances[nearest])};assert abs(math.fsum(beta[i])-1)<1e-12;rows.append(row)
out={'scope':'Retrospective frozen-model weight audit; no targets read, fitting, prediction changes or new validation','sourceHashes':{str(p/f):h(p/f) for f in ['inputs/Parse.json','native/Parse.json']},'studies':summary,'queries':rows};(r/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
