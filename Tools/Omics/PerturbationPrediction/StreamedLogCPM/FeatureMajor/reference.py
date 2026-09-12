from pathlib import Path
import json,hashlib,numpy as np
p=Path('/Users/n/numivivo-logcpm-feature-major-20260912');source=Path('/Users/n/numivivo-count-store-final-20260909/store');sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((p/'protocol.json').read_text());assert sha(p/'plan.json')==protocol['planSHA256'];plan=json.loads((p/'plan.json').read_text());receipt=json.loads((source/'receipt.json').read_text());groups=np.array(plan['rowGroups']);totals=np.array(plan['rowTotals'],dtype=np.uint64);g=len(plan['featureIDs']);k=len(plan['groupIDs']);sums=np.zeros(k*g);observed=np.zeros(len(totals),dtype=np.uint64);dtype=np.dtype([('row','<u4'),('feature','<u4'),('count','<u8')]);digest=hashlib.sha256();entries=0;previous=-1
with (source/'counts.bin').open('rb') as f:
 while True:
  b=f.read(1048576)
  if not b:break
  digest.update(b);a=np.frombuffer(b,dtype=dtype);row=a['row'];col=a['feature'];key=col.astype(np.uint64)*len(totals)+row;assert key[0]>previous and np.all(key[1:]>key[:-1]);previous=int(key[-1]);np.add.at(observed,row,a['count']);values=np.log1p(a['count'].astype(float)/totals[row]*1e6)
  # Indexed accumulation avoids a fresh group-by-feature allocation per chunk.
  np.add.at(sums,groups[row]*g+col,values);entries+=len(a)
assert digest.hexdigest()==bytes(receipt['counts']['bytes']).hex() and entries==receipt['entries'];assert np.array_equal(observed,totals)
sizes=np.bincount(groups,minlength=k);reference=sums.reshape(k,g)/sizes[:,None];native=json.loads((p/'native.json').read_text());assert native['featureIDs']==plan['featureIDs'] and native['groupIDs']==plan['groupIDs'];assert native['cellCounts']==sizes.tolist();assert native['zeroCellCounts']==np.bincount(groups[totals==0],minlength=k).tolist();actual=np.array(native['means']);assert actual.shape==reference.shape;error=float(np.max(abs(actual-reference)));assert error<1e-10
np.save(p/'reference.npy',reference,allow_pickle=False)
r=dict(status='PASS-all-condition-means',cells=len(totals),features=g,groups=k,entries=entries,values=actual.size,maximumAbsoluteError=error,countsSHA256=digest.hexdigest(),nativeSHA256=sha(p/'native.json'),referenceSHA256=sha(p/'reference.npy'),scriptSHA256=sha(Path(__file__)),scope='Numerical full-matrix normalization; no biological prediction or independent replication claim.')
(p/'verification.json').write_text(json.dumps(r,indent=2));print(r)
