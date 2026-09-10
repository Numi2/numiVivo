#!/usr/bin/env python3
"""Freeze training-only count inputs and separately sealed held-out rows."""
import argparse,copy,gzip,hashlib,itertools,json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','source_report','product_report']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.root.mkdir(exist_ok=True)
def load(p):return json.loads(gzip.decompress(p.read_bytes()))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):
 raw=json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
 p.write_bytes(gzip.compress(raw,mtime=0) if p.suffix=='.gz' else raw+b'\n')
report=load(a.source_report);product=load(a.product_report)
assert report['pseudobulk']==product['pseudobulk'] and report['metadata']==product['metadata']
metadata=report['metadata'];bulk=report['pseudobulk'];m=bulk['matrix'];groups=bulk['groups']
matrix=csr_matrix((np.array(m['counts'],dtype=np.int64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount']))
assert matrix.shape==(63,11076) and int(matrix.sum())==67349705
write(a.root/'metadata.json.gz',metadata)
write(a.root/'probe.json.gz',dict(bulk=bulk,request=product['contrasts'][0]['request']))
def selected(rows):
 v=matrix[rows].tocsr();return dict(bulk,groups=[groups[i] for i in rows],matrix=dict(cellCount=len(rows),featureCount=v.shape[1],rowOffsets=v.indptr.tolist(),featureIndices=v.indices.tolist(),counts=v.data.tolist()))
records=[];unavailable=[]
for group in sorted({g['cellGroup'] for g in groups}):
 rows=[i for i,g in enumerate(groups) if g['cellGroup']==group and len(g['sourceCellIndices'])>=10]
 arms={c:[i for i in rows if groups[i]['condition']==c] for c in ['Vehicle','LPS']}
 if min(map(len,arms.values()))<4:
  unavailable.append(dict(cellGroup=group,reason='Fewer than four independent animals per arm; cannot retain three training animals',replicates={k:len(v) for k,v in arms.items()}));continue
 assert all(len(v)==4 for v in arms.values())
 for control,treated in itertools.product(arms['Vehicle'],arms['LPS']):
  test_rows=[control,treated];train_rows=[i for i in rows if i not in test_rows]
  name=f'{len(records):03d}';d=a.root/name;d.mkdir(exist_ok=False)
  request=copy.deepcopy(product['contrasts'][0]['request']);request.update(id='crowell-heldout-'+name,cellGroup=group,sizeFactors='librarySize',includedDonorIDs=sorted(groups[i]['donorID'] for i in train_rows))
  request['negativeBinomialOptions']=dict(trend='gammaParametric',effectPriorEstimation='weightedUpperQuantile')
  training=selected(train_rows)
  assert set(g['donorID'] for g in training['groups']).isdisjoint(groups[i]['donorID'] for i in test_rows)
  write(d/'training.json.gz',dict(bulk=training,request=request))
  write(d/'test.json.gz',dict(sampleIDs=[groups[i]['sampleIDs'][0] for i in test_rows],conditions=[groups[i]['condition'] for i in test_rows],counts=matrix[test_rows].toarray().tolist()))
  records.append(dict(id=name,cellGroup=group,trainingSourceRows=train_rows,testSourceRows=test_rows,trainingAnimals=request['includedDonorIDs'],testAnimals=[groups[i]['donorID'] for i in test_rows],trainingSHA256=sha(d/'training.json.gz'),testSHA256=sha(d/'test.json.gz')))
assert len(records)==112
write(a.root/'protocol.json',dict(sourceReport=str(a.source_report),sourceReportSHA256=sha(a.source_report),productReport=str(a.product_report),productReportSHA256=sha(a.product_report),metadataSHA256=sha(a.root/'metadata.json.gz'),sourceCommit='e5ec73f80c69114553fbe249ee400aeab6dec7db',protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),cases=records,unavailable=unavailable,qualification='All 112 predeclared overlapping folds; eight independent animals; conditional on observed library depth; previously inspected source; no FDR or interval coverage claim'))
print(json.dumps(dict(cases=len(records),unavailable=unavailable)))
