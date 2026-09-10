#!/usr/bin/env python3
"""Native selected aggregation, source binding and excluded-row isolation checks."""
import argparse,copy,hashlib,json,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);checks=[]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run(name,source,plan,error=None):
    path=a.out/(name+'.json');path.write_text(json.dumps(plan));dest=a.out/name
    r=subprocess.run([str(a.binary),'singlecell-h5ad-pseudobulk',str(source),'--plan',str(path),'--output',str(dest)],capture_output=True,text=True)
    (a.out/(name+'.log')).write_text(r.stdout+r.stderr)
    assert r.returncode==(65 if error else 0),(name,r.returncode,r.stderr)
    if error:assert error.lower() in r.stderr.lower() and not dest.exists(),(name,r.stderr)
    else:assert sha(dest/'original.h5ad')==sha(source)
    checks.append(name)
    return json.loads((dest/'report.json').read_text()) if not error else None
samples=[dict(id=s,condition=s,biologicalReplicateID=s,batchID='b',organism='NCBITaxon:9606') for s in ['s0','s1']]
mapping=dict(schemaVersion=1,id='selection-fixture',evidence='synthetic',sourceDescription='Cell selection interoperability only',countUnit='umiCount',matrixPath='X',samples=samples,sampleColumn='sample',mitochondrialFeatureIDs=['g0'])
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[])
counts=np.array([[2,0,3],[0,71,0],[5,0,7],[0,0,0]],dtype=np.int64)
obs=pd.DataFrame({'sample':['s0','s1','s0','s1']},index=['c0','c1','c2','c3'])
var=pd.DataFrame(index=['g0','g1','g2'])
def make(name,values,layout):
    source=a.out/(name+'.h5ad')
    matrix=values if layout=='dense' else sparse.csr_matrix(values).asformat(layout)
    ad.AnnData(X=matrix,obs=obs,var=var).write_h5ad(source);return source

def selected(source,indices=[2,0,3]):
    return dict(plan,cellSelection=dict(source=dict(bytes=list(bytes.fromhex(sha(source)))),observationIndices=indices,provenance='Fixture explicit cohort'))
for layout in ['csr','csc','dense']:
    source=make(layout,counts,layout);sp=selected(source)
    r=run(layout+'-selected',source,sp)
    assert r['sourceCellCount']==4 and r['sourceObservationIndices']==[2,0,3]
    assert [c['barcode'] for c in r['metadata']['cells']]==['c2','c0','c3']
    assert [v['totalCounts'] for v in r['quality']]==[12,5,0]
    assert [v['mitochondrialCounts'] for v in r['quality']]==[5,2,0]
    bulk=r['pseudobulk'];m=bulk['matrix']
    actual=sparse.csr_matrix((m['counts'],m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray()
    assert np.array_equal(actual,np.array([[7,0,10],[0,0,0]]))
    assert [g['sourceCellIndices'] for g in bulk['groups']]==[[0,1],[2]]
    assert r['canonicalNonzeros']==4
    v=subprocess.run([str(a.binary),'singlecell-h5ad-pseudobulk-verify',str(a.out/(layout+'-selected'))],capture_output=True,text=True)
    (a.out/(layout+'-verify.log')).write_text(v.stdout+v.stderr);assert v.returncode==0,v.stderr;checks.append(layout+'-verify')
    original=run(layout+'-all',source,plan)
    assert 'sourceObservationIndices' not in original and 'sourceCellCount' not in original
    assert [v['totalCounts'] for v in original['quality']]==[5,71,12,0]
    mutated=counts.copy();mutated[1]=[100001,222222,333333]
    changed=make(layout+'-excluded-mutated',mutated,layout)
    rr=run(layout+'-excluded-mutated',changed,selected(changed))
    assert rr==r,'Excluded counts affected selected report'
    invalid=counts.astype(float);invalid[1,1]=-1
    bad=make(layout+'-invalid-excluded',invalid,layout)
    run(layout+'-reject-invalid-excluded',bad,selected(bad),'count')
    run(layout+'-reject-stale-source',changed,sp,'selection source changed')
source=a.out/'csr.h5ad'
for name,indices,message in [('duplicate',[0,0],'indices or provenance'),('negative',[-1],'indices or provenance'),('empty',[],'indices or provenance'),('outside',[4],'outside source axis')]:
    run('reject-'+name,source,selected(source,indices),message)
sp=selected(source);sp['cellSelection']['unexpected']=True
run('reject-unknown-selection-key',source,sp,'unknown')
sp=selected(source);sp['cellSelection']['provenance']=' '
run('reject-empty-provenance',source,sp,'indices or provenance')
sp=selected(source);sp['reduction']={}
run('reject-selection-with-reduction',source,sp,'aggregation and contrasts only')
# Replay cannot accept a reordered source-row mapping with the old receipt.
path=a.out/'csr-selected/plan.json';x=json.loads(path.read_text());x['cellSelection']['observationIndices']=[0,2,3];path.write_text(json.dumps(x))
v=subprocess.run([str(a.binary),'singlecell-h5ad-pseudobulk-verify',str(path.parent)],capture_output=True,text=True)
(a.out/'reject-tampered-selection.log').write_text(v.stdout+v.stderr)
assert v.returncode==65 and 'changed' in v.stderr;checks.append('reject-tampered-selection')
assert not list(a.out.glob('.numivivo-stream-*'))
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',binarySHA256=sha(a.binary),checks=checks),indent=2)+'\n')
print(json.dumps(dict(status='passed',checks=len(checks))))
