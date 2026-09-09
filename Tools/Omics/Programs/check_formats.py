#!/usr/bin/env python3
"""Controlled numerical fixture: signed scores, missing features and sparse orientation."""
import argparse,json,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
x=np.array([[2,1,7],[0,0,5],[0,0,0]],dtype=np.int32)
d=dict(id='signed',organism='fixture',featureNamespace='fixture',sourceURI='urn:numivivo:numerical-fixture',sourceVersion='1',sourceDescription='Numerical fixture only',
       members=[dict(featureID='A',weight=2),dict(featureID='B',weight=-1),dict(featureID='missing',weight=1)],minimumWeightCoverage=0.75)
mapping=dict(schemaVersion=1,id='program-fixture',evidence='synthetic',countUnit='umiCount',matrixPath='X',sourceDescription='Numerical fixture only',sampleColumn='sample',
    samples=[dict(id='s',biologicalReplicateID='d',condition='c',batchID='unreported',organism='fixture')])
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[],programs=dict(definitions=[d]));(a.out/'plan.json').write_text(json.dumps(plan))
results=[]
for encoding in ['csr','csc']:
    obj=ad.AnnData(getattr(sparse,encoding+'_matrix')(x),obs=pd.DataFrame({'sample':['s']*3},index=['c0','c1','c2']),var=pd.DataFrame(index=['A','B','C']))
    path=a.out/(encoding+'.h5ad');obj.write_h5ad(path,compression='gzip')
    r=subprocess.run([str(a.binary.resolve()),'singlecell-h5ad-pseudobulk',str(path),'--plan',str(a.out/'plan.json'),'--output',str(a.out/encoding)],capture_output=True,text=True)
    (a.out/(encoding+'.stderr')).write_text(r.stderr);assert r.returncode==0,r.stderr
    result=json.loads((a.out/encoding/'report.json').read_text())['programs'];results.append(result)
assert results[0]==results[1]
r=results[0];expected=(2*np.log1p(2000)-np.log1p(1000))/3
assert abs(r['scores'][0][0]-expected)<1e-14
assert r['scores'][1:]==[[0],[None]] and r['detectedMembers']==[[2],[0],[0]]
assert r['programs'][0]['missingFeatureIDs']==['missing'] and r['programs'][0]['weightCoverage']==0.75
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',checks=['exact CSR/CSC program result','signed L1 score with all-gene denominator','zero-expression score distinct from missing empty-library score','explicit missing-feature coverage'],qualification='Numerical fixture only'),indent=2)+'\n')
print('Program orientation and missing/zero handling passed')
