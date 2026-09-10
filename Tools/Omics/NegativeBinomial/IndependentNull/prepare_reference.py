#!/usr/bin/env python3
"""Prepare exact aggregate counts and independently reconstructed native designs/offsets."""
import argparse, csv, gzip, hashlib, json
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
root=p.parse_args().root;sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
audit=json.loads((root/'count-audit.json').read_text());check=json.loads((root/'source-check.json').read_text())
assert check['status']=='passed' and sha(root/'reference-counts.npz')==audit['referenceCountsSHA256']
assert sha(root/'source-report.json.gz')==check['compressedReportSHA256']
# Only this previously hash-verified, locally generated audit archive uses object string axes.
with np.load(root/'reference-counts.npz',allow_pickle=True) as archive:
    counts=archive['counts'];donors=archive['donors'].tolist();features=archive['features'].tolist()
report=json.loads(gzip.decompress((root/'source-report.json.gz').read_bytes()));groups=report['pseudobulk']['groups']
manifest=json.loads((root/'requests/manifest.json').read_text());out=root/'reference-inputs';out.mkdir(exist_ok=False)
def table(path,header,rows):
    with path.open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)
for case in manifest:
    if case['method']!='wald':continue
    request=json.loads((root/'requests'/case['path']).read_text());selected=set(request['includedDonorIDs'])
    indices=[i for i,g in enumerate(groups) if g['donorID'] in selected and len(g['sourceCellIndices'])>=request['minimumCellsPerPseudobulk']]
    observations=[groups[i] for i in indices];ids=[g['donorID'] for g in observations]
    assert len(ids)==12 and set(ids)==selected
    y=counts[[donors.index(d) for d in ids]].T
    batches=sorted({g['batchIDs'][0] for g in observations})
    names=['intercept','treatment-minus-control']+['batch:'+b for b in batches[1:]]
    x=np.array([[1,int(g['condition']=='shamB')]+[int(g['batchIDs']==[b]) for b in batches[1:]] for g in observations])
    assert np.linalg.matrix_rank(x)==x.shape[1] and sum(x[:,1])==6
    reference=np.flatnonzero(np.all(y>0,axis=1));assert len(reference)>=10
    logged=np.log(y[reference].astype(float));ratios=np.exp(logged-logged.mean(axis=1,keepdims=True))
    logf=np.log(np.median(ratios,axis=0));factors=np.exp(logf-logf.mean())
    libraries=y.sum(axis=0);eligible=(y.sum(axis=1)>=10)&((y>0).sum(axis=1)>=3)
    contrast=[0,1]+[0]*(len(names)-2)
    directory=out/case['cohort'];directory.mkdir()
    table(directory/'counts.tsv',['featureID']+ids,([g]+row.tolist() for g,row in zip(features,y)))
    table(directory/'design.tsv',['sampleID']+names,([d]+row.tolist() for d,row in zip(ids,x)))
    table(directory/'samples.tsv',['sampleID','donor','condition','batch','libraryCounts','sizeFactor'],
        ([d,d,g['condition'],g['batchIDs'][0],int(l),float(f)] for d,g,l,f in zip(ids,observations,libraries,factors)))
    meta=dict(cohort=case['cohort'],sourceCountsSHA256=audit['referenceCountsSHA256'],sourceReportSHA256=check['compressedReportSHA256'],
        protocolSHA256=check['protocolSHA256'],minimumFeatureCounts=10,minimumExpressingPseudobulks=3,
        eligibleFeatures=int(eligible.sum()),attemptedFeatures=len(features),contrast=contrast,
        columnNames=names,design=x.tolist(),sourcePseudobulkIndices=indices,donorIDs=ids,
        sizeFactors=factors.tolist(),libraryCounts=libraries.tolist(),referenceFeatureIndices=reference.tolist(),
        normalization='Independent NumPy reconstruction of native all-positive median-ratio factors; native output checked separately')
    meta['files']={n:sha(directory/n) for n in ['counts.tsv','design.tsv','samples.tsv']}
    (directory/'input.json').write_text(json.dumps(meta,sort_keys=True,indent=2)+'\n')
    print(case['cohort'],'eligible',int(eligible.sum()),'reference',len(reference),'residualDF',len(ids)-len(names),flush=True)
