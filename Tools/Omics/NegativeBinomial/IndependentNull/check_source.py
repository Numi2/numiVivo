#!/usr/bin/env python3
"""Check every native cell and pseudobulk against the independently audited source."""
import argparse, gzip, hashlib, json
from pathlib import Path
import h5py, numpy as np
from anndata.io import read_elem

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
root=p.parse_args().root
sha=lambda b:hashlib.sha256(b).hexdigest()
read=lambda n:json.loads((root/n).read_text())
audit=read('count-audit.json');plan=read('stream-plan.json');receipt=read('native-receipt.json')
cohorts=read('cohorts.json')
assert cohorts['protocolSHA256']==sha(Path(__file__).with_name('PROTOCOL.md').read_bytes())
compressed=(root/'source-report.json.gz').read_bytes();raw=gzip.decompress(compressed);report=json.loads(raw)
assert bytes(receipt['source']['bytes']).hex()==audit['sourceSHA256']==cohorts['source']['sha256']
assert bytes(receipt['report']['bytes']).hex()==sha(raw)
assert bytes(receipt['plan']['bytes']).hex()==sha(json.dumps(plan,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode())
assert sha((root/'reference-counts.npz').read_bytes())==audit['referenceCountsSHA256']
assert sha((root/'reference-qc.npz').read_bytes())==audit['referenceQCSHA256']
# These are our local audit products, bound to exact hashes above. Pandas string
# axes were serialized as object arrays; never enable pickle for source downloads.
with np.load(root/'reference-counts.npz',allow_pickle=True) as archive:
    counts={key:archive[key] for key in archive.files}
with np.load(root/'reference-qc.npz',allow_pickle=True) as archive:
    qc={key:archive[key] for key in archive.files}
with h5py.File(root/'source.h5ad','r') as h:
    obs=read_elem(h['obs']);var=read_elem(h['raw/var'])
metadata=report['metadata'];bulk=report['pseudobulk'];mapping=plan['mapping']
assert report['contrasts']==[] and report['canonicalNonzeros']==audit['canonicalNonzeros']
assert report.get('sourceObservationIndices') is None and report.get('sourceCellCount') is None
assert metadata['countUnit']==bulk['countUnit']=='umiCount' and metadata['evidence']=='measured'
assert metadata['samples']==mapping['samples']
assert len(metadata['cells'])==len(report['quality'])==len(obs)==audit['cells']
assert bulk['featureIDs']==counts['features'].tolist()==var.index.astype(str).tolist()
expected_features=[dict(id=str(i),name=str(row.gene_symbols),mitochondrial=bool(row.mito)) for i,row in var.iterrows()]
assert metadata['features']==expected_features
donors=counts['donors'].tolist();donor_index={d:i for i,d in enumerate(donors)}
sample_by_id={s['id']:s for s in mapping['samples']};mt=int(var.mito.sum())
assert qc['barcodes'].tolist()==obs.index.astype(str).tolist()
for i,(cell,q) in enumerate(zip(metadata['cells'],report['quality'])):
    donor=str(obs.donor_id.iloc[i]);barcode=str(obs.index[i])
    assert cell==dict(barcode=barcode,sampleID=donor,group='B cell')
    assert qc['sample'][i]==donor_index[donor]
    assert q['sampleID']==donor and q['barcode']==barcode
    for field,key in [('totalCounts','totals'),('detectedFeatures','detected'),('mitochondrialCounts','mitochondrialCounts')]:
        assert q[field]==int(qc[key][i]),(i,field)
    assert q['mitochondrialFeatureCount']==mt
    expected=float(qc['mitochondrialCounts'][i])/int(qc['totals'][i]) if qc['totals'][i] and mt else None
    assert q.get('mitochondrialFraction')==expected
matrix=bulk['matrix'];seen=set();membership=[]
assert matrix['cellCount']==len(bulk['groups'])==len(donors)==108 and matrix['featureCount']==len(var)
assert len(matrix['counts'])==len(matrix['featureIndices'])==audit['aggregateNonzeros']
assert len(matrix['rowOffsets'])==109 and matrix['rowOffsets'][0]==0 and matrix['rowOffsets'][-1]==audit['aggregateNonzeros']
for row,group in enumerate(bulk['groups']):
    donor=group['donorID'];assert donor not in seen;seen.add(donor);sample=sample_by_id[donor]
    assert group['biologicalReplicateID']==donor and group['condition']==sample['condition']
    assert group['batchIDs']==[sample['batchID']] and group['sampleIDs']==[donor]
    assert group['cellGroup']=='B cell' and group['organism']==sample['organism']
    expected=np.flatnonzero(qc['sample']==donor_index[donor]).tolist()
    assert group['sourceCellIndices']==expected;membership.extend(expected)
    lo,hi=matrix['rowOffsets'][row:row+2];cols=matrix['featureIndices'][lo:hi]
    reference=counts['counts'][donor_index[donor]]
    assert cols==np.flatnonzero(reference).tolist()
    assert matrix['counts'][lo:hi]==reference[cols].tolist()
assert sorted(membership)==list(range(len(obs))) and seen==set(donors)
result=dict(status='passed',sourceSHA256=audit['sourceSHA256'],reportSHA256=sha(raw),compressedReportSHA256=sha(compressed),
    protocolSHA256=cohorts['protocolSHA256'],planFileSHA256=sha((root/'stream-plan.json').read_bytes()),
    sourceCells=len(obs),sourceFeatures=len(var),donors=len(donors),canonicalNonzeros=report['canonicalNonzeros'],
    aggregateNonzeros=len(matrix['counts']),totalUMIs=sum(matrix['counts']),
    checked='Every original cell identity, QC value, raw feature, sample annotation, group membership and donor aggregate count; no cell-by-gene dense array')
(root/'source-check.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,indent=2))
