#!/usr/bin/env python3
"""Independently reconstruct author QC decisions and verify native selected counts."""
import argparse,csv,gzip,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy import sparse
from prepare_ingestion import sha,write
from prepare_cohort import RESTORED,AUTHOR,AUTHOR_URL
from restore_geo_identities import BARCODES,GUIDES

p=argparse.ArgumentParser(description=__doc__)
for key in ['cohort','barcodes','guides','full_report','full_reference','descriptors','out']:
    p.add_argument('--'+key.replace('_','-'),type=Path,required=True)
a=p.parse_args();assert not a.out.exists()
assert sha(a.barcodes)==BARCODES and sha(a.guides)==GUIDES
with gzip.open(a.barcodes,'rt') as f:barcodes=[row[0] for row in csv.reader(f,delimiter='\t')]
with gzip.open(a.guides,'rt') as f:
    reader=csv.DictReader(f);index=reader.fieldnames[0];records={r[index]:r for r in reader}
selected=[];excluded=[];labels=[]
for i,barcode in enumerate(barcodes):
    r=records.get(barcode)
    if r is None:reasons=['missing-original-guide-record'];label=None
    else:
        label=r['guide identity'];reasons=[]
        if int(r['number of cells'])!=1:reasons.append('not-one-assigned-identity')
        assert r['good coverage'] in ['TRUE','FALSE']
        if r['good coverage']=='FALSE':reasons.append('not-good-guide-coverage')
        if label=='*':reasons.append('ambiguous-guide-star')
    if reasons:excluded.append(dict(sourceRow=i,barcode=barcode,guide=label,reasons=reasons))
    else:selected.append(i);labels.append(label)
c=json.loads((a.cohort/'cohort.json').read_text());plan=json.loads((a.cohort/'plan.json').read_text())
assert c['selectedRows']==selected and c['excluded']==excluded
assert c['barcodes']==[barcodes[i] for i in selected]
assert plan['cellSelection']['observationIndices']==selected and bytes(plan['cellSelection']['source']['bytes']).hex()==RESTORED
ref=np.load(a.cohort/'reference.npz');assert np.array_equal(ref['selectedRows'],selected)
r=json.loads((a.cohort/'native/report.json').read_text());full=json.loads(a.full_report.read_text())
assert r['sourceCellCount']==len(barcodes) and r['sourceObservationIndices']==selected
assert r['metadata']['cells']==[full['metadata']['cells'][i] for i in selected]
assert r['metadata']['features']==full['metadata']['features']
assert r['metadata']['samples']==full['metadata']['samples']
assert r['quality']==[full['quality'][i] for i in selected]
mat=r['pseudobulk']['matrix'];groups=r['pseudobulk']['groups']
actual=sparse.csr_matrix((np.array(mat['counts'],dtype=np.int64),mat['featureIndices'],mat['rowOffsets']),shape=(mat['cellCount'],mat['featureCount']))
assert np.array_equal(actual.toarray(),ref['counts'])
assert r['pseudobulk']['featureIDs']==c['features']
assert [g['condition'] for g in groups]==c['groups']
for g in groups:
    members=[i for i,label in enumerate(labels) if label==g['condition']]
    assert g['sourceCellIndices']==members
    assert g['sampleIDs']==sorted({r['metadata']['cells'][i]['sampleID'] for i in members})
    assert g['batchIDs']==sorted({'GSM2406681-GEM-'+barcodes[selected[i]].rsplit('-',1)[1] for i in members})
    assert g['biologicalReplicateID']=='K562-pooled-replication-unresolved' and 'donorID' not in g
assert sorted(i for g in groups for i in g['sourceCellIndices'])==list(range(len(selected)))
q=np.load(a.full_reference);totals=q['totalCounts']
assert sum(x['totalCounts'] for x in r['quality'])==int(ref['counts'].sum())==int(totals[selected].sum())
assert sum(x['detectedFeatures'] for x in r['quality'])==r['canonicalNonzeros']
assert int(totals[selected].sum())+int(totals[[x['sourceRow'] for x in excluded]].sum())==1039857798
candidates=json.loads((a.descriptors/'candidates.json').read_text())['candidates']
retained=[x['target'] for x in candidates if set(x['sourceGuides'])&set(labels)]
removed=[x['target'] for x in candidates if not set(x['sourceGuides'])&set(labels)]
write(a.out,dict(status='passed-native-author-cohort',sourceCells=len(barcodes),selectedCells=len(selected),excludedCells=len(excluded),
    selectedGuideGroups=len(groups),selectedUMIs=int(ref['counts'].sum()),excludedUMIs=int(totals[[x['sourceRow'] for x in excluded]].sum()),
    selectedNonzeros=r['canonicalNonzeros'],allSelectedMetadataQualityCountsAndMembershipExact=True,
    independentlyReconstructedExclusions=True,retainedCandidatePrefixes=retained,removedCandidatePrefixes=removed,
    originalSourceSHA256=RESTORED,reportSHA256=sha(a.cohort/'native/report.json'),planSHA256=sha(a.cohort/'plan.json'),
    authorCodeURL=AUTHOR_URL,authorCodeSHA256=AUTHOR,controlsVerified=False,predictionsFitted=False))
print(a.out.read_text())
