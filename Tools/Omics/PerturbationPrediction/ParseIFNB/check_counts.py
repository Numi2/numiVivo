import os
from pathlib import Path
import json,hashlib,numpy as np,sys
ROOT=Path(os.environ['NUMIVIVO_PARSE_STUDY']);checks={};receipt=json.loads((ROOT/'native/receipt.json').read_text());report=json.loads((ROOT/'native/report.json').read_text());plan=json.loads((ROOT/'native/plan.json').read_text());reference=np.load(ROOT/'ingest/independent-counts.npz');freeze=json.loads((ROOT/'prepared/freeze.json').read_text())
for name in ['plan','report']:
 digest=hashlib.sha256((ROOT/f'native/{name}.json').read_bytes()).digest();assert list(digest)==receipt[name]['bytes']
assert plan==json.loads((ROOT/'prepared/plan.json').read_text())
metadata=plan['metadata'];quality=report['quality'];assert len(quality)==len(metadata['cells'])==freeze['cells']==725031
for i,(cell,qc) in enumerate(zip(metadata['cells'],quality)):
 assert cell['barcode']==qc['barcode'] and cell['sampleID']==qc['sampleID']
 assert int(reference['cellTotals'][i])==qc['totalCounts'] and plan['rowNonzeros'][i]==qc['detectedFeatures']
 assert qc['mitochondrialFeatureCount']==0 and qc.get('mitochondrialFraction') is None
bulk=report['pseudobulk'];matrix=bulk['matrix'];assert bulk['featureIDs']==[x['id'] for x in metadata['features']];assert len(bulk['groups'])==24
sample=np.array([x['sampleID'] for x in metadata['cells']]);assert report['canonicalNonzeros']==freeze['nonzeros']==1373870697
for i,group in enumerate(bulk['groups']):
 name=group['biologicalReplicateID']+'_'+group['condition'];idx=reference['groups'].tolist().index(name);start,end=matrix['rowOffsets'][i:i+2];values=np.zeros(40352,dtype=np.uint64);values[matrix['featureIndices'][start:end]]=matrix['counts'][start:end]
 assert np.array_equal(values,reference['counts'][idx]);assert np.array_equal(group['sourceCellIndices'],np.flatnonzero(sample==name))
counts=np.zeros(5,dtype=np.int64);bytes_read=0;requests=0
for i in range(3456):
 check=json.loads((ROOT/f'ingest/run-{i:04d}.json').read_text());counts+=np.array([check[k] for k in ['sourceGeneCountMismatch','sourceTscpCountMismatch','sourceMinusMatrixTotal','sourceMinusMatrixGeneCount','sortedRows']]);bytes_read+=sum(x['stop']-x['start']+1 for x in check['ranges']);requests+=len(check['ranges'])
result=dict(status='passed',scope='All selected cells and retained features; source version and fetched ranges, not an entire-source-file SHA256 or validation of unselected count values',cells=725031,features=40352,donorConditionGroups=24,canonicalNonzeros=report['canonicalNonzeros'],totalRetainedCounts=int(reference['cellTotals'].sum()),allCellIdentitiesTotalsAndNonzerosMatch=True,all24By40352AggregateCountsMatch=True,allGroupMembershipsMatch=True,sourceGeneCountMismatchedCells=int(counts[0]),sourceTscpCountMismatchedCells=int(counts[1]),sourceMinusRetainedTotalCounts=int(counts[2]),sourceMinusRetainedGeneCounts=int(counts[3]),rowsCanonicalized=int(counts[4]),countHTTPBytes=bytes_read,countHTTPRequests=requests,predictionFitOrScoring=False,mitochondrialFractions='unavailable: no source annotation admitted')
if (ROOT/'replay/execution.json').exists():
 replay=json.loads((ROOT/'replay/execution.json').read_text());assert replay['status']=='passed';other=np.load(ROOT/'replay/independent-counts.npz');assert set(other.files)==set(reference.files)
 for key in reference.files:assert np.array_equal(reference[key],other[key])
 assert replay['streamSHA256']==json.loads((ROOT/'ingest/execution.json').read_text())['streamSHA256'];result['completeNativeReplay']='passed'
(ROOT/'counts-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
