"""Compare every file-backed QC row and group against source, prior counts and the resident native owner."""
import json,os
import numpy as np
from parse_support import *

def compare(source,work,independent):
 new=work/'native-file';old=work/'native-resident';report=json.loads((new/'report.json').read_text())
 receipt=json.loads((new/'receipt.json').read_text());legacy=json.loads((old/'receipt.json').read_text())
 assert bytes(receipt['stream']['bytes'])==bytes(legacy['stream']['bytes']) and receipt['streamBytes']==legacy['streamBytes']==RECORDS*16
 for path,key in [('report.json','report'),('quality.bin','quality')]:assert sha(new/path)==bytes(receipt[key]['bytes']).hex()
 for path,key in [('plan.json','plan'),('report.json','report')]:assert sha(old/path)==bytes(legacy[key]['bytes']).hex()
 assert report['cellCount']==CELLS and report['canonicalNonzeros']==RECORDS and len(report['groups'])==24
 assert report['matrix']==object_value(old/'report.json','matrix')
 assert report['featureIDs']==list(array_items(old/'report.json','featureIDs'))
 group_ids=[g['biologicalReplicateID']+'_'+g['condition'] for g in report['groups']]
 assert len(set(group_ids))==24
 with np.load(independent) as measured,np.load(source/'donor-complete-counts.npz') as before:
  names=measured['groups'].tolist();counts=measured['counts'];previous_names=before['groups'].tolist();previous_counts=before['counts']
  assert counts.shape==previous_counts.shape==(24,FEATURES)
  matrix=report['matrix']
  for i,name in enumerate(group_ids):
   start,end=matrix['rowOffsets'][i:i+2];actual=np.zeros(FEATURES,dtype=np.uint64)
   actual[matrix['featureIndices'][start:end]]=matrix['counts'][start:end]
   assert np.array_equal(actual,counts[names.index(name)]) and np.array_equal(actual,previous_counts[previous_names.index(name)])
 legacy_quality=array_items(old/'report.json','quality');members=[0]*24;total_count=0;cells=0
 with (new/'quality.bin').open('rb') as f:
  for expected in source_rows(source):
   raw=f.read(QC.size);assert len(raw)==QC.size;total,mitochondrial,nnz,group=QC.unpack(raw);prior=next(legacy_quality)
   assert total==expected['totalCounts']==prior['totalCounts'] and nnz==expected['nonzeros']==prior['detectedFeatures']
   assert prior['barcode'].encode()==expected['barcode'].encode() and prior['sampleID'].encode()==expected['sampleID'].encode()
   assert mitochondrial==prior['mitochondrialCounts']==prior['mitochondrialFeatureCount']==0 and prior.get('mitochondrialFraction') is None
   assert group_ids[group]==expected['sampleID'];members[group]+=1;total_count+=total;cells+=1
  assert f.read(1)==b'' and next(legacy_quality,None) is None
 # The legacy report's per-group membership list is read one group at a time;
 # check every index against the file-backed ordinal without a whole-cohort map.
 keys=['biologicalReplicateID','donorID','condition','organism','cellGroup','sampleIDs','batchIDs'];groups_checked=0;indices_checked=0
 with (new/'quality.bin').open('rb') as f:
  for i,group in enumerate(array_items(old/'report.json','groups',maximum_item_chars=2097152)):
   assert all(group.get(k)==report['groups'][i].get(k) for k in keys)
   indices=group['sourceCellIndices'];assert len(indices)==members[i]==report['groups'][i]['sourceCellCount'];previous=-1
   for row in indices:
    assert previous<row<CELLS;previous=row
    assert QC.unpack(os.pread(f.fileno(),QC.size,row*QC.size))[3]==i
    indices_checked+=1
   groups_checked+=1
 assert cells==indices_checked==CELLS and groups_checked==24 and total_count==3070817047
 return dict(status='passed',cells=CELLS,groups=24,features=FEATURES,records=RECORDS,totalCounts=total_count,
             allQCAndOriginalIdentityBytesMatch=True,allGroupMembersMatch=True,allAggregateCoordinatesMatchPriorAndFreshIndependentCounts=True,
             fileAndResidentStreamHashesMatch=True,predictionFitOrScoring=False)
