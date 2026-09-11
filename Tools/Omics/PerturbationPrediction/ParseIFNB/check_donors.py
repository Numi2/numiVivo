"""Verify complete disjoint source coverage and all native ingest/replay outputs."""
from pathlib import Path
import os,json,hashlib,numpy as np
from run_donors import check_bundle
ROOT=Path(os.environ['NUMIVIVO_PARSE_STUDY']);frozen=json.loads((ROOT/'donor-plans/manifest.json').read_text());parent=json.loads((ROOT/'prepared/plan.json').read_text());seen=np.zeros(725031,dtype=bool);totals=np.zeros(725031,dtype=np.uint64);groups=[];counts=[];run_ids=[];qc_sums={k:0 for k in ['sourceGeneCountMismatch','sourceTscpCountMismatch','sourceMinusMatrixTotal','sourceMinusMatrixGeneCount','sortedRows']};receipts=[]
for phase in ['ingest','verify']:
 complete=json.loads((ROOT/'donor-execution'/phase/'complete.json').read_text());assert complete['status']=='passed' and complete['cells']==725031 and complete['records']==1373870697
for part in frozen['parts']:
 donor=part['donor'];plans=ROOT/'donor-plans'/donor;mapping=np.load(plans/'global-selected-rows.npy');assert len(mapping)==part['cells'] and len(np.unique(mapping))==len(mapping) and not seen[mapping].any();seen[mapping]=True
 plan=json.loads((plans/'plan.json').read_text());assert plan['metadata']['features']==parent['metadata']['features']
 assert all(plan['metadata']['cells'][i]==parent['metadata']['cells'][int(global_row)] and plan['rowNonzeros'][i]==parent['rowNonzeros'][int(global_row)] for i,global_row in enumerate(mapping))
 entries=[]
 for phase in ['ingest','verify']:
  entry=json.loads((ROOT/'donor-execution'/phase/donor/'published.json').read_text());assert entry['status']=='passed' and entry['planSHA256']==part['planSHA256'];receipt=check_bundle(ROOT/entry['bundle'],plans/'plan.json',ROOT/entry['independentCounts']);assert hashlib.sha256((ROOT/entry['bundle']/'receipt.json').read_bytes()).hexdigest()==entry['receiptSHA256'];assert hashlib.sha256((ROOT/entry['independentCounts']).read_bytes()).hexdigest()==entry['independentCountsSHA256'];assert bytes(receipt['stream']['bytes']).hex()==entry['streamSHA256'];entries.append(entry)
 assert entries[0]['streamSHA256']==entries[1]['streamSHA256'] and entries[0]['binarySHA256']==entries[1]['binarySHA256']
 first=np.load(ROOT/entries[0]['independentCounts']);second=np.load(ROOT/entries[1]['independentCounts']);assert set(first.files)==set(second.files)
 for key in first.files:assert np.array_equal(first[key],second[key])
 totals[mapping]=first['cellTotals'];groups+=first['groups'].tolist();counts += list(first['counts']);receipts.append(dict(donor=donor,streamSHA256=entries[0]['streamSHA256'],streamBytes=entries[0]['streamBytes'],ingestPeakRSS=entries[0]['maximumNativeRSS'],verifyPeakRSS=entries[1]['maximumNativeRSS']))
 for run in json.loads((plans/'runs.json').read_text()):
  run_ids.append(run['run']);old=json.loads((ROOT/entries[0]['attempt']/f'run-{run["run"]:04d}.json').read_text());new=json.loads((ROOT/entries[1]['attempt']/f'run-{run["run"]:04d}.json').read_text())
  assert [(x['start'],x['stop'],x['SHA256']) for x in old['ranges']]==[(x['start'],x['stop'],x['SHA256']) for x in new['ranges']]
  for key in qc_sums:qc_sums[key]+=old[key]
assert seen.all() and sorted(run_ids)==list(range(3456)) and len(set(groups))==24 and qc_sums['sourceGeneCountMismatch']==30634
np.savez_compressed(ROOT/'donor-complete-counts.npz',groups=np.array(groups),counts=np.array(counts,dtype=np.uint64),cellTotals=totals)
result=dict(status='passed',cells=725031,features=40352,groups=24,donors=12,records=1373870697,totalRetainedCounts=int(totals.sum()),allCellAndAggregateCountsMatch=True,allSourceRowsCoveredExactlyOnce=True,allDonorNativeVerificationsPassed=True,allRangeHashesAndDonorStreamHashesMatch=True,sourceQCDiscrepancies=qc_sums,streams=receipts,sourceFingerprintScope='Pinned ETag and exact fetched ranges; twelve donor stream SHA256 values, not a whole-source-file or monolithic-stream SHA256',predictionFitOrScoring=False)
(ROOT/'donor-counts-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
