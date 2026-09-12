from pathlib import Path
import os,sys,json,hashlib,time
root=Path('/Users/n/numivivo-parse-bcell-counts-20260912')
sys.path.insert(0,str(root));os.environ['NUMIVIVO_PARSE_STUDY']=str(root)
from run_donors import check_bundle
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
manifest=read(Path('/Users/n/numivivo-bcell-preparation-manifest.json'))
for name,h in manifest['members'].items():assert sha(root/name)==h,name
f=read(root/'donor-plans/manifest.json')
for name,h in f['boundSourceFiles'].items():assert sha(root/name)==h,name
checked=[];receipts=[]
for part in f['parts']:
 donor=part['donor']
 r=read(root/'execution/ingest'/donor/'published.json');p=root/'donor-plans'/donor
 assert read(root/r['attempt']/'execution.json')==r
 assert sha(root/r['independentCounts'])==r['independentCountsSHA256']
 assert sha(root/r['bundle']/'receipt.json')==r['receiptSHA256']
 receipt=check_bundle(root/r['bundle'],p/'plan.json',root/r['independentCounts'])
 assert receipt['streamBytes']==r['streamBytes']==r['records']*16
 assert bytes(receipt['stream']['bytes']).hex()==r['streamSHA256']
 keys=['selectedCells','selectedNonzeros','selectedTotalCounts','sourceGeneCountMismatch','sourceTscpCountMismatch','sourceMinusMatrixTotal','sourceMinusMatrixGeneCount']
 totals=dict.fromkeys(keys,0);ranges=0
 runs=read(p/'runs.json')
 for run in runs:
  qc=read(root/r['attempt']/('run-%04d.json'%run['run']))
  assert qc['selectedCells']==run['selectedHi']-run['selectedLo']
  assert qc['ranges']
  ranges+=len(qc['ranges'])
  for k in keys:totals[k]+=qc[k]
 assert totals['selectedCells']==r['cells'] and totals['selectedNonzeros']==r['records'] and totals['selectedTotalCounts']==r['totalCounts']
 assert r['status']=='passed' and r['phase']=='ingest' and r['donor']==donor
 assert r['cells']==part['cells'] and r['records']==part['records']
 assert sha(p/'plan.json')==r['planSHA256']==part['planSHA256']
 assert sha(p/'runs.json')==part['runsSHA256']
 assert r['binarySHA256']==f['binarySHA256']
 checked.append(dict(donor=donor,totals=totals,sourceRuns=len(runs),rangeReceipts=ranges,streamSHA256=r['streamSHA256']))
 receipts.append(r)
completion=read(root/'ingest-complete.json')
assert completion==dict(status='passed',parts=receipts,cells=72446,records=124909573,manifestSHA256=sha(root/'donor-plans/manifest.json'))
assert len(checked)==12 and len({x['donor'] for x in checked})==12
aggregate={k:sum(x['totals'][k] for x in checked) for k in keys}
assert aggregate['selectedCells']==72446 and aggregate['selectedNonzeros']==124909573
out=dict(status='passed-all-donor-ingestion-offline-review',createdUnix=time.time(),donors=checked,totals=aggregate,sourceRuns=sum(x['sourceRuns'] for x in checked),rangeReceipts=sum(x['rangeReceipts'] for x in checked),completeIngestionCohort=True,sourceReplay=False,predictionScored=False,preparationManifestSHA256=sha(Path('/Users/n/numivivo-bcell-preparation-manifest.json')))
Path('/Users/n/numivivo-bcell-all-ingest-review.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='donors'},indent=2))
