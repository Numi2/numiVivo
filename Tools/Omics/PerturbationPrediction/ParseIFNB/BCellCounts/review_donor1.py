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
r=read(root/'execution/ingest/Donor1/published.json');p=root/'donor-plans/Donor1'
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
out=dict(status='passed-Donor1-ingest-offline-review',createdUnix=time.time(),donor='Donor1',phase='ingest',sourceRuns=len(runs),rangeReceipts=ranges,totals=totals,streamSHA256=r['streamSHA256'],preparationManifestSHA256=sha(Path('/Users/n/numivivo-bcell-preparation-manifest.json')),completeCohort=False,sourceReplay=False,predictionScored=False)
Path('/Users/n/numivivo-bcell-donor1-review.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

