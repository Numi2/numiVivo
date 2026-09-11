"""Run/replay complete native donor streams; resume only validated donor receipts."""
from pathlib import Path
import argparse,json,hashlib,concurrent.futures,subprocess,time,os,platform
import run_counts as adapter
ROOT=adapter.ROOT;np=adapter.np

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def check_bundle(bundle,plan_path,counts_path):
 receipt=json.loads((bundle/'receipt.json').read_text());plan=json.loads((bundle/'plan.json').read_text());report=json.loads((bundle/'report.json').read_text());assert plan==json.loads(plan_path.read_text())
 for key in ['plan','report']:assert sha(bundle/(key+'.json'))==bytes(receipt[key]['bytes']).hex()
 with np.load(counts_path) as stored:expected={key:stored[key] for key in stored.files}
 assert len(report['quality'])==len(plan['metadata']['cells'])
 for i,(cell,qc) in enumerate(zip(plan['metadata']['cells'],report['quality'])):
  assert cell['sampleID']==qc['sampleID'] and cell['barcode']==qc['barcode'] and int(expected['cellTotals'][i])==qc['totalCounts'] and plan['rowNonzeros'][i]==qc['detectedFeatures']
 bulk=report['pseudobulk'];assert bulk['featureIDs']==[f['id'] for f in plan['metadata']['features']];matrix=bulk['matrix'];assert len(bulk['groups'])==2
 sample=np.array([c['sampleID'] for c in plan['metadata']['cells']])
 for i,g in enumerate(bulk['groups']):
  name=g['biologicalReplicateID']+'_'+g['condition'];j=expected['groups'].tolist().index(name);lo,hi=matrix['rowOffsets'][i:i+2];values=np.zeros(adapter.FEATURES,dtype=np.uint64);values[matrix['featureIndices'][lo:hi]]=matrix['counts'][lo:hi]
  assert np.array_equal(values,expected['counts'][j]) and np.array_equal(g['sourceCellIndices'],np.flatnonzero(sample==name))
 assert report['canonicalNonzeros']==sum(plan['rowNonzeros']);return receipt

def main():
 p=argparse.ArgumentParser();p.add_argument('--phase',choices=['ingest','verify'],required=True);a=p.parse_args();frozen=json.loads((ROOT/'donor-plans/manifest.json').read_text());assert frozen['status']=='frozen-disjoint-complete-donor-partitions' and len(frozen['parts'])==12
 binary=ROOT/'build-pooled/numivivo-omics';binary_sha=sha(binary);phase_root=ROOT/'donor-execution'/a.phase;phase_root.mkdir(parents=True,exist_ok=True);all_completed=[]
 for part in frozen['parts']:
  donor=part['donor'];plans=ROOT/'donor-plans'/donor;plan_path=plans/'plan.json';assert sha(plan_path)==part['planSHA256'];runs=json.loads((plans/'runs.json').read_text());directory=phase_root/donor;directory.mkdir(exist_ok=True);pointer=directory/'published.json'
  baseline=None
  if a.phase=='verify':
   baseline=json.loads((ROOT/'donor-execution/ingest'/donor/'published.json').read_text());assert baseline['status']=='passed' and baseline['binarySHA256']==binary_sha
  if pointer.exists():
   completed=json.loads(pointer.read_text());assert completed['status']=='passed' and completed['binarySHA256']==binary_sha and completed['planSHA256']==part['planSHA256']
   bundle=ROOT/completed['bundle'];receipt=check_bundle(bundle,plan_path,ROOT/completed['independentCounts']);assert bytes(receipt['stream']['bytes']).hex()==completed['streamSHA256'] and sha(bundle/'receipt.json')==completed['receiptSHA256'] and sha(ROOT/completed['independentCounts'])==completed['independentCountsSHA256'];all_completed.append(completed);print('RESUMED',donor,a.phase,flush=True);continue
  attempt=directory/f'attempt-{len(list(directory.glob("attempt-*")))+1:03d}';attempt.mkdir(exist_ok=False);bundle=ROOT/baseline['bundle'] if baseline else attempt/'native';command=[str(binary),'singlecell-count-stream-verify',str(bundle)] if baseline else [str(binary),'singlecell-count-stream-pseudobulk','--plan',str(plan_path),'--output',str(bundle)]
  result=dict(status='running',phase=a.phase,donor=donor,binarySHA256=binary_sha,planSHA256=part['planSHA256'],command=command,platform=platform.platform(),cells=part['cells'],records=part['records'],bundle=str(bundle.relative_to(ROOT)),attempt=str(attempt.relative_to(ROOT)));write(attempt/'execution.json',result)
  groups=sorted({r['sample'] for r in runs});matrix=np.zeros((2,adapter.FEATURES),dtype=np.uint64);totals=np.zeros(part['cells'],dtype=np.uint64);digest=hashlib.sha256();received=0;before=time.monotonic();stdout=(attempt/'native.stdout').open('wb');stderr=(attempt/'native.stderr').open('wb');process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=stdout,stderr=stderr)
  try:
   with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
    pending={};next_submit=0
    for wanted in range(len(runs)):
     while next_submit<min(len(runs),wanted+16):pending[next_submit]=pool.submit(adapter.one,runs[next_submit]);next_submit+=1
     run,data,bulk,rowtotals,qc=pending.pop(wanted).result();index=run['run']
     prior=(ROOT/baseline['attempt']/f'run-{index:04d}.json') if baseline else ROOT/f'ingest/run-{index:04d}.json'
     if prior.exists():
      old=json.loads(prior.read_text());assert [(x['start'],x['stop'],x['SHA256']) for x in qc['ranges']]==[(x['start'],x['stop'],x['SHA256']) for x in old['ranges']]
     elif baseline:raise AssertionError('missing original donor source-range receipt')
     process.stdin.write(data);digest.update(data);received+=len(data);matrix[groups.index(run['sample'])]+=bulk;totals[run['lo']:run['hi']]=rowtotals
     write(attempt/f'run-{index:04d}.json',qc)
     if wanted%25==0:print('DONOR-STREAM',a.phase,donor,wanted,len(runs),run['hi'],received,flush=True)
   process.stdin.close();_,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status);stdout.close();stderr.close();assert process.returncode==0,(attempt/'native.stderr').read_text()[-2000:]
   counts_path=attempt/'independent-counts.npz';np.savez_compressed(counts_path,groups=np.array(groups),counts=matrix,cellTotals=totals);receipt=check_bundle(bundle,plan_path,counts_path)
   assert received==part['records']*16==receipt['streamBytes'] and digest.hexdigest()==bytes(receipt['stream']['bytes']).hex()
   if baseline:assert digest.hexdigest()==baseline['streamSHA256']
   result.update(status='passed',seconds=time.monotonic()-before,streamSHA256=digest.hexdigest(),streamBytes=received,totalCounts=int(totals.sum()),maximumNativeRSS=usage.ru_maxrss,nativeExit=0,allCellAndAggregateCountsVerified=True,independentCounts=str(counts_path.relative_to(ROOT)),independentCountsSHA256=sha(counts_path),receiptSHA256=sha(bundle/'receipt.json'))
   write(attempt/'execution.json',result);temporary=pointer.with_suffix('.tmp');write(temporary,result);os.replace(temporary,pointer);all_completed.append(result);print('DONOR-COMPLETE',a.phase,donor,result['cells'],result['records'],flush=True)
  except BaseException as error:
   result.update(status='failed',errorType=type(error).__name__,error=str(error),streamedBytes=received,seconds=time.monotonic()-before)
   if process.stdin and not process.stdin.closed:process.stdin.close()
   process.wait();write(attempt/'execution.json',result);raise
  finally:stdout.close();stderr.close()
 assert len(all_completed)==12 and sum(x['cells'] for x in all_completed)==725031 and sum(x['records'] for x in all_completed)==1373870697
 write(phase_root/'complete.json',dict(status='passed',phase=a.phase,binarySHA256=binary_sha,sourceManifestSHA256=sha(ROOT/'donor-plans/manifest.json'),donors=all_completed,cells=725031,records=1373870697,scope='Twelve disjoint complete native donor streams. Each retains its own canonical stream SHA256, full cell QC and independent aggregate verification.'))
 print('COMPLETE COHORT',a.phase,725031,1373870697,flush=True)
if __name__=='__main__':main()
