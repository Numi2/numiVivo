from pathlib import Path
import os,sys,json,hashlib,subprocess,time,fcntl,concurrent.futures
import numpy as np
source=Path('/Users/n/numivivo-parse-bcell-counts-20260912')
os.environ['NUMIVIVO_PARSE_STUDY']=str(source);sys.path.insert(0,str(source))
import run_bcells
s=Path('/Users/n/numivivo-parse-logcpm-20260912')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
lock=(source/'controller.lock').open('r');fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
binary=Path('/Users/n/numivivo-streamed-logcpm-20260912/stream-check')
manifest=json.loads((source/'donor-plans/manifest.json').read_text())
for name,h in manifest['boundSourceFiles'].items():assert sha(source/name)==h
freeze={'driverSHA256':sha(s/'run.py'),'maximumAbsoluteErrorTolerance':1e-10,'binarySHA256':sha(binary),'manifestSHA256':sha(source/'donor-plans/manifest.json'),'adapterSHA256':{name:sha(source/name) for name in ['run_bcells.py','run_counts.py']},'endpoint':'mean per-cell log1p CPM, all 40352 features; pooled frozen B-cell membership','predictionFitted':False,'predictionScored':False}
write(s/'protocol.json',freeze)
completed=[]
for part in manifest['parts']:
 donor=part['donor'];d=s/donor;d.mkdir()
 old=json.loads((source/'execution/ingest'/donor/'published.json').read_text())
 planpath=source/'donor-plans'/donor/'plan.json';assert sha(planpath)==part['planSHA256']
 plan=json.loads(planpath.read_text());ind=source/old['independentCounts'];assert sha(ind)==old['independentCountsSHA256']
 totals=np.load(ind)['cellTotals']
 samples={x['id']:x['condition'] for x in plan['metadata']['samples']}
 groups=sorted(set(samples.values()));assert groups==['IFN-beta','PBS']
 assignment=np.array([groups.index(samples[c['sampleID']]) for c in plan['metadata']['cells']])
 p=dict(featureIDs=[x['id'] for x in plan['metadata']['features']],groupIDs=groups,rowGroups=assignment.tolist(),rowTotals=totals.tolist())
 write(d/'plan.json',p)
 runsfile=source/'donor-plans'/donor/'runs.json';assert sha(runsfile)==part['runsSHA256'];runs=json.loads(runsfile.read_text())
 sums=np.zeros((2,40352));digest=hashlib.sha256();received=0
 log=(d/'native.log').open('wb');proc=subprocess.Popen([str(binary),str(d/'plan.json'),str(d/'native.json')],stdin=subprocess.PIPE,stdout=log,stderr=subprocess.STDOUT)
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
   pending={};submitted=0
   for wanted in range(len(runs)):
    while submitted<min(len(runs),wanted+8):
     pending[submitted]=pool.submit(run_bcells.one,runs[submitted]);submitted+=1
    run,data,bulk,rowtotals,qc=pending.pop(wanted).result()
    prior=json.loads((source/old['attempt']/("run-%04d.json"%run['run'])).read_text())
    assert [(v['start'],v['stop'],v['SHA256']) for v in qc['ranges']]==[(v['start'],v['stop'],v['SHA256']) for v in prior['ranges']]
    np.testing.assert_array_equal(rowtotals,totals[run['selectedLo']:run['selectedHi']])
    records=np.frombuffer(data,dtype='<u8').reshape(-1,2);rows=(records[:,0]&np.uint64(0xffffffff)).astype(int);cols=(records[:,0]>>np.uint64(32)).astype(int)
    vals=np.log1p(records[:,1].astype(float)/totals[rows]*1e6)
    for g in range(2):
     take=assignment[rows]==g;sums[g]+=np.bincount(cols[take],weights=vals[take],minlength=40352)
    proc.stdin.write(data);digest.update(data);received+=len(data)
    write(d/("run-%04d.json"%run['run']),qc)
    if wanted%25==0:print(donor,wanted+1,len(runs),flush=True)
  proc.stdin.close();assert proc.wait()==0
  assert digest.hexdigest()==old['streamSHA256'] and received==old['streamBytes']
  a=json.loads((d/'native.json').read_text());sizes=np.bincount(assignment,minlength=2);ref=sums/sizes[:,None]
  assert a['featureIDs']==p['featureIDs'] and a['groupIDs']==groups and a['cellCounts']==sizes.tolist()
  assert a['zeroCellCounts']==np.bincount(assignment[totals==0],minlength=2).tolist()
  np.testing.assert_allclose(a['means'],ref,rtol=0,atol=1e-10)
  np.save(d/'reference.npy',ref)
  rec={'donor':donor,'status':'PASS','cells':len(totals),'maximumAbsoluteError':float(np.max(abs(ref-np.array(a['means'])))),'streamSHA256':digest.hexdigest(),'nativeSHA256':sha(d/'native.json'),'planSHA256':sha(d/'plan.json'),'referenceSHA256':sha(d/'reference.npy')}
  write(d/'verification.json',rec);completed.append(rec);print(rec,flush=True)
 except BaseException as error:
  if proc.poll() is None:proc.kill()
  proc.wait();write(d/'failure.json',{'type':type(error).__name__,'error':str(error),'streamBytes':received});raise
 finally:log.close()
write(s/'complete.json',{'status':'PASS-all-donors','donors':completed,'cells':sum(x['cells'] for x in completed),'predictionFitted':False,'predictionScored':False})
