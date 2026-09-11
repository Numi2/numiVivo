"""Tee the complete real Parse count stream to file-backed and prior resident native owners."""
from pathlib import Path
import argparse,concurrent.futures,hashlib,importlib.util,json,os,platform,subprocess,sys,tarfile,time
import numpy as np
from parse_support import *
from compare_parse_counts import compare

OLD_BINARY_SHA='20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516'
def references(repository,work):
 archive=repository/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-counts/results.tar.gz'
 members=archive_members(archive,COUNTS_SHA);selected={}
 for name,identity in members.items():
  path=Path(name)
  if name.startswith('donor-execution/ingest/') and path.name.startswith('run-') and path.suffix=='.json':
   i=int(path.stem[4:]);assert i not in selected;selected[i]=identity
 assert sorted(selected)==list(range(3456))
 target=work/'expected-ranges'
 if target.exists():
  assert {p.name for p in target.iterdir()}=={f'run-{i:04d}.json' for i in selected}
  assert all(sha(target/f'run-{i:04d}.json')==item['SHA256'] for i,item in selected.items());return
 target.mkdir();byhash={value['SHA256']:i for i,value in selected.items()};assert len(byhash)==3456
 with tarfile.open(archive,'r|gz') as t:
  for item in t:
   if not item.name.startswith('objects/'):continue
   digest=item.name.split('/')[1]
   if digest not in byhash:continue
   assert item.size<1048576;raw=t.extractfile(item).read();assert hashlib.sha256(raw).hexdigest()==digest
   (target/f'run-{byhash[digest]:04d}.json').write_bytes(raw)
 assert len(list(target.iterdir()))==3456

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source-study',type=Path,required=True);parser.add_argument('--evidence-repo',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);parser.add_argument('--phase',choices=['ingest','verify'],required=True);a=parser.parse_args()
 source=a.source_study.resolve();repository=a.evidence_repo.resolve();work=a.work.resolve();phase=work/a.phase;phase.mkdir(exist_ok=False)
 state=phase/'execution.json';record=dict(status='binding-inputs',phase=a.phase,pid=os.getpid(),startedUnix=time.time(),predictionFitOrScoring=False);write(state,record)
 binding=bind_source(source,repository);references(repository,work)
 axis_run=json.loads((work/'axis-execution.json').read_text());assert axis_run['status']=='passed-complete-cell-axis'
 binary=work/'build/numivivo-omics';old_binary=source/'build-pooled/numivivo-omics';assert sha(binary)==axis_run['binarySHA256'] and sha(old_binary)==OLD_BINARY_SHA
 freeze=json.loads((source/'donor-recipe-v2/source-freeze.json').read_text())
 adapter_path=source/'donor-recipe-v2/run_counts.py';assert sha(adapter_path)==next(x['SHA256'] for x in freeze['sourceFiles'] if x['path']=='run_counts.py')
 os.environ['NUMIVIVO_PARSE_STUDY']=str(source)
 spec=importlib.util.spec_from_file_location('frozen_parse_count_adapter',adapter_path);adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter);adapter.VERIFY=False
 runs=json.loads((source/'prepared/runs.json').read_text());groups=sorted({r['sample'] for r in runs});group_index={name:i for i,name in enumerate(groups)}
 assert len(groups)==24 and len(runs)==3456 and RECORDS*((1<<32)-1)<2**64
 matrix=np.zeros((24,FEATURES),dtype=np.uint64);digest=hashlib.sha256();received=0
 if a.phase=='ingest':commands={
  'file':[str(binary),'singlecell-file-count-stream','--axis',str(work/'axis'),'--output',str(work/'native-file')],
  'resident':[str(old_binary),'singlecell-count-stream-pseudobulk','--plan',str(source/'prepared/plan.json'),'--output',str(work/'native-resident')]}
 else:commands={
  'file':[str(binary),'singlecell-file-count-stream-verify',str(work/'native-file')],
  'resident':[str(old_binary),'singlecell-count-stream-verify',str(work/'native-resident')]}
 children={};handles=[];record.update(status='running',commands=commands,binaries=dict(file=sha(binary),resident=OLD_BINARY_SHA),sourceBindings=binding,platform=platform.platform());write(state,record)
 try:
  for name,command in commands.items():
   output=(phase/(name+'.stdout')).open('xb');error=(phase/(name+'.stderr')).open('xb');handles.extend([output,error]);process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=output,stderr=error)
   children[name]=dict(process=process,bytesSent=0,started=time.monotonic(),pid=process.pid)
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
   pending={};submitted=0
   for wanted in range(3456):
    while submitted<min(3456,wanted+16):pending[submitted]=pool.submit(adapter.one,runs[submitted]);submitted+=1
    run,data,bulk,totals,qc=pending.pop(wanted).result();assert run['run']==wanted
    before=json.loads((work/f'expected-ranges/run-{wanted:04d}.json').read_text())
    fields=['nonzeros','totalCounts','sourceGeneCountMismatch','sourceTscpCountMismatch','sourceMinusMatrixTotal','sourceMinusMatrixGeneCount','sortedRows']
    assert all(qc[k]==before[k] for k in fields)
    assert [(r['start'],r['stop'],r['SHA256']) for r in qc['ranges']]==[(r['start'],r['stop'],r['SHA256']) for r in before['ranges']]
    if a.phase=='verify':
     first=json.loads((work/f'ingest/run-{wanted:04d}.json').read_text())
     assert [(r['start'],r['stop'],r['SHA256']) for r in qc['ranges']]==[(r['start'],r['stop'],r['SHA256']) for r in first['ranges']]
    for child in children.values():
     n=child['process'].stdin.write(data);assert n==len(data);child['bytesSent']+=n
    digest.update(data);received+=len(data);matrix[group_index[run['sample']]]+=bulk
    write(phase/f'run-{wanted:04d}.json',qc)
    if wanted%25==0:print('PAIRED-STREAM',a.phase,wanted,run['hi'],received,flush=True)
  assert received==RECORDS*16
  for child in children.values():child['process'].stdin.close()
  results={}
  for name,child in children.items():
   process=child['process'];pid,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status)
   results[name]=dict(nativeExit=process.returncode,maximumNativeRSS=usage.ru_maxrss,seconds=time.monotonic()-child['started'],bytesSent=child['bytesSent'],pid=pid)
  for handle in handles:handle.close()
  record.update(native=results,streamSHA256=digest.hexdigest(),streamBytes=received,status='checking-complete-results');write(state,record)
  assert all(r['nativeExit']==0 for r in results.values()),results
  for name in children:
   receipt=json.loads((phase/(name+'.stdout')).read_text());assert receipt['streamBytes']==received and bytes(receipt['stream']['bytes']).hex()==digest.hexdigest()
  if a.phase=='verify':assert digest.hexdigest()==json.loads((work/'ingest/execution.json').read_text())['streamSHA256']
  independent=phase/'independent-counts.npz';np.savez_compressed(independent,groups=np.array(groups),counts=matrix)
  checked=compare(source,work,independent);write(phase/'independent-check.json',checked)
  assert bind_source(source,repository)==binding
  assert sha(binary)==axis_run['binarySHA256'] and sha(old_binary)==OLD_BINARY_SHA
  record.update(status='passed',cells=CELLS,records=RECORDS,independentCheck=checked,finishedUnix=time.time());write(state,record);print(json.dumps(record),flush=True)
 except BaseException as error:
  for child in children.values():
   process=child['process']
   if process.poll() is None:process.terminate();process.wait()
  for handle in handles:
   if not handle.closed:handle.close()
  record.update(status='failed',errorType=type(error).__name__,error=str(error),confirmedStreamBytes=received,finishedUnix=time.time());write(state,record);raise
if __name__=='__main__':main()
