"""Wait for an owned live count pipeline, then retain and restore its completed evidence."""
import argparse,gzip,hashlib,json,os,subprocess,sys,time
from pathlib import Path
from parse_support import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--source-study',type=Path,required=True);p.add_argument('--implementation-repo',type=Path,required=True);p.add_argument('--dependencies-repo',type=Path,required=True);p.add_argument('--evidence-out',type=Path,required=True);p.add_argument('--restoration-out',type=Path,required=True);a=p.parse_args()
 study=a.study.resolve();recipe=Path(__file__).resolve().parent;state=study/'count-retention-pipeline.json';assert not state.exists()
 freeze=json.loads((recipe/'source-freeze.json').read_text())
 for item in freeze['sourceFiles']:assert sha(recipe/item['name'])==item['SHA256']
 def read_pipeline():
  for attempt in range(10):
   try:return json.loads((study/'count-pipeline.json').read_text())
   except json.JSONDecodeError:
    if attempt==9:raise
    time.sleep(0.1)
 pipeline=read_pipeline();expected_pid=pipeline['pid']
 def owner():return subprocess.check_output(['ps','-p',str(expected_pid),'-o','lstart=','-o','command='],text=True).strip()
 identity=owner();assert str(study/'run_count_pipeline.py') in identity
 record=dict(status='waiting-for-complete-count-pipeline',pid=os.getpid(),sourcePipelinePID=expected_pid,sourcePipelineIdentity=identity,startedUnix=time.time(),recipe=str(recipe),predictionFitOrScoring=False)
 def save():
  temp=state.with_name(state.name+'.next');temp.write_text(json.dumps(record,indent=2)+'\n');os.replace(temp,state)
 save()
 try:
  while True:
   pipeline=read_pipeline();assert pipeline['pid']==expected_pid
   if pipeline['status']=='passed-paired-full-count-ingestion-and-replay':break
   assert pipeline['status'] in ('ingest-running','verify-running'),pipeline
   assert owner()==identity,'Source controller is absent or its identity changed; do not restart automatically'
   record.update(lastSourceStatus=pipeline['status'],lastLiveObservationUnix=time.time());save();time.sleep(30)
  record.update(status='retaining-completed-counts',sourcePipelineTerminal=pipeline);save()
  command=[sys.executable,str(recipe/'retain_counts.py'),'--study',str(study),'--source-study',str(a.source_study),'--implementation-repo',str(a.implementation_repo),'--dependencies-repo',str(a.dependencies_repo),'--out',str(a.evidence_out)]
  with (study/'retain-paired-counts.log').open('x') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  record.update(status='restoring-and-offline-checking');save()
  restored=study/'restored-paired-counts'
  command=[sys.executable,str(recipe/'restore_counts.py'),'--evidence',str(a.evidence_out),'--implementation-repo',str(a.implementation_repo),'--dependencies-repo',str(a.dependencies_repo),'--source-study',str(a.source_study),'--out',str(restored),'--reuse-from',str(study)]
  with (study/'restore-paired-counts.log').open('x') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
  restoration=json.loads((restored/'restoration.json').read_text());assert restoration['status']=='passed-count-restoration-and-offline-check'
  assert not a.restoration_out.exists();a.restoration_out.mkdir(parents=True,exist_ok=False)
  for name,path in [('restoration.json',restored/'restoration.json'),('offline-check.json',restored/'restoration-check.json'),('offline-check.log',restored/'restoration-check.log'),('retention.log',study/'retain-paired-counts.log'),('restoration.log',study/'restore-paired-counts.log')]:
   (a.restoration_out/name).write_bytes(path.read_bytes())
  sources=[]
  for item in freeze['sourceFiles']:
   f=recipe/item['name'];assert sha(f)==item['SHA256'];sources.append(dict(path=item['path'],bytes=f.stat().st_size,SHA256=item['SHA256'],rawUTF8=f.read_text()))
  raw=(json.dumps(dict(schemaVersion=1,sourceFiles=sources),sort_keys=True,separators=(',',':'))+'\n').encode();(a.restoration_out/'recipe.json.gz').write_bytes(gzip.compress(raw,mtime=0))
  records=[]
  for f in sorted(a.restoration_out.iterdir()):
   raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
   records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
  write(a.restoration_out/'manifest.json',dict(schemaVersion=1,records=records))
  subprocess.run([sys.executable,str(a.dependencies_repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(a.restoration_out)],check=True)
  record.update(status='retained-restored-and-offline-checked',finishedUnix=time.time(),countEvidence=str(a.evidence_out),restorationEvidence=str(a.restoration_out),restoration=restoration);save();print(json.dumps(record),flush=True)
 except BaseException as error:
  record.update(status='failed',errorType=type(error).__name__,error=str(error),finishedUnix=time.time());save();raise

if __name__=='__main__':main()
