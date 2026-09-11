"""Full Parse aggregate-interface comparison against the frozen prior native statistical owner."""
import argparse,concurrent.futures,gzip,hashlib,json,os,platform,subprocess,time
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def root_fields(raw):
 # Retain exact serialized values as well as decoded values, so all statistical
 # fields can be compared without weakening equality through float re-encoding.
 text=raw.decode();decoder=json.JSONDecoder();i=0;result={}
 while text[i].isspace():i+=1
 assert text[i]=='{';i+=1
 while True:
  while text[i].isspace():i+=1
  if text[i]=='}':return result
  key,end=decoder.raw_decode(text,i);i=end
  while text[i].isspace():i+=1
  assert text[i]==':';i+=1
  while text[i].isspace():i+=1
  start=i;value,i=decoder.raw_decode(text,i);assert key not in result;result[key]=(value,text[start:i])
  while text[i].isspace():i+=1
  if text[i]=='}':return result
  assert text[i]==',';i+=1

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();study=a.study.resolve()
 protocol=json.loads((study/'protocol.json').read_text());freeze=json.loads((study/'source-freeze.json').read_text());state=study/'native-pipeline.json';assert not state.exists()
 record=dict(status='starting',pid=os.getpid(),protocolSHA256=sha(study/'protocol.json'),startedUnix=time.time(),predictionFitOrScoring=False)
 write(state,record)
 for path,digest in freeze['files'].items():assert sha(Path(path))==digest,path
 counts=Path(protocol['countStudy']);source_state=json.loads((counts/'ingest/execution.json').read_text());assert source_state['status']=='passed' and source_state['records']==1373870697
 for owner,name in [('file','native-file'),('resident','native-resident')]:assert sha(counts/name/'receipt.json')==protocol['sourceReceiptSHA256'][owner]
 commands={
  'file':[str(study/'build/numivivo-omics'),'singlecell-file-expression',str(counts/'native-file'),'--plan',str(study/'file-plan.json'),'--output',str(study/'native-file')],
  'resident':[str(study/'native-reference'),str(counts/'native-resident'),str(study/'contrast.json'),'-']}
 def execute(owner):
  command=commands[owner];started=time.monotonic();result=dict(status='running',command=command,startedUnix=time.time(),platform=platform.platform());path=study/(owner+'-execution.json');write(path,result)
  with (study/(owner+'.stderr')).open('xb') as errors:
   if owner=='resident':
    process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=errors);result['pid']=process.pid;write(path,result)
    digest=hashlib.sha256();size=0
    with gzip.open(study/'resident-report.json.gz','xb',compresslevel=6) as target:
     for chunk in iter(lambda:process.stdout.read(1048576),b''):digest.update(chunk);size+=len(chunk);target.write(chunk)
    process.stdout.close();result.update(uncompressedReportSHA256=digest.hexdigest(),uncompressedReportBytes=size)
   else:
    with (study/'file.stdout').open('xb') as output:
     process=subprocess.Popen(command,stdout=output,stderr=errors);result['pid']=process.pid;write(path,result)
   pid,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status)
  result.update(status='passed' if process.returncode==0 else 'failed',nativeExit=process.returncode,maximumNativeRSS=usage.ru_maxrss,seconds=time.monotonic()-started,finishedUnix=time.time());write(path,result);return result
 try:
  record.update(status='native-comparison-running',commands=commands);write(state,record)
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   futures={owner:pool.submit(execute,owner) for owner in commands};results={owner:future.result() for owner,future in futures.items()}
  record.update(native=results,status='checking-native-results');write(state,record)
  assert all(result['nativeExit']==0 for result in results.values()),'One or both native owners failed; preserve both outcomes without changing the declared model'
  with (study/'file-verify.stdout').open('xb') as output,(study/'file-verify.stderr').open('xb') as errors:
   subprocess.run([str(study/'build/numivivo-omics'),'singlecell-file-expression-verify',str(study/'native-file')],stdout=output,stderr=errors,check=True)
  actual_raw=(study/'native-file/report.json').read_bytes();expected_raw=gzip.decompress((study/'resident-report.json.gz').read_bytes())
  actual=root_fields(actual_raw);expected=root_fields(expected_raw);assert actual.keys()==expected.keys()
  for key in actual:
   if key!='design':assert actual[key][1]==expected[key][1],key
  design=root_fields(actual['design'][1].encode());old_design=root_fields(expected['design'][1].encode());assert design.keys()==old_design.keys()
  for key in design:
   if key!='observations':assert design[key][1]==old_design[key][1],key
  observations=design['observations'][0];old=old_design['observations'][0];assert len(observations)==len(old)==24
  for i,(row,previous) in enumerate(zip(observations,old)):
   reference=row['fileMembership'];assert 'sourceCellIndices' not in row and 'fileMembership' not in previous
   assert reference['sourceCellCount']==len(previous['sourceCellIndices']) and reference['groupIndex']==design['sourcePseudobulkIndices'][0][i]
   assert bytes(reference['countBundleReceipt']['bytes']).hex()==protocol['sourceReceiptSHA256']['file']
   assert {k:v for k,v in row.items() if k!='fileMembership'}=={k:v for k,v in previous.items() if k!='sourceCellIndices'}
  checked=dict(status='passed-full-parse-native-aggregate-handoff',cells=725031,features=40352,observations=24,allNonDesignJSONValuesByteExact=True,allNonMembershipDesignJSONValuesByteExact=True,allOriginalObservationIdentitiesAndCellCountsExact=True,allFileMembershipReferencesBoundToInputReceipt=True,nativeReconstructionPassed=True,testedFeatures=actual['testedFeatures'][0],native=results,sourceReceiptSHA256=protocol['sourceReceiptSHA256'],
   memoryScope='Each native process including input loading, inference and output writing; baseline stdout is compressed by a separate Python adapter. Prior native core runner versus new product command, not an end-to-end product speed comparison.',predictionFitOrScoring=False)
  write(study/'native-comparison.json',checked)
  for path,digest in freeze['files'].items():assert sha(Path(path))==digest,path
  record.update(status=checked['status'],comparison=checked,finishedUnix=time.time());write(state,record);print(json.dumps(checked),flush=True)
 except BaseException as error:
  record.update(status='failed',errorType=type(error).__name__,error=str(error),finishedUnix=time.time());write(state,record);raise

if __name__=='__main__':main()
