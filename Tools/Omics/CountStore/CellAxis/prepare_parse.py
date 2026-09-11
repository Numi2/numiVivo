"""Stream every qualified Parse cell identity and exact matrix total into the native file axis."""
from pathlib import Path
import argparse,hashlib,json,os,platform,subprocess,time
from parse_support import *

def main():
 p=argparse.ArgumentParser();p.add_argument('--source-study',type=Path,required=True);p.add_argument('--evidence-repo',type=Path,required=True);p.add_argument('--work',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args()
 source=a.source_study.resolve();work=a.work.resolve();binary=a.binary.resolve();work.mkdir(parents=True,exist_ok=True)
 state=work/'axis-execution.json';assert not state.exists()
 binding=bind_source(source,a.evidence_repo.resolve());write(work/'source-bindings.json',binding)
 header=dict(schemaVersion=1,metadata=metadata_without_cells(source/'prepared/plan.json'),cellCount=CELLS,nonzeros=RECORDS,hasRowTotals=True,annotations=[],
  sourceDeclaration=json.dumps(dict(sourcePlanSHA256=PLAN_SHA,preparationArchiveSHA256=PREPARATION_SHA,countArchiveSHA256=COUNTS_SHA,retainedMatrixTotalsSHA256=binding['independentTotalsSHA256'],scope='All original selected count-axis rows; exact prior matrix totals, not historical QC. No additional filtering or prediction.',license='Parse Biosciences CC BY-NC 4.0'),sort_keys=True))
 write(work/'header.json',header);command=[str(binary),'singlecell-cell-axis-import','--header',str(work/'header.json'),'--output',str(work/'axis')]
 record=dict(status='running',command=command,binarySHA256=sha(binary),platform=platform.platform(),startedUnix=time.time(),sourceBindings=binding,predictionFitOrScoring=False);write(state,record)
 with (work/'axis.stdout').open('xb') as output,(work/'axis.stderr').open('xb') as errors:
  process=subprocess.Popen(command,stdin=subprocess.PIPE,stdout=output,stderr=errors);digest=hashlib.sha256();count=0;started=time.monotonic()
  try:
   for row in source_rows(source):
    raw=(json.dumps(row,separators=(',',':'),ensure_ascii=False)+'\n').encode();process.stdin.write(raw);digest.update(raw);count+=1
   process.stdin.close();pid,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status)
   assert process.returncode==0,(process.returncode,(work/'axis.stderr').read_text())
   receipt=json.loads((work/'axis.stdout').read_text());assert bytes(receipt['inputJSONL']['bytes']).hex()==digest.hexdigest()
   record.update(status='native-import-passed',cells=count,inputJSONLSHA256=digest.hexdigest(),nativeExit=process.returncode,maximumNativeRSS=usage.ru_maxrss,seconds=time.monotonic()-started);write(state,record)
  except BaseException as error:
   if process.poll() is None:process.terminate();process.wait()
   record.update(status='failed',errorType=type(error).__name__,error=str(error),cellsWritten=count);write(state,record);raise
 check=check_axis(source,work/'axis',header);write(work/'axis-independent-check.json',check)
 with (work/'axis-verify.stdout').open('xb') as output,(work/'axis-verify.stderr').open('xb') as errors:
  subprocess.run([str(binary),'singlecell-cell-axis-verify',str(work/'axis')],stdout=output,stderr=errors,check=True)
 assert (work/'axis.stdout').read_bytes()==(work/'axis-verify.stdout').read_bytes()
 record.update(status='passed-complete-cell-axis',independentCheck=check,nativeReopenVerified=True,finishedUnix=time.time());write(state,record);print(json.dumps(record),flush=True)
if __name__=='__main__':main()
