"""Measure actual product publish/reconstruction on the complete qualified Parse source."""
import argparse,hashlib,json,os,subprocess,tarfile,time
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(p,d):
 t=p.with_name(p.name+'.next');t.write_text(json.dumps(d,indent=2)+'\n');os.replace(t,p)
def measured(command,s,name):
 start=time.perf_counter()
 with (s/(name+'.stdout')).open('x') as out,(s/(name+'.stderr')).open('x') as err:
  process=subprocess.Popen(command,stdout=out,stderr=err);_,status,usage=os.wait4(process.pid,0);process.returncode=os.waitstatus_to_exitcode(status)
 result=dict(command=command,pid=process.pid,nativeExit=process.returncode,maximumNativeRSS=usage.ru_maxrss,seconds=time.perf_counter()-start)
 save(s/(name+'.json'),result);assert result['nativeExit']==0,result;return result

def compare_archive(s,repo,protocol):
 command=['git','-C',str(repo),'cat-file','blob',protocol['baselineCommit']+':'+protocol['baselineArchivePath']]
 process=subprocess.Popen(command,stdout=subprocess.PIPE)
 class Reader:
  def __init__(self):self.h=hashlib.sha256();self.size=0
  def read(self,n=-1):
   b=process.stdout.read(n);self.h.update(b);self.size+=len(b);return b
 reader=Reader();compared=0;target='objects/'+protocol['expectedReportSHA256'];found=False
 with tarfile.open(fileobj=reader,mode='r|gz') as archive:
  for item in archive:
   if item.name!=target:continue
   assert not found and item.isfile() and item.size==protocol['expectedReportBytes'];found=True
   with archive.extractfile(item) as old,(s/'native-file/report.json').open('rb') as new:
    for b in iter(lambda:old.read(1048576),b''):assert new.read(len(b))==b;compared+=len(b)
    assert not new.read(1)
 while reader.read(1048576):pass
 assert process.wait()==0 and found and reader.h.hexdigest()==protocol['baselineArchiveSHA256'] and compared==protocol['expectedReportBytes']
 return dict(allOriginalNativeReportBytesComparedDirectly=True,reportBytes=compared,reportSHA256=sha(s/'native-file/report.json'),baselineArchiveSHA256=reader.h.hexdigest())

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--repo',type=Path,required=True);a=p.parse_args();s=a.study.resolve();r=a.repo.resolve();protocol=json.loads((s/'protocol.json').read_text());freeze=json.loads((s/'source-freeze.json').read_text())
 def check():
  for name,digest in freeze['files'].items():assert sha(Path(name))==digest,name
 check();state=s/'pipeline.json';assert not state.exists();record=dict(status='running',pid=os.getpid(),startedUnix=time.time(),predictionFitOrScoring=False,protocolSHA256=sha(s/'protocol.json'));save(state,record)
 try:
  binary=str(s/'build/numivivo-omics');output=s/'native-file'
  record['publish']=measured([binary,'singlecell-file-expression',protocol['sourceBundle'],'--plan',str(s/'plan.json'),'--output',str(output)],s,'publish')
  record['comparison']=compare_archive(s,r,protocol);assert record['comparison']['reportSHA256']==protocol['expectedReportSHA256'];record['status']='published-and-all-report-bytes-exact';save(state,record)
  record['verify']=measured([binary,'singlecell-file-expression-verify',str(output)],s,'verify');check()
  record['memoryGatePassed']=record['publish']['maximumNativeRSS']<=protocol['maximumPublishRSS']
  record['elapsedGatePassed']=record['publish']['seconds']<=protocol['maximumPublishSeconds']
  record['status']='passed-complete-native-reconstruction-and-comparison';record['finishedUnix']=time.time();save(state,record);print(json.dumps(record))
 except BaseException as error:
  record.update(status='failed',errorType=type(error).__name__,error=str(error),finishedUnix=time.time());save(state,record);raise
if __name__=='__main__':main()
