"""Run a pinned native paired analysis, streaming input/output across SSH."""
from pathlib import Path
import gzip,hashlib,json,os,shlex,subprocess,sys,threading,time
root=Path(sys.argv[1]);origin=sys.argv[2];remote='/Users/n/numivivo-joint-counts-20260911'
freeze=json.loads((root/'runtime-freeze.json').read_bytes());expected=freeze['jointBinarySHA256']
native=remote+'/build/joint-counts'
# Verify the exact binary before execution, without changing it during a run.
script='import hashlib,pathlib,sys; p=pathlib.Path('+repr(native)+'); assert hashlib.sha256(p.read_bytes()).hexdigest()=='+repr(expected)
command='python3 -c '+shlex.quote(script)+' && /usr/bin/time -l '+shlex.quote(native)
ssh=['ssh','-4','-C','-o','ConnectTimeout=10','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=1','macmini',command]
started=time.time();errors=[];record={'status':'running','startedUnix':started,'origin':origin,'controllerPID':os.getpid(),'binarySHA256':expected,'command':ssh}
def save(): (root/(origin+'-pipeline.json')).write_text(json.dumps(record,indent=2)+'\n')
save()
with (root/(origin+'-native.stderr')).open('xb') as stderr:
 p=subprocess.Popen(ssh,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=stderr)
 record['sshPID']=p.pid;save()
 def send():
  try:
   h=hashlib.sha256()
   with gzip.open(root/(origin+'-input.json.gz'),'rb') as f:
    while data:=f.read(1<<20):h.update(data);p.stdin.write(data)
   expected_input=json.loads((root/(origin+'-prepare-state.json')).read_bytes())['inputSHA256']
   assert h.hexdigest()==expected_input
  except BaseException as e:errors.append(repr(e))
  finally:p.stdin.close()
 writer=threading.Thread(target=send);writer.start();output_hash=hashlib.sha256();size=0
 with (root/(origin+'-native.jsonl.gz')).open('xb') as raw,gzip.GzipFile(filename='',fileobj=raw,mode='wb',mtime=0) as f:
  while data:=p.stdout.read(1<<20):output_hash.update(data);size+=len(data);f.write(data)
 code=p.wait();writer.join()
record.update(status='passed' if code==0 and not errors else 'failed',exitCode=code,errors=errors,seconds=time.time()-started,outputSHA256=output_hash.hexdigest(),outputBytes=size,finishedUnix=time.time());save();print(json.dumps(record));sys.exit(code or bool(errors))
