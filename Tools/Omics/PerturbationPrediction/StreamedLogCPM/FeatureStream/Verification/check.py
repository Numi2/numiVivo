from pathlib import Path
import subprocess,json,hashlib,shutil,struct
p=Path('/Users/n/numivivo-logcpm-feature-verify-20260912');source=Path('/Users/n/numivivo-count-store-final-20260909/store');prior=Path('/Users/n/numivivo-logcpm-feature-output-20260912/cli-check');o=p/'checks';o.mkdir();binary=p/'cli-build/numivivo-omics';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def command(bundle):return [str(binary),'singlecell-logcpm-feature-verify',str(bundle)]
with (source/'counts.bin').open('rb') as f:r=subprocess.run(['/usr/bin/time','-l',*command(prior/'norman')],stdin=f,capture_output=True)
(o/'norman.log').write_bytes(r.stdout+r.stderr);r.check_returncode()
good=struct.pack('<IIQ',0,1,1)+struct.pack('<IIQ',0,3,2)
r=subprocess.run(command(prior/'small'),input=good,capture_output=True);r.check_returncode();checks={}
for kind in ['mean','summary','truncated','trailing','symlink','input']:
 bundle=o/kind;shutil.copytree(prior/'small',bundle);receipt=json.loads((bundle/'receipt.json').read_text())
 if kind=='mean':
  f=bundle/'means.bin';data=bytearray(f.read_bytes());data[8:16]=struct.pack('<d',999.);f.write_bytes(data);receipt['meansSHA256']=sha(f)
 elif kind=='summary':
  f=bundle/'summary.json';v=json.loads(f.read_text());v['cellCounts']=[999];f.write_text(json.dumps(v));receipt['summarySHA256']=sha(f)
 elif kind=='truncated':
  f=bundle/'means.bin';f.write_bytes(f.read_bytes()[:-8]);receipt['meansSHA256']=sha(f);receipt['meansBytes']=str(f.stat().st_size)
 elif kind=='trailing':
  f=bundle/'means.bin';f.write_bytes(f.read_bytes()+b'12345678');receipt['meansSHA256']=sha(f);receipt['meansBytes']=str(f.stat().st_size)
 elif kind=='symlink':
  f=bundle/'means.bin';f.unlink();f.symlink_to(prior/'small/means.bin')
 (bundle/'receipt.json').write_text(json.dumps(receipt))
 r=subprocess.run(command(bundle),input=good[:-1] if kind=='input' else good,capture_output=True);assert r.returncode!=0,kind;checks[kind]=r.returncode;(o/(kind+'.log')).write_bytes(r.stdout+r.stderr)
r=dict(status='PASS',fullNormanReconstructed=True,smallReconstructed=True,rehashedTamperingRejected=checks,binarySHA256=sha(binary),cliSHA256=sha(p/'source/Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift'),scriptSHA256=sha(Path(__file__)))
(o/'verification.json').write_text(json.dumps(r,indent=2));print(r)
