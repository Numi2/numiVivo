from pathlib import Path
import subprocess,hashlib,json,struct
p=Path('/Users/n/numivivo-logcpm-feature-major-20260912');s=Path('/Users/n/numivivo-count-store-final-20260909/store');o=p/'cli-check';o.mkdir();binary=p/'cli-build/numivivo-omics';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();receipt=json.loads((s/'receipt.json').read_text());digest=bytes(receipt['counts']['bytes']).hex()
def cmd(plan,out,h):return [str(binary),'singlecell-logcpm-stream','--plan',str(plan),'--stream-sha256',h,'--output',str(out)]
with (s/'counts.bin').open('rb') as f:r=subprocess.run(['/usr/bin/time','-l',*cmd(p/'plan.json',o/'native.json',digest)],stdin=f,capture_output=True)
(o/'native.log').write_bytes(r.stdout+r.stderr);r.check_returncode();assert json.loads((o/'native.json').read_text())==json.loads((p/'native.json').read_text())
plan=dict(featureIDs=['a','b'],groupIDs=['x'],rowGroups=[0,0,0],rowTotals=[3,4,0],ordering='featureMajor');path=o/'plan.json';path.write_text(json.dumps(plan));data=b''.join(struct.pack('<IIQ',*x) for x in [(0,0,1),(1,0,4),(0,1,2)])
r=subprocess.run(cmd(path,o/'small.json',hashlib.sha256(data).hexdigest()),input=data,capture_output=True);r.check_returncode()
checks={}
for name,ordering in [('wrong-layout','rowMajor'),('unknown-layout','unknown')]:
 plan['ordering']=ordering;path=o/(name+'-plan.json');path.write_text(json.dumps(plan));out=o/(name+'.json');r=subprocess.run(cmd(path,out,hashlib.sha256(data).hexdigest()),input=data,capture_output=True);assert r.returncode!=0 and not out.exists();checks[name]=r.returncode
result=dict(status='PASS-feature-major-cli',fullNormanExact=True,layoutFailures=checks,binarySHA256=sha(binary),sourceSHA256=sha(p/'source/Sources/NumiVivoKit/Omics/VivoStreamedLogCPM.swift'),cliSHA256=sha(p/'source/Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift'),planSHA256=sha(p/'plan.json'),inputSHA256=digest,outputSHA256=sha(o/'native.json'),scriptSHA256=sha(Path(__file__)))
(o/'receipt.json').write_text(json.dumps(result,indent=2));print(result)
