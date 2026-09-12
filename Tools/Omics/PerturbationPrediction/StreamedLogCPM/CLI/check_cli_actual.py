from pathlib import Path
import hashlib,json,subprocess,struct
s=Path('/Users/n/numivivo-logcpm-consumer-20260912');base=Path('/Users/n/numivivo-streamed-logcpm-20260912');o=s/'cli-check-actual';o.mkdir()
binary=s/'cli-build/numivivo-omics';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
stream=base/'counts.bin';digest=sha(stream)
def command(plan,output,h):return [str(binary),'singlecell-logcpm-stream','--plan',str(plan),'--stream-sha256',h,'--output',str(output)]
with stream.open('rb') as f:r=subprocess.run(command(base/'plan.json',o/'real.json',digest),stdin=f,capture_output=True)
(o/'real.log').write_bytes(r.stdout+r.stderr);r.check_returncode();assert json.loads((o/'real.json').read_text())==json.loads((base/'native.json').read_text())
plan=o/'plan.json';plan.write_text(json.dumps(dict(featureIDs=['g'],groupIDs=['c'],rowGroups=[0,0],rowTotals=[2,0])))
valid=struct.pack('<IIQ',0,0,2);cases={'truncated':valid[:-1],'duplicate':valid+valid,'wrong-total':struct.pack('<IIQ',0,0,1),'out-of-axis':struct.pack('<IIQ',0,1,2),'empty':b''}
checks={}
for name,data in cases.items():
 out=o/(name+'.json');r=subprocess.run(command(plan,out,hashlib.sha256(data).hexdigest()),input=data,capture_output=True);assert r.returncode!=0 and not out.exists();checks[name]=r.returncode;(o/(name+'.log')).write_bytes(r.stdout+r.stderr)
out=o/'wrong-hash.json';r=subprocess.run(command(plan,out,'0'*64),input=valid,capture_output=True);assert r.returncode!=0 and not out.exists();checks['wrong-hash']=r.returncode
out=o/'valid.json';r=subprocess.run(command(plan,out,hashlib.sha256(valid).hexdigest()),input=valid,capture_output=True);r.check_returncode();v=json.loads(out.read_text());assert v['cellCounts']==[2] and v['zeroCellCounts']==[1]
before=out.read_bytes();r=subprocess.run(command(plan,out,hashlib.sha256(valid).hexdigest()),input=valid,capture_output=True);assert r.returncode!=0 and out.read_bytes()==before;checks['overwrite']=r.returncode
receipt=dict(status='PASS',realMatrixExact=True,failureCases=checks,zeroCellsRetained=True,binarySHA256=sha(binary),sourceSHA256=sha(s/'source/Sources/NumiVivoKit/Omics/VivoStreamedLogCPM.swift'),cliSHA256=sha(s/'source/Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift'),inputSHA256=digest,outputSHA256=sha(o/'real.json'),scriptSHA256=sha(Path(__file__)))
(o/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
