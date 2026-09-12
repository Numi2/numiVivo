from pathlib import Path
import subprocess,json,hashlib,struct,resource,signal,numpy as np
p=Path('/Users/n/numivivo-logcpm-feature-output-20260912');source=Path('/Users/n/numivivo-count-store-final-20260909/store');old=Path('/Users/n/numivivo-logcpm-feature-major-20260912');o=p/'cli-check';o.mkdir();binary=p/'cli-build/numivivo-omics';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def cmd(plan,out,h):return [str(binary),'singlecell-logcpm-feature-stream','--plan',str(plan),'--stream-sha256',h,'--output',str(out)]
def verify(out):
 r=json.loads((out/'receipt.json').read_text());s=json.loads((out/'summary.json').read_text())
 assert r['format']=='feature-major-group-f64-le-v1'
 for n in ['means','summary','plan']:assert sha(out/(n+('.bin' if n=='means' else '.json')))==r[n+'SHA256']
 assert (out/'means.bin').stat().st_size==int(r['meansBytes'])==len(s['featureIDs'])*len(s['groupIDs'])*8
 return s
h=bytes(json.loads((source/'receipt.json').read_text())['counts']['bytes']).hex()
with (source/'counts.bin').open('rb') as f:r=subprocess.run(['/usr/bin/time','-l',*cmd(old/'plan.json',o/'norman',h)],stdin=f,capture_output=True)
(o/'norman.log').write_bytes(r.stdout+r.stderr);r.check_returncode();s=verify(o/'norman');assert sha(o/'norman/means.bin')==sha(p/'means.bin');assert s==json.loads((p/'summary.json').read_text())
plan=dict(featureIDs=['zero-before','a','zero-between','b','zero-after'],groupIDs=['x'],rowGroups=[0,0],rowTotals=[3,0],ordering='featureMajor');path=o/'plan.json';path.write_text(json.dumps(plan));good=struct.pack('<IIQ',0,1,1)+struct.pack('<IIQ',0,3,2)
checks={}
for name,data in [('duplicate',good[:16]*2),('decreasing',good[16:]+good[:16]),('missing-total',good[:16]),('partial',good[:-1]),('excess',struct.pack('<IIQ',0,1,4))]:
 out=o/name;r=subprocess.run(cmd(path,out,hashlib.sha256(data).hexdigest()),input=data,capture_output=True);assert r.returncode!=0 and not out.exists();checks[name]=r.returncode
out=o/'wrong-hash';r=subprocess.run(cmd(path,out,'0'*64),input=good,capture_output=True);assert r.returncode!=0 and not out.exists();checks['wrong-hash']=r.returncode
out=o/'small';r=subprocess.run(cmd(path,out,hashlib.sha256(good).hexdigest()),input=good,capture_output=True);r.check_returncode();s=verify(out);v=np.fromfile(out/'means.bin',dtype='<f8');assert s['cellCounts']==[2] and s['zeroCellCounts']==[1] and np.all(v[[0,2,4]]==0)
expected=np.log1p(np.array([0,1,0,2,0])/3*1e6)/2;np.testing.assert_allclose(v,expected,rtol=1e-15)
before=sha(out/'means.bin');r=subprocess.run(cmd(path,out,hashlib.sha256(good).hexdigest()),input=good,capture_output=True);assert r.returncode!=0 and sha(out/'means.bin')==before
checks['overwrite']=r.returncode
def limit():signal.signal(signal.SIGXFSZ,signal.SIG_IGN);resource.setrlimit(resource.RLIMIT_FSIZE,(1,1))
out=o/'write-failure';r=subprocess.run(cmd(path,out,hashlib.sha256(good).hexdigest()),input=good,capture_output=True,preexec_fn=limit);assert r.returncode!=0 and not out.exists();checks['write-failure']=r.returncode;assert not list(o.glob('.numivivo-logcpm-features-*'))
r=dict(status='PASS',fullNormanBinaryExact=True,metadataExact=True,allZeroFeaturesPreserved=True,zeroCellsPreserved=True,failureChecks=checks,stagingClean=True,binarySHA256=sha(binary),scriptSHA256=sha(Path(__file__)),sourceSHA256=sha(p/'source/Sources/NumiVivoKit/Omics/VivoStreamedLogCPM.swift'),cliSHA256=sha(p/'source/Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift'),bundleReceiptSHA256=sha(o/'norman/receipt.json'))
(o/'verification.json').write_text(json.dumps(r,indent=2));print(r)
