from pathlib import Path
import json,hashlib,numpy as np
s=Path('/Users/n/numivivo-streamed-logcpm-20260912')
a=json.loads((s/'native.json').read_text());p=json.loads((s/'plan.json').read_text())
assert a['featureIDs']==p['featureIDs'] and a['groupIDs']==p['groupIDs']
assert a['cellCounts']==[904,904,903] and a['zeroCellCounts']==[0,0,0]
ref=np.load(s/'reference.npy');error=float(np.max(abs(np.array(a['means'])-ref)))
np.testing.assert_allclose(a['means'],ref,rtol=0,atol=1e-11)
report={'status':'PASS','cells':2711,'features':36601,'nonzeros':5218473,'counts':11786194,'maximumAbsoluteMeanError':error,'tolerance':1e-11,'partitionPurpose':'deterministic numerical verification only; not biological groups','biologicalPrediction':False,'fullProductBuild':False,'nativeSourceSHA256':hashlib.sha256((s/'VivoStreamedLogCPM.swift').read_bytes()).hexdigest(),'sourceSHA256':'5fbff5a4d85e0df345f6502e966ec787a8a4c429fd6b88a8772c43fd915cf3ff','files':{n:hashlib.sha256((s/n).read_bytes()).hexdigest() for n in ['Main.swift','check','prepare.py','plan.json','counts.bin','native.json','reference.npy']}}
(s/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
