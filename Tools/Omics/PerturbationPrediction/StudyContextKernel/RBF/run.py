from pathlib import Path
import json,hashlib,subprocess,time
import numpy as np
s=Path(__file__).parent;source=s.parent/'numivivo-parse-context-evaluation-20260912';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
meta=json.loads((source/'source-metadata.json').read_text());a=np.load(source/'source-cohorts.npz');old=json.loads((source/'protocol.json').read_text())
for n in ['source-cohorts.npz','source-metadata.json']:
 expected=[v for k,v in old['dependencies'].items() if Path(k).name==n];assert expected==[sha(source/n)]
(s/'inputs').mkdir();(s/'native').mkdir();bindings={}
for query in sorted(meta['studies']):
 train=[v for v in sorted(meta['studies']) if v!=query]
 obj=dict(featureIDs=meta['featureIDs'],queryStudy=query,trainingDonorIDs=[d for t in train for d in meta['studies'][t]],trainingStudies=[t for t in train for d in meta['studies'][t]],queryDonorIDs=meta['studies'][query],queryControls=a[query+'_controls'].tolist())
 for key in ['controls','treated']:obj[key]=np.concatenate([a[t+'_'+key] for t in train]).tolist()
 p=s/'inputs'/(query+'.json');p.write_text(json.dumps(obj,separators=(',',':')));bindings[p.name]=sha(p)
execution=dict(inputs=bindings,protocolSHA256=sha(s/'protocol.json'),sourceSHA256=sha(s/'ContextKernel.swift'),binarySHA256=sha(s/'context-kernel'),driverSHA256=sha(Path(__file__)))
(s/'execution.json').write_text(json.dumps(execution,indent=2)+'\n');records={}
for n in bindings:
 start=time.time()
 with (s/'native'/(n+'.log')).open('xb') as f:subprocess.run([str(s/'context-kernel'),str(s/'inputs'/n),str(s/'native'/n)],stdout=f,stderr=f,check=True)
 records[n]=dict(SHA256=sha(s/'native'/n),seconds=time.time()-start);print(n,'complete',flush=True)
(s/'prediction-freeze.json').write_text(json.dumps(dict(createdUnix=time.time(),outputs=records,executionSHA256=sha(s/'execution.json'),scored=False),indent=2)+'\n')
