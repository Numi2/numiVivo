from pathlib import Path
import json,hashlib,tarfile
p=Path(__file__).parent;m=json.loads((p/'manifest.json').read_text())
assert hashlib.sha256((p/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(p/'evidence.tar.gz','r:gz') as t:
 files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 assert set(files)==set(m['members'])
 for name,b in files.items():assert hashlib.sha256(b).hexdigest()==m['members'][name],name
 def read(n):return json.loads(files[n])
 protocol=read('protocol.json');freeze=read('prediction-freeze.json');scores=read('scores.json')
 assert hashlib.sha256(files['ContextKernel.swift']).hexdigest()==protocol['files']['ContextKernel.swift']
 assert files['ContextKernel.swift']==(p/'ContextKernel.swift').read_bytes()
 for name,h in protocol['inputs'].items():assert hashlib.sha256(files['inputs/'+name]).hexdigest()==h
 for study,f in freeze['folds'].items():assert hashlib.sha256(files['native/'+study+'.json']).hexdigest()==f['outputSHA256']
 assert read('score-verification.json')['donorBaselineScores']==225
 assert hashlib.sha256(files['scores.json']).hexdigest()==read('score-verification.json')['scoresSHA256']
 assert sum(x['predictions'] for x in read('verification.json')['folds'].values())==885000
 assert [k for k,v in scores['folds'].items() if v['fivePercentGainOverBoth']]==['HIRISA']
 assert read('summary.json')==json.loads((p/'summary.json').read_text())
 assert read('tests.json')['status']=='pass'
print('PASS complete archive, frozen input/output bindings and retained three-study verdict')
