from pathlib import Path
import tarfile,json,hashlib
p=Path(__file__).parent;m=json.loads((p/'manifest.json').read_text())
assert hashlib.sha256((p/'evidence.tar.gz').read_bytes()).hexdigest()==m['archiveSHA256']
with tarfile.open(p/'evidence.tar.gz','r:gz') as t:
 files={x.name:t.extractfile(x).read() for x in t.getmembers() if x.isfile()}
 assert set(files)==set(m['members'])
 for n,b in files.items():assert hashlib.sha256(b).hexdigest()==m['members'][n],n
 a=json.loads(files['audit.json']);v=json.loads(files['verification.json']);f=json.loads(files['freeze.json'])
 assert a['bCellRows']==v['selectedBCells']==72446
 assert a['exactCorrespondences']==v['exactGenes']==11600 and a['missingFeatures']==v['missingGenes']==200
 assert hashlib.sha256(files['source-roster.json']).hexdigest()==f['sourceRosterSHA256']
 assert hashlib.sha256(files['source-roster-codes.npz']).hexdigest()==f['sourceCodesSHA256']
 assert hashlib.sha256(files['model-metadata.json']).hexdigest()==f['modelFeatureMetadataSHA256']
 assert hashlib.sha256(files['admit.py']).hexdigest()==f['scriptSHA256']
 assert hashlib.sha256(files['selected-labels.npz']).hexdigest()==a['rowCodesSHA256']
 assert not a['predictionFitted'] and not a['predictionScored'] and not a['modelPanelChanged']
print('PASS complete metadata archive and input bindings; prediction remains unexecuted')
