"""Freeze every qualified training moment and its donor/condition identity."""
from pathlib import Path
import gzip, hashlib, json, sys

root=Path(sys.argv[1]);parent=Path('/Users/home/numivivo-count-calibration-20260911')
parent_manifest=json.loads((Path(__file__).resolve().parents[4]/'Tools/Omics/CountObservation/Calibration/evidence/2026-09-11/manifest.json').read_bytes())
def sha(data): return hashlib.sha256(data).hexdigest()
def write_gzip(path,data):
 with path.open('xb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(data)
freeze={'parentCommit':'b8cb39cc9a639fcb319322be0f14f4e109d65d6b','protocolSHA256':sha((root/'protocol.json').read_bytes()),'origins':{}}
for origin in ['Kang','HIRISA']:
 compressed=(parent/(origin+'-native.json.gz')).read_bytes()
 report=json.loads(gzip.decompress(compressed));metadata=(parent/'inputs'/(origin+'.json')).read_bytes();meta=json.loads(metadata)
 for relative,data in [(origin+'-native.json.gz',compressed),('inputs/'+origin+'.json',metadata)]:
  entry=parent_manifest['files']['study/'+relative];assert len(data)==entry['bytes'] and sha(data)==entry['SHA256']
 verification=json.loads((parent/(origin+'-verification.json')).read_bytes())
 assert verification['status']=='passed' and verification['sourceBinding']==report['bindingSHA256']
 assert report['metadataSHA256']==sha(metadata)
 assert len(meta['groups'])==len(report['moments'])
 assert sum(m['cells'] for m in report['moments'])==report['cells']
 for g,m in zip(meta['groups'],report['moments']):assert len(g['sourceRows'])==m['cells']
 assert all([f['featureID'] for f in model['features']]==meta['featureIDs'] for model in report['models'])
 binding={'parentReportCompressedSHA256':sha(compressed),'parentMetadataSHA256':sha(metadata),'parentRawCellStreamSHA256':report['streamSHA256'],'parentTrainingBindingSHA256':report['bindingSHA256']}
 obj={'origin':origin,'featureIDs':meta['featureIDs'],'controlConditionID':'control','treatedConditionID':'IFNB','sourceBindings':binding,
      'groups':[dict(donorID=g['donorID'],conditionID=g['conditionID'],moments=m) for g,m in zip(meta['groups'],report['moments'])]}
 raw=json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode();p=root/(origin+'-input.json.gz');write_gzip(p,raw)
 freeze['origins'][origin]=dict(binding,inputSHA256=sha(raw),compressedInputSHA256=sha(p.read_bytes()),features=len(meta['featureIDs']),cells=report['cells'],groups=len(meta['groups']))
(root/'input-freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
print(json.dumps(freeze,indent=2))
