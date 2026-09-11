"""Retain every input case/result, recipe and failure log in a restorable archive."""
from pathlib import Path
import gzip,hashlib,io,json,sys,tarfile
study=Path(sys.argv[1]);repo=Path(sys.argv[2]);out=Path(sys.argv[3]);out.mkdir(parents=True,exist_ok=False)
def sha(b):return hashlib.sha256(b).hexdigest()
items={}
for p in sorted(study.iterdir()):
 if p.is_file() and p.suffix in ['.json','.jsonl','.log','.stdout','.stderr','.py']:items['study/'+p.name]=p.read_bytes()
for p in sorted((study/'inputs').iterdir()):items['study/inputs/'+p.name]=p.read_bytes()
for name in ['sources.sha256']:
 items['study/build/'+name]=(study/'build'/name).read_bytes()
for name in ['VivoCountObservation.swift']:
 items['source/'+name]=(repo/'Sources/NumiVivoKit/Omics'/name).read_bytes()
for name in ['prepare.py','verify.py','Check.swift','test.sh','retain.py']:
 items['recipes/'+name]=(repo/'Tools/Omics/CountObservation'/name).read_bytes()
items['source/SingleCellCountObservationTests.swift']=(repo/'Tests/NumiVivoIntegrationTests/SingleCellCountObservationTests.swift').read_bytes()
manifest={'schemaVersion':1,'purpose':'Complete conditional count posterior numerical qualification; no biological calibration claim','files':{name:{'bytes':len(raw),'SHA256':sha(raw)} for name,raw in items.items()}}
archive=out/'results.tar.gz'
with archive.open('wb') as f,gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as gz,tarfile.open(fileobj=gz,mode='w') as tar:
 for name,raw in sorted(items.items()):
  info=tarfile.TarInfo(name);info.size=len(raw);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(raw))
manifest['archive']={'bytes':archive.stat().st_size,'SHA256':sha(archive.read_bytes())}
# Read every archived byte back; compact retention is not a summary-only receipt.
with tarfile.open(archive,'r:gz') as tar:
 assert set(tar.getnames())==set(items)
 for name in items:assert tar.extractfile(name).read()==items[name]
(out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
(out/'results.json').write_bytes((study/'results.json').read_bytes())
print(json.dumps({'files':len(items),'expandedBytes':sum(map(len,items.values())),'archive':manifest['archive'],'allArchivedBytesVerified':True}))
