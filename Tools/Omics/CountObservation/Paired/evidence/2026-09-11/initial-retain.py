"""Content-addressed, lossless retention with aliases for identical complete replays."""
from pathlib import Path,PurePosixPath
import gzip,hashlib,io,json,os,sys,tarfile

def digest(b):return hashlib.sha256(b).hexdigest()
def restore(evidence,destination):
 manifest=json.loads((evidence/'manifest.json').read_bytes());archive=evidence/'objects.tar.gz'
 assert digest(archive.read_bytes())==manifest['archive']['SHA256']
 destination.mkdir(parents=True,exist_ok=False)
 with tarfile.open(archive,'r:gz') as tar:
  assert set(tar.getnames())=={'objects/'+h for h in manifest['objects']}
  for name,entry in manifest['files'].items():
   path=PurePosixPath(name);assert not path.is_absolute() and '..' not in path.parts
   member=tar.getmember('objects/'+entry['SHA256']);assert member.isfile()
   raw=tar.extractfile(member).read();assert len(raw)==entry['bytes'] and digest(raw)==entry['SHA256']
   output=destination/str(path);output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(raw);output.chmod(entry['mode'])
 print(json.dumps({'restoredLogicalFiles':len(manifest['files']),'verifiedUniqueObjects':len(manifest['objects'])}))

if sys.argv[1]=='restore':restore(Path(sys.argv[2]),Path(sys.argv[3]));sys.exit()
study,repo,out=map(Path,sys.argv[1:4]);out.mkdir(parents=True,exist_ok=False)
files={};objects={}
def retain(name,p):
 raw=p.read_bytes();h=digest(raw);objects[h]=raw
 files[name]={'SHA256':h,'bytes':len(raw),'mode':0o755 if os.access(p,os.X_OK) else 0o644}
for p in sorted(study.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(study)
 if rel.parts[0] in ['publication','restored-evidence']:continue
 retain('study/'+str(rel),p)
for p in sorted((repo/'Tools/Omics/CountObservation/Paired').iterdir()):
 if p.is_file():retain('recipes/'+p.name,p)
for name in ['VivoCountObservation.swift','VivoCountObservationCalibration.swift','VivoPairedCountMoments.swift']:
 retain('source/'+name,repo/'Sources/NumiVivoKit/Omics'/name)
retain('source/SingleCellPairedCountTests.swift',repo/'Tests/NumiVivoIntegrationTests/SingleCellPairedCountTests.swift')
parent=Path('/Users/home/numivivo-count-calibration-20260911')
for origin in ['Kang','HIRISA']:
 retain('reference/'+origin+'-reference-moments.json.gz',parent/(origin+'-reference-moments.json.gz'))
 retain('reference/inputs/'+origin+'.json',parent/'inputs'/(origin+'.json'))
parent_manifest=repo/'Tools/Omics/CountObservation/Calibration/evidence/2026-09-11/manifest.json'
retain('reference/parent-manifest.json',parent_manifest)
external={'parentCalibrationArchive':{'commit':'b8cb39cc9a639fcb319322be0f14f4e109d65d6b','path':'Tools/Omics/CountObservation/Calibration/evidence/2026-09-11/objects.tar.gz','SHA256':json.loads(parent_manifest.read_bytes())['archive']['SHA256'],'purpose':'Original full cell-stream qualification and source-data lineage; not needed to recompute paired results from retained moments'}}
archive=out/'objects.tar.gz'
with archive.open('wb') as f,gzip.GzipFile(fileobj=f,filename='',mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w') as tar:
 for h,raw in sorted(objects.items()):
  info=tarfile.TarInfo('objects/'+h);info.size=len(raw);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(raw))
manifest={'schemaVersion':1,'format':'sha256-objects-with-logical-file-aliases','files':files,'objects':{h:{'bytes':len(b)} for h,b in objects.items()},'externalSources':external,'archive':{'bytes':archive.stat().st_size,'SHA256':digest(archive.read_bytes())},'scope':'Complete paired reports and all donor omissions, independent references, full source moments and identities, exact executable, source snapshots, logs and failures; original H5AD qualification remains bound to the published parent archive.'}
(out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n');(out/'results.json').write_bytes((study/'results.json').read_bytes())
print(json.dumps({'logicalFiles':len(files),'uniqueObjects':len(objects),'archive':manifest['archive']}))
