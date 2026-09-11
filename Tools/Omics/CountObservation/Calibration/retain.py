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
for p in sorted((repo/'Tools/Omics/CountObservation/Calibration').iterdir()):
 if p.is_file():retain('recipes/'+p.name,p)
retain('recipes/CountObservationReference.py',repo/'Tools/Omics/CountObservation/verify.py')
for name in ['VivoCountObservation.swift','VivoCountObservationCalibration.swift']:
 retain('source/'+name,repo/'Sources/NumiVivoKit/Omics'/name)
retain('source/SingleCellCountCalibrationTests.swift',repo/'Tests/NumiVivoIntegrationTests/SingleCellCountCalibrationTests.swift')
external={}
for origin in ['Kang','HIRISA']:
 meta=json.loads((study/'inputs'/(origin+'.json')).read_bytes())
 external[meta['sourcePath']]={'bytes':meta['sourceBytes'],'SHA256':meta['sourceSHA256'],'purpose':'Original complete H5AD needed to reconstruct raw training cell streams'}
prior=Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs')
paths=[prior/x for x in ['feature-mapping.json','source-bindings.json','Kang-cohort.json','HIRISA-cohort.json','Kang-source-counts.npz','HIRISA-source-counts.npz']]
paths += [Path('/Users/home/numivivo-kang-gamma-product-20260909/report.json'),Path('/Users/home/numivivo-hiris-20260910/inference-inputs/cohorts/02.json.gz'),Path('/Users/home/numivivo-count-observation-20260911/inputs/cases.json'),Path('/Users/home/numivivo-count-observation-20260911/inputs/freeze.json')]
for p in paths:external[str(p)]={'bytes':p.stat().st_size,'SHA256':digest(p.read_bytes()),'purpose':'Parent source qualification or original control-cell case linkage'}
archive=out/'objects.tar.gz'
with archive.open('wb') as f,gzip.GzipFile(fileobj=f,filename='',mode='wb',mtime=0) as gz,tarfile.open(fileobj=gz,mode='w') as tar:
 for h,raw in sorted(objects.items()):
  info=tarfile.TarInfo('objects/'+h);info.size=len(raw);info.mode=0o644;info.mtime=0;tar.addfile(info,io.BytesIO(raw))
manifest={'schemaVersion':1,'format':'sha256-objects-with-logical-file-aliases','files':files,'objects':{h:{'bytes':len(b)} for h,b in objects.items()},'externalSources':external,'archive':{'bytes':archive.stat().st_size,'SHA256':digest(archive.read_bytes())},'scope':'Complete native training and replay reports, source moments, fitted models, query input/output and all failures; original full H5AD sources remain external.'}
(out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n');(out/'results.json').write_bytes((study/'results.json').read_bytes())
print(json.dumps({'logicalFiles':len(files),'uniqueObjects':len(objects),'archive':manifest['archive']}))
