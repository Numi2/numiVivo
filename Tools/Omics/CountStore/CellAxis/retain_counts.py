"""Retain only a complete paired native count ingestion/replay qualification."""
import argparse,gzip,hashlib,io,json,os,subprocess,sys,tarfile
from pathlib import Path
from parse_support import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--source-study',type=Path,required=True);p.add_argument('--implementation-repo',type=Path,required=True);p.add_argument('--dependencies-repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 study=a.study.resolve();repo=a.implementation_repo.resolve();dependencies_repo=a.dependencies_repo.resolve()
 # Refuse incomplete work before creating any output or changing an active artifact.
 assert json.loads((study/'count-pipeline.json').read_text())['status']=='passed-paired-full-count-ingestion-and-replay','Both complete native phases are required'
 for phase in ['ingest','verify']:
  result=json.loads((study/phase/'execution.json').read_text());assert result['status']=='passed' and result['records']==1373870697
 assert not a.out.exists() and not a.out.with_name(a.out.name+'.incomplete').exists()
 freeze=json.loads((study/'count-recipe/source-freeze.json').read_text())
 for item in freeze['nativeBuildInputs']:assert sha(repo/item['path'])==item['SHA256'],item['path']
 assert sha(study/'build/numivivo-omics')==freeze['binarySHA256']
 recipe=study/'count-retention-recipe';assert recipe.is_dir()
 recipe_freeze=json.loads((recipe/'source-freeze.json').read_text())
 for item in recipe_freeze['sourceFiles']:assert sha(recipe/item['name'])==item['SHA256']==sha(repo/item['path']),item['path']
 assert sha(Path(__file__))==next(x['SHA256'] for x in recipe_freeze['sourceFiles'] if x['name']=='retain_counts.py')
 checked=study/'paired-counts-verification.json';log=study/'paired-counts-offline-check.log'
 with log.open('x') as output:
  subprocess.run([sys.executable,str(recipe/'check_counts.py'),'--source-study',str(a.source_study),'--work',str(study),'--dependencies-repo',str(dependencies_repo),'--out',str(checked)],stdout=output,stderr=subprocess.STDOUT,check=True)
 result=json.loads(checked.read_text());assert result['status']=='passed-complete-paired-count-evidence'
 axis_relative='Tools/Omics/CountStore/CellAxis/evidence/2026-09-11-axis'
 axis_archive=repo/axis_relative/'results.tar.gz';axis_sha='37ea8f865d6b3fea1c8785abcc61e419f0f2443983c07aa5895f8b58bb26643e';assert sha(axis_archive)==axis_sha
 dependencies=[dict(role='axis-and-native-executable',path=axis_relative,archive='results.tar.gz',archiveSHA256=axis_sha)]
 for role,name,filename,digest in [
  ('original-preparation','2026-09-11-preparation','preparation.tar.gz','f08c78ca7a74ab01eab6c3f7748615f45e217cbf95446ba35ef830eaf9c75d02'),
  ('original-complete-counts','2026-09-11-counts','results.tar.gz','ae94b3d586d2212fc09b1d451f893e7b3badc05ad57eb04550bece351af12b97')]:
  relative='Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/'+name;assert sha(dependencies_repo/relative/filename)==digest
  dependencies.append(dict(role=role,path=relative,archive=filename,archiveSHA256=digest))
 files=[]
 for folder in ['native-file','native-resident','ingest','verify','count-retention-recipe']:
  files.extend(f for f in (study/folder).rglob('*') if f.is_file() and '__pycache__' not in f.parts)
 files.extend(study/name for name in ['count-pipeline.json','run_count_pipeline.py','ingest.log','verify.log','paired-counts-verification.json','paired-counts-offline-check.log','count-retention-preflight.json'])
 files=sorted(set(files));assert all(not f.is_symlink() for f in files)
 members={str(f.relative_to(study)):dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode & 0o777) for f in files};objects={}
 for f in files:objects.setdefault(members[str(f.relative_to(study))]['SHA256'],f)
 temporary=a.out.with_name(a.out.name+'.incomplete');temporary.mkdir(parents=True,exist_ok=False);archive=temporary/'results.tar.gz'
 raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  item=tarfile.TarInfo('members.json');item.size=len(raw);t.addfile(item,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   item=tarfile.TarInfo('objects/'+h);item.size=f.stat().st_size
   with f.open('rb') as source:t.addfile(item,source)
 subprocess.run([sys.executable,str(dependencies_repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py'),'verify','--archive',str(archive)],check=True)
 write(temporary/'summary.json',result)
 write(temporary/'contents.json',dict(status='verified',members=members,archiveSHA256=sha(archive),logicalBytes=sum(v['bytes'] for v in members.values()),scope='Both complete same-cohort native count consumers, all3456range/QC records in each phase, exact original canonical input identities, all cell QC/membership and independent aggregates, per-child RSS, frozen checker and controller. No biological prediction. Parse derivatives: CC BY-NC 4.0.'))
 write(temporary/'dependencies.json',dict(schemaVersion=1,required=dependencies,restore='restore_counts.py combines this result with the exact axis/executable dependency in a fresh directory, then runs the frozen offline checker against a source-study bound to the original preparation and count archives. Source-study may be restored with ParseIFNB/restore_donors.py. Original-count dependencies transitively retain the previous native executable and source adapter. Restoration does not repeat native count execution or network replay.'))
 sources=[dict(path=item['path'],bytes=(recipe/item['name']).stat().st_size,SHA256=item['SHA256'],rawUTF8=(recipe/item['name']).read_text()) for item in recipe_freeze['sourceFiles']]
 (temporary/'recipe.json.gz').write_bytes(gzip.compress((json.dumps(dict(schemaVersion=1,sourceFiles=sources),sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
 records=[]
 for f in sorted(temporary.iterdir()):
  raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(temporary/'manifest.json',dict(schemaVersion=1,records=records))
 subprocess.run([sys.executable,str(dependencies_repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(temporary)],check=True)
 assert all(sha(f)==members[str(f.relative_to(study))]['SHA256'] for f in files),'Completed artifact changed during packing'
 for item in freeze['nativeBuildInputs']:assert sha(repo/item['path'])==item['SHA256'],item['path']
 assert sha(axis_archive)==axis_sha
 os.rename(temporary,a.out)
 print(json.dumps(dict(status='retained-complete-paired-count-evidence',members=len(members),archiveSHA256=sha(a.out/'results.tar.gz'),archiveBytes=(a.out/'results.tar.gz').stat().st_size,predictionFitOrScoring=False)))

if __name__=='__main__':main()
