"""Retain completed donor ingestion/replay, referencing immutable preparation."""
import argparse,gzip,hashlib,io,json,os,subprocess,sys,tarfile
from pathlib import Path
from retain import sha,write
from verify_donor_preparation import verify as verify_preparation


def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();study=a.study.resolve();repo=Path(__file__).resolve().parents[4]
 # Fail before creating any result directory unless both complete phases exist.
 for phase in ['ingest','verify']:
  complete=json.loads((study/'donor-execution'/phase/'complete.json').read_text())
  assert complete['status']=='passed' and complete['cells']==725031 and complete['records']==1373870697 and len(complete['donors'])==12
 preparation=repo/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-preparation/preparation.tar.gz'
 bindings=verify_preparation(study,preparation);write(study/'preparation-bindings.json',bindings)
 recipe=study/'donor-recipe-v2';freeze=json.loads((recipe/'source-freeze.json').read_text())
 for item in freeze['sourceFiles']:assert sha(recipe/item['path'])==item['SHA256']
 assert sha(study/'build-pooled/numivivo-omics')==freeze['binarySHA256']
 assert freeze['binarySHA256']=='20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516','This dated retainer requires the qualified executable retained in its dependency archive'
 subprocess.run([sys.executable,str(recipe/'check_donors.py')],env=dict(os.environ,NUMIVIVO_PARSE_STUDY=str(study)),check=True)
 checked=json.loads((study/'donor-counts-verification.json').read_text());assert checked['status']=='passed' and checked['predictionFitOrScoring'] is False
 dependencies=[]
 for relative,filename,expected in [
  ('Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-preparation','preparation.tar.gz','f08c78ca7a74ab01eab6c3f7748615f45e217cbf95446ba35ef830eaf9c75d02'),
  ('Tools/Omics/CountStore/Stream/Memory/evidence/2026-09-11','results.tar.gz','fa98901c2f0ce867620f7cd37494f74fca4c238dd32159c4a7751e55cbe5abe4')]:
  directory=repo/relative;assert sha(directory/filename)==expected
  subprocess.run([sys.executable,str(repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(directory)],check=True)
  dependencies.append(dict(path=relative,manifestSHA256=sha(directory/'manifest.json'),archive=filename,archiveSHA256=expected))
 files=[]
 for folder in ['donor-plans','donor-execution','donor-recipe','donor-recipe-v2','ingest']:
  files += [f for f in (study/folder).rglob('*') if f.is_file() and '__pycache__' not in f.parts]
 for name in ['preparation-bindings.json','preparation-bindings-attempt1.json','preparation-binding-mutations.json','donor-complete-counts.npz','donor-counts-verification.json','donor-plan.log','donor-ingest.log','donor-ingest-v2.log','donor-verify.log','donor1-recovery.log','recover_donor1.py','counts.log','transport-tests.log','build-pooled/sources.sha256','pooled-source-build-verification.json','tests-pooled-final.log','host.json','adapter-environment.json']:
  f=study/name;assert f.is_file(),name;files.append(f)
 files=sorted(set(files));assert all(not f.is_symlink() for f in files)
 members={str(f.relative_to(study)):dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode & 0o777) for f in files};objects={}
 for f in files:objects.setdefault(members[str(f.relative_to(study))]['SHA256'],f)
 assert not a.out.exists();temporary=a.out.with_name(a.out.name+'.incomplete');temporary.mkdir(parents=True,exist_ok=False)
 archive=temporary/'results.tar.gz';raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  item=tarfile.TarInfo('members.json');item.size=len(raw);t.addfile(item,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   item=tarfile.TarInfo('objects/'+h);item.size=f.stat().st_size
   with f.open('rb') as source:t.addfile(item,source)
 subprocess.run([sys.executable,str(repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py'),'verify','--archive',str(archive)],check=True)
 write(temporary/'contents.json',dict(status='verified',members=members,archiveSHA256=sha(archive),logicalBytes=sum(v['bytes'] for v in members.values()),scope='Complete twelve-donor native ingestion and source replay, all independent cell and aggregate counts, original failed monolithic attempt and checker recovery. Preparation axes, chunk maps and executables are retained in the exact required dependency archives. No prediction fit or scoring. Parse Biosciences data and derivatives: CC BY-NC 4.0.'))
 write(temporary/'dependencies.json',dict(schemaVersion=1,required=dependencies,restore='Restore the preparation archive into a fresh study directory, then this archive into a separate directory and merge nonconflicting paths. The corrected executable is runtime/corrected/numivivo-omics in the memory-control dependency and must be placed at build-pooled/numivivo-omics; retain its exact hash. Read source-freeze for the executed recipe.'))
 checked['preparationBinding']={k:v for k,v in bindings.items() if k!='inputSHA256'}
 write(temporary/'summary.json',checked)
 sources=[Path(__file__).parent/'restore_donors.py',Path(__file__).parent/'check_preparation_mutations.py',Path(__file__).parent/'verify_donor_preparation.py',Path(__file__),Path(__file__).parent/'retain.py',repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py',repo/'Tools/Omics/PerturbationPrediction/Duration/common.py']
 recipe_data=dict(schemaVersion=1,sourceFiles=[dict(path=str(f.relative_to(repo)),bytes=f.stat().st_size,SHA256=sha(f),rawUTF8=f.read_text()) for f in sources]);(temporary/'recipe.json.gz').write_bytes(gzip.compress((json.dumps(recipe_data,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
 records=[]
 for f in sorted(temporary.iterdir()):
  raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(temporary/'manifest.json',dict(schemaVersion=1,records=records))
 subprocess.run([sys.executable,str(repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(temporary)],check=True)
 assert all(sha(f)==members[str(f.relative_to(study))]['SHA256'] for f in files),'Source artifact changed during packing'
 assert all(sha(study/name)==digest for name,digest in bindings['inputSHA256'].items()),'Frozen preparation input changed during packing'
 assert sha(preparation)==bindings['preparationArchiveSHA256']
 os.rename(temporary,a.out);print(json.dumps(dict(status='retained',files=len(files),archiveSHA256=sha(a.out/'results.tar.gz'),archiveBytes=(a.out/'results.tar.gz').stat().st_size)))

if __name__=='__main__':main()
