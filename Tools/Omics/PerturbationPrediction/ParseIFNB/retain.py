"""Retain exact admitted axes/count results; reuse the verified archive format."""
import argparse,gzip,hashlib,io,json,subprocess,tarfile
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--preparation-only',action='store_true');p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();repo=Path(__file__).resolve().parents[4]
 if a.preparation_only:
  assert json.loads((a.study/'sources/chunk-index-verification.json').read_text())['status']=='passed'
  assert json.loads((a.study/'prepared/freeze.json').read_text())['status']=='count-admission-only'
 else:
  for mode in ['ingest','replay']:assert json.loads((a.study/mode/'execution.json').read_text())['status']=='passed'
  assert json.loads((a.study/'counts-verification.json').read_text())['completeNativeReplay']=='passed'
 a.out.mkdir(parents=True,exist_ok=False);files=[]
 for name in (['sources','axes','chunks','prepared'] if a.preparation_only else ['sources','axes','chunks','prepared','native','ingest','replay']):
  files += [x for x in (a.study/name).rglob('*') if x.is_file()]
 files += [x for x in a.study.iterdir() if x.is_file() and x.suffix in ['.json','.log','.py'] and (not a.preparation_only or x.name not in ['counts.log','replay.log'])]
 files += [a.study/'build/numivivo-omics',a.study/'build/sources.sha256'];files=sorted(set(files));assert all(not f.is_symlink() for f in files)
 members={str(f.relative_to(a.study)):dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode & 0o777) for f in files};objects={}
 for f in files:objects.setdefault(members[str(f.relative_to(a.study))]['SHA256'],f)
 archive=a.out/('preparation.tar.gz' if a.preparation_only else 'results.tar.gz');raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  item=tarfile.TarInfo('members.json');item.size=len(raw);t.addfile(item,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   item=tarfile.TarInfo('objects/'+h);item.size=f.stat().st_size
   with f.open('rb') as data:t.addfile(item,data)
 subprocess.run(['python3',str(repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py'),'verify','--archive',str(archive)],check=True)
 write(a.out/'contents.json',dict(status='verified',members=members,logicalBytes=sum(v['bytes'] for v in members.values()),uniqueBytes=sum(f.stat().st_size for f in objects.values()),archiveSHA256=sha(archive),scope=('PREPARATION ONLY: Complete selected cell axes, HDF5-enumerated chunk maps, source range identities, tested executable and failed/abandoned preparation attempts. Complete native count scanning, independent aggregate comparison and full replay remain pending.' if a.preparation_only else 'Complete selected cell axes, native reports and receipts, independent aggregates, every requested source range identity, exact executable and failed/abandoned attempts.')+' The 227 GB source and 21.98 GB canonical stream are not retained. Original data and derivatives: Parse Biosciences, CC BY-NC 4.0.'))
 recipefiles=[x for x in Path(__file__).parent.iterdir() if x.is_file()]
 recipefiles += [repo/x for x in ['Sources/NumiVivoKit/Omics/VivoCountStreamPseudobulk.swift','Sources/NumiVivoKit/Omics/VivoH5ADPseudobulk.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift','Tests/NumiVivoIntegrationTests/SingleCellCountStreamTests.swift','Tools/Omics/H5AD/build.sh','Tools/Omics/CountStore/Stream/test.sh','Tools/Omics/CountStore/Stream/README.md','Tools/Omics/PerturbationPrediction/Duration/archive.py','Tools/Omics/PerturbationPrediction/Duration/common.py']]
 sourcefiles=[dict(path=str(f.relative_to(repo)),bytes=f.stat().st_size,SHA256=sha(f),rawUTF8=f.read_text()) for f in sorted(set(recipefiles))]
 (a.out/'recipe.json.gz').write_bytes(gzip.compress((json.dumps(dict(schemaVersion=1,sourceFiles=sourcefiles),sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
 records=[]
 for f in sorted(a.out.iterdir()):
  raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records));print(json.dumps(dict(files=len(files),archiveBytes=archive.stat().st_size,archiveSHA256=sha(archive))))
if __name__=='__main__':main()
