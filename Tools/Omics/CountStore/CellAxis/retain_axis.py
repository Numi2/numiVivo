"""Retain the dated complete Parse cell-axis qualification, excluding active count jobs."""
import argparse,gzip,hashlib,io,json,os,subprocess,sys,tarfile
from pathlib import Path
from parse_support import sha,write,PREPARATION_SHA,COUNTS_SHA

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--dependencies-repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 study=a.study.resolve();repo=Path(__file__).resolve().parents[4];dependencies_repo=a.dependencies_repo.resolve()
 execution=json.loads((study/'axis-execution.json').read_text());freeze=json.loads((study/'count-recipe/source-freeze.json').read_text())
 assert execution['status']=='passed-complete-cell-axis' and execution['cells']==725031 and execution['nativeReopenVerified']
 assert execution['independentCheck']['allOriginalIdentityBytesAndRowsExact'] and execution['independentCheck']['allDeclaredMatrixTotalsAndCardinalitiesExact']
 assert freeze['unitTests']==14 and freeze['unitSuites']==2 and freeze['completeCellAxis']==execution
 assert sha(study/'build/numivivo-omics')==freeze['binarySHA256']==execution['binarySHA256']
 assert '14 tests in 2 suites passed' in (study/'tests-unicode-final.log').read_text()
 assert len(freeze['nativeBuildInputs'])==105
 for item in freeze['nativeBuildInputs']:assert sha(repo/item['path'])==item['SHA256'],item['path']
 for item in freeze['recipe']:assert sha(study/'count-recipe'/item['path'])==item['SHA256'],item['path']
 assert (study/'axis.stdout').read_bytes()==(study/'axis-verify.stdout').read_bytes()
 dependencies=[]
 for name,filename,expected in [('2026-09-11-preparation','preparation.tar.gz',PREPARATION_SHA),('2026-09-11-counts','results.tar.gz',COUNTS_SHA)]:
  relative='Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/'+name
  assert sha(dependencies_repo/relative/filename)==expected
  dependencies.append(dict(path=relative,archive=filename,archiveSHA256=expected))
 names=['header.json','source-bindings.json','axis-execution.json','axis-independent-check.json','axis.stdout','axis.stderr','axis-verify.stdout','axis-verify.stderr','prepare-parse.log','build-unicode-final.log','tests-unicode-final.log','build/numivivo-omics','build/sources.sha256']
 for folder in ['axis','axis-recipe','count-recipe']:
  names.extend(str(f.relative_to(study)) for f in (study/folder).iterdir() if f.is_file())
 names=sorted(set(names));files={name:study/name for name in names}
 assert all(f.is_file() and not f.is_symlink() for f in files.values())
 members={name:dict(SHA256=sha(f),bytes=f.stat().st_size,mode=f.stat().st_mode & 0o777) for name,f in files.items()}
 objects={}
 for name,f in files.items():objects.setdefault(members[name]['SHA256'],f)
 assert not a.out.exists();temporary=a.out.with_name(a.out.name+'.incomplete');temporary.mkdir(parents=True,exist_ok=False)
 archive=temporary/'results.tar.gz';raw=(json.dumps(members,sort_keys=True,indent=2)+'\n').encode()
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:
  item=tarfile.TarInfo('members.json');item.size=len(raw);t.addfile(item,io.BytesIO(raw))
  for h,f in sorted(objects.items()):
   item=tarfile.TarInfo('objects/'+h);item.size=f.stat().st_size
   with f.open('rb') as source:t.addfile(item,source)
 verifier=dependencies_repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py'
 subprocess.run([sys.executable,str(verifier),'verify','--archive',str(archive)],check=True)
 write(temporary/'contents.json',dict(status='verified',members=members,archiveSHA256=sha(archive),logicalBytes=sum(v['bytes'] for v in members.values()),scope='Complete cell-axis import and native reopen only. Frozen full-count recipe retained, but no running count output is included. Parse data and derivatives: CC BY-NC 4.0.'))
 write(temporary/'summary.json',dict(schemaVersion=1,cellAxis=execution,unitTests=14,unitSuites=2,nativeBuildInputCount=105,fullPairedCountQualification='pending-separate-evidence',predictionFitOrScoring=False))
 write(temporary/'dependencies.json',dict(schemaVersion=1,required=dependencies,restore='Use Duration/archive.py restore for this archive. It restores the executed axis recipe and native executable. Reproducing source-independent checks also requires the exact preparation and count dependencies; use ParseIFNB/restore_donors.py to reconstruct those sources. No source count replay is implied by restoring bytes.'))
 source_names=set(item['path'] for item in freeze['nativeBuildInputs'])
 source_names.update(['Tests/NumiVivoIntegrationTests/SingleCellFileAxisTests.swift','Tools/Omics/H5AD/build.sh'])
 source_names.update(str(f.relative_to(repo)) for f in Path(__file__).parent.iterdir() if f.suffix in ('.py','.sh'))
 sources=[dict(path=name,bytes=(repo/name).stat().st_size,SHA256=sha(repo/name),rawUTF8=(repo/name).read_text()) for name in sorted(source_names)]
 recipe=dict(schemaVersion=1,baseCommit=freeze['baseCommit'],sourceFiles=sources,executedAxisRecipe='axis-recipe inside archive; this precedes the optional count-checker helper extension',pendingCountRecipe='count-recipe inside archive; qualification not yet claimed')
 (temporary/'recipe.json.gz').write_bytes(gzip.compress((json.dumps(recipe,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
 records=[]
 for f in sorted(temporary.iterdir()):
  raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(temporary/'manifest.json',dict(schemaVersion=1,records=records))
 subprocess.run([sys.executable,str(dependencies_repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(temporary)],check=True)
 assert all(sha(files[name])==item['SHA256'] for name,item in members.items()),'Artifact changed during packing'
 assert all(sha(repo/item['path'])==item['SHA256'] for item in sources),'Source changed during packing'
 os.rename(temporary,a.out)
 print(json.dumps(dict(status='retained-complete-cell-axis',members=len(members),archiveSHA256=sha(a.out/'results.tar.gz'),archiveBytes=(a.out/'results.tar.gz').stat().st_size,predictionFitOrScoring=False)))

if __name__=='__main__':main()
