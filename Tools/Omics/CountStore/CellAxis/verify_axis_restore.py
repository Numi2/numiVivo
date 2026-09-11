"""Fresh extraction and native reopen of the exact dated cell-axis archive."""
import argparse,gzip,hashlib,json,subprocess,sys
from pathlib import Path
from parse_support import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--archive',type=Path,required=True);p.add_argument('--dependencies-repo',type=Path,required=True);p.add_argument('--restore-to',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 expected='37ea8f865d6b3fea1c8785abcc61e419f0f2443983c07aa5895f8b58bb26643e'
 assert sha(a.archive)==expected and not a.restore_to.exists() and not a.out.exists()
 verifier=a.dependencies_repo/'Tools/Omics/PerturbationPrediction/Duration/archive.py'
 a.out.mkdir(parents=True,exist_ok=False)
 with (a.out/'extraction.log').open('xb') as log:
  subprocess.run([sys.executable,str(verifier),'restore','--archive',str(a.archive),'--out',str(a.restore_to)],stdout=log,stderr=subprocess.STDOUT,check=True)
 binary=a.restore_to/'build/numivivo-omics';execution=json.loads((a.restore_to/'axis-execution.json').read_text())
 assert sha(binary)==execution['binarySHA256']
 with (a.out/'native-reopen.stdout').open('xb') as out,(a.out/'native-reopen.stderr').open('xb') as err:
  subprocess.run([str(binary),'singlecell-cell-axis-verify',str(a.restore_to/'axis')],stdout=out,stderr=err,check=True)
 assert (a.out/'native-reopen.stdout').read_bytes()==(a.restore_to/'axis.stdout').read_bytes()
 import tarfile
 with tarfile.open(a.archive,'r:gz') as t:members=json.load(t.extractfile('members.json'))
 assert {str(f.relative_to(a.restore_to)) for f in a.restore_to.rglob('*') if f.is_file()}==set(members)
 for name,item in members.items():assert sha(a.restore_to/name)==item['SHA256'],name
 assert sha(a.archive)==expected
 write(a.out/'summary.json',dict(status='passed-fresh-extraction-and-native-axis-reopen',archiveSHA256=expected,members=len(members),binarySHA256=sha(binary),cells=execution['cells'],allRestoredBytesExact=True,receiptByteIdentical=True,countStreamReplayed=False,predictionFitOrScoring=False))
 repo=Path(__file__).resolve().parents[4]
 sources=[Path(__file__),Path(__file__).parent/'parse_support.py']
 recipe=dict(schemaVersion=1,sourceFiles=[dict(path=str(f.relative_to(repo)),bytes=f.stat().st_size,SHA256=sha(f),rawUTF8=f.read_text()) for f in sources])
 (a.out/'recipe.json.gz').write_bytes(gzip.compress((json.dumps(recipe,sort_keys=True,separators=(',',':'))+'\n').encode(),mtime=0))
 records=[]
 for f in sorted(a.out.iterdir()):
  raw=f.read_bytes();compressed=f.name=='recipe.json.gz';decoded=gzip.decompress(raw) if compressed else raw
  records.append(dict(sourcePath=f.name.removesuffix('.gz') if compressed else f.name,sourceBytes=len(decoded),sourceSHA256=hashlib.sha256(decoded).hexdigest(),storedPath=f.name,storedBytes=len(raw),storedSHA256=sha(f),gzipEncoded=compressed,bundledSourceFiles=compressed))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records))
 subprocess.run([sys.executable,str(a.dependencies_repo/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'),str(a.out)],check=True)
 print((a.out/'summary.json').read_text())

if __name__=='__main__':main()
