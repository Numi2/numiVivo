#!/usr/bin/env python3
"""Keep all new native bytes and complete outcomes; bind retained large source inputs."""
import argparse,gzip,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'ReferenceMapping/Logistic'))
from bundles import pack,restore
# Import explicitly to avoid the reference recipe's unrelated prepare module.
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['pack','restore']);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--archive',type=Path);a=p.parse_args()
 if a.action=='restore':
  info=json.loads((a.archive/'native-dependencies.json').read_text());assert sha(a.archive/'native-results.tar.gz')==info['archive']['SHA256'];restore(a.archive/'native-results.tar.gz',a.out)
  for name,record in info['externalH5AD'].items():
   source=a.study/record['studyPath'];assert source.resolve().is_relative_to(a.study.resolve()) and sha(source)==record['SHA256'] and source.stat().st_size==record['bytes'];target=a.out/name;assert target.resolve().is_relative_to(a.out.resolve()) and not target.exists();target.parent.mkdir(parents=True,exist_ok=True)
   copied=subprocess.run(['/bin/cp','-c',str(source),str(target)],capture_output=True).returncode==0
   if not copied:assert not target.exists();shutil.copyfile(source,target)
   assert sha(target)==record['SHA256']
  print('Restored complete native result directories');return
 a.out.mkdir(parents=True,exist_ok=False);stage=Path(tempfile.mkdtemp(prefix='native-archive-stage-',dir=a.study));dependencies={};candidates={}
 for folder in ['prepared','prediction-inputs']:
  for file in (a.study/folder).glob('*.h5ad'):candidates.setdefault(sha(file),file)
 for folder in ['native-aggregation','native-prediction']:
  for file in sorted((a.study/folder).rglob('*')):
   if not file.is_file():continue
   assert not file.is_symlink();name=str(file.relative_to(a.study))
   if file.suffix=='.h5ad':
    digest=sha(file);source=candidates[digest];dependencies[name]=dict(studyPath=str(source.relative_to(a.study)),SHA256=digest,bytes=file.stat().st_size)
   else:
    dest=stage/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(file,dest)
 record=pack(stage,a.out/'native-results.tar.gz');write(a.out/'native-dependencies.json',dict(archive=record,externalH5AD=dependencies,scope='All new native non-H5AD bytes stored once by content; every source H5AD retained separately and hash-bound.'))
 opened=subprocess.run(['/usr/sbin/lsof','+D',str(stage)],capture_output=True,text=True);assert opened.returncode==1 and not opened.stdout.strip() and not opened.stderr.strip();shutil.rmtree(stage)
 records=[];groups={}
 def add(group,file,name):
  raw=file.read_bytes();groups.setdefault(group,[]).append(dict(path=name,bytes=len(raw),SHA256=hashlib.sha256(raw).hexdigest(),rawUTF8=raw.decode()))
 def store(name,raw,bundled=False):
  data=gzip.compress(raw,compresslevel=9,mtime=0);stored=name.replace('/','--')+'.gz';(a.out/stored).write_bytes(data);records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=hashlib.sha256(raw).hexdigest(),storedPath=stored,storedBytes=len(data),storedSHA256=hashlib.sha256(data).hexdigest(),gzipEncoded=True,bundledSourceFiles=bundled))
 for f in sorted(a.study.iterdir()):
  if f.is_file() and f.suffix in {'.json','.log','.py','.soft'}:add('execution',f,'study/'+f.name)
 for folder in ['prepared','prediction-inputs','scores','scores-repeat']:
  for f in sorted((a.study/folder).iterdir()):
   if f.is_file() and f.suffix in {'.json','.log'}:add('execution',f,'study/'+str(f.relative_to(a.study)))
   elif f.is_file() and f.suffix=='.npz' and folder!='scores-repeat':store(str(f.relative_to(a.study)),f.read_bytes())
 for f in sorted(Path(__file__).resolve().parent.iterdir()):
  if f.is_file():add('recipe',f,'recipe/'+f.name)
 external=[]
 for folder in ['sources','prepared']:
  for f in sorted((a.study/folder).rglob('*')):
   if f.is_file() and (f.suffix=='.h5ad' or f.suffix=='.h5' or f.parent.name=='qc'):external.append(dict(path=str(f.relative_to(a.study)),bytes=f.stat().st_size,SHA256=sha(f)))
 write(a.out/'source-dependencies.json',dict(files=external,scope='Complete downloaded sources, every-barcode QC arrays and prepared H5ADs remain local; no downsampling or deletion of unique data.'))
 for name,items in groups.items():store(name+'.json',(json.dumps(dict(schemaVersion=1,sourceFiles=items),sort_keys=True,separators=(',',':'))+'\n').encode(),True)
 for name in ['native-results.tar.gz','native-dependencies.json','source-dependencies.json']:
  f=a.out/name;records.append(dict(sourcePath=name,sourceBytes=f.stat().st_size,sourceSHA256=sha(f),storedPath=name,storedBytes=f.stat().st_size,storedSHA256=sha(f),gzipEncoded=False,bundledSourceFiles=False))
 write(a.out/'manifest.json',dict(schemaVersion=1,records=records));print(json.dumps(dict(status='archived',members=len(records),storedBytes=sum(r['storedBytes'] for r in records),nativeLogicalBytes=record['decodedBytes'],nativeH5ADDependencies=len(dependencies))))
if __name__=='__main__':main()
