"""Restore complete paired count artifacts and recheck all evidence without source I/O."""
import argparse,ctypes,importlib.util,json,os,shutil,subprocess,sys,tarfile,uuid
from pathlib import Path
from parse_support import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--evidence',type=Path,required=True);p.add_argument('--implementation-repo',type=Path,required=True);p.add_argument('--dependencies-repo',type=Path,required=True);p.add_argument('--source-study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--reuse-from',type=Path,help='Optional completed study; APFS copies require exact archive hashes and separate inodes');a=p.parse_args()
 evidence=a.evidence.resolve();repository=a.implementation_repo.resolve();dependencies=a.dependencies_repo.resolve()
 result=json.loads((evidence/'summary.json').read_text());assert result['status']=='passed-complete-paired-count-evidence' and result['records']==1373870697
 assert not a.out.exists() and not a.out.is_symlink()
 outer=dependencies/'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'
 subprocess.run([sys.executable,str(outer),str(evidence)],check=True)
 archive_module=dependencies/'Tools/Omics/PerturbationPrediction/Duration/archive.py';sys.path.insert(0,str(archive_module.parent))
 spec=importlib.util.spec_from_file_location('qualified_count_archive',archive_module);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 declared=json.loads((evidence/'dependencies.json').read_text())['required'];assert len(declared)==3
 archives=[]
 for item in declared:
  root=repository if item['role']=='axis-and-native-executable' else dependencies
  archive=root/item['path']/item['archive'];assert sha(archive)==item['archiveSHA256']
  if item['role']=='axis-and-native-executable':
   subprocess.run([sys.executable,str(outer),str(archive.parent)],check=True)
   archives.append((archive,item['archiveSHA256'],module.verify(archive)))
 assert len(archives)==1
 archive=evidence/'results.tar.gz';expected_digest=json.loads((evidence/'contents.json').read_text())['archiveSHA256'];assert sha(archive)==expected_digest
 archives.append((archive,expected_digest,module.verify(archive)))
 expected={}
 for _,_,members in archives:
  for name,item in members.items():assert name not in expected or expected[name]==item,('Conflicting archive path',name);expected[name]=item
 staging=a.out.absolute().with_name(a.out.name+'.incomplete-'+uuid.uuid4().hex);staging.mkdir(parents=True,exist_ok=False)
 clonefile=getattr(ctypes.CDLL(None,use_errno=True),'clonefile',None)
 if clonefile is not None:clonefile.argtypes=[ctypes.c_char_p,ctypes.c_char_p,ctypes.c_int];clonefile.restype=ctypes.c_int
 objects={};extracted=cloned=reused=0
 for archive,_,members in archives:
  byhash={}
  for name,item in members.items():byhash.setdefault(item['SHA256'],[]).append((name,item))
  with tarfile.open(archive,'r|gz') as t:
   for entry in t:
    digest=entry.name.removeprefix('objects/')
    for name,item in byhash.get(digest,[]):
     target=staging/name
     if target.exists():assert sha(target)==digest;continue
     target.parent.mkdir(parents=True,exist_ok=True);previous=objects.get(digest);copied=False
     candidate=previous or (a.reuse_from/name if a.reuse_from else None)
     if candidate is not None and clonefile is not None and candidate.is_file() and not candidate.is_symlink() and sha(candidate)==digest:
      copied=clonefile(os.fsencode(candidate),os.fsencode(target),0)==0
      if copied:
       assert candidate.stat().st_ino!=target.stat().st_ino
       if previous is None:reused+=1
       else:cloned+=1
      elif target.exists():target.unlink()
     if not copied:
      assert shutil.disk_usage(staging).free>item['bytes']+200000000,'Insufficient capacity for fresh restore with reserve'
      source=previous.open('rb') if previous is not None else t.extractfile(entry)
      with source,target.open('xb') as destination:shutil.copyfileobj(source,destination,1048576)
      extracted+=1
     assert target.stat().st_size==item['bytes'] and sha(target)==digest;target.chmod(item['mode']);objects[digest]=target
 assert {str(f.relative_to(staging)) for f in staging.rglob('*') if f.is_file()}==set(expected)
 checker=staging/'count-retention-recipe/check_counts.py';freeze=json.loads((staging/'count-retention-recipe/source-freeze.json').read_text())
 for item in freeze['sourceFiles']:assert sha(staging/'count-retention-recipe'/item['name'])==item['SHA256']
 with (staging/'restoration-check.log').open('x') as log:
  subprocess.run([sys.executable,str(checker),'--source-study',str(a.source_study),'--work',str(staging),'--dependencies-repo',str(dependencies),'--out',str(staging/'restoration-check.json')],stdout=log,stderr=subprocess.STDOUT,check=True)
 assert json.loads((staging/'restoration-check.json').read_text())==result
 for name,item in expected.items():assert sha(staging/name)==item['SHA256'],name
 assert all(sha(archive)==digest for archive,digest,_ in archives)
 record=dict(status='passed-count-restoration-and-offline-check',members=len(expected),directExtractions=extracted,withinDestinationAPFSCopies=cloned,sourceReuseAPFSCopies=reused,allRestoredHashesExact=True,completeOfflineResultIdentical=True,archives=[dict(path=str(archive),SHA256=digest) for archive,digest,_ in archives],nativeCountReexecuted=False,networkSourceReplayed=False,predictionFitOrScoring=False)
 write(staging/'restoration.json',record);os.rename(staging,a.out);print(json.dumps(record))

if __name__=='__main__':main()
