#!/usr/bin/env python3
"""Execute both native Visium formats, reconstruct, and retain restorable bundles."""
import argparse,hashlib,json,os,re,shutil,subprocess,tarfile,time
from pathlib import Path
from prepare import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 a.out.mkdir(parents=True,exist_ok=False);freeze=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,h in freeze['files'].items():assert sha(a.inputs/name)==h,name
 commands=[];archives=[]
 def run(label,args):
  start=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True)
  (a.out/(label+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
  commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-start,maximumResidentBytes=int(peak[1]) if peak else None));write(a.out/'commands.json',commands)
  assert r.returncode==0,(label,r.stderr[-1500:])
 for label in ['legacy','header']:
  free=shutil.disk_usage(a.out).free
  assert free>1000000000,('capacity-before-native',label,free)
  root=a.out/('scratch-'+label)
  run(label+'-import',['multiassay-visium-import',a.inputs/'outs','--plan',a.inputs/(label+'-plan.json'),'--output',root])
  run(label+'-verify',['multiassay-verify',root])
  assert {p.name for p in root.iterdir()}=={'original.h5','positions.csv','plan.json','dataset.json','dataset.h5mu','receipt.json'}
  assert all(p.is_file() and not p.is_symlink() for p in root.iterdir())
  members={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in sorted(root.iterdir())};archive=a.out/(label+'.tar.gz')
  with tarfile.open(archive,'w:gz',compresslevel=6) as t:
   for name in members:t.add(root/name,arcname=name,recursive=False)
  with tarfile.open(archive,'r:gz') as t:
   assert {m.name for m in t.getmembers()}==set(members)
   for m in t.getmembers():
    assert m.isfile();h=hashlib.sha256();n=0
    with t.extractfile(m) as f:
     for b in iter(lambda:f.read(1048576),b''):h.update(b);n+=len(b)
    assert n==members[m.name]['bytes'] and h.hexdigest()==members[m.name]['SHA256']
  opened=subprocess.run(['/usr/sbin/lsof','+D',str(root)],capture_output=True,text=True)
  assert opened.returncode==1 and not opened.stdout.strip() and not opened.stderr.strip()
  archives.append(dict(label=label,path=archive.name,bytes=archive.stat().st_size,SHA256=sha(archive),members=members,archiveVerified=True,uncompressedBytes=sum(x['bytes'] for x in members.values()),scratchPath=str(root)))
  write(a.out/'archive-progress.json',archives)
  # Delete only this finished invocation's exact scratch after a complete
  # restorable archive check and no open handles. Original inputs stay intact.
  shutil.rmtree(root)
  print(json.dumps(dict(completed=label,archiveBytes=archive.stat().st_size,availableBytes=shutil.disk_usage(a.out).free)),flush=True)
 for name in ['dataset.json','dataset.h5mu']:assert archives[0]['members'][name]==archives[1]['members'][name]
 write(a.out/'execution.json',dict(status='passed',inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),binarySHA256=sha(a.binary),runnerSHA256=sha(__file__),hdf5LibrarySHA256=sha(os.environ['NUMIVIVO_HDF5_LIBRARY']),commands=commands,archives=archives,formatDatasetBytesExact=True))
if __name__=='__main__':main()
