#!/usr/bin/env python3
"""Run every predeclared Baron donor fold with original source snapshots."""
import argparse,gzip,json,os,re,shutil,subprocess,time
from pathlib import Path
from prepare import sha,write
from bundles import pack

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,h in freeze['files'].items():assert sha(a.inputs/name)==h,name
 commands=[];archives=[]
 write(a.out/'execution-start.json',dict(binarySHA256=sha(a.binary),inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),runnerSHA256=sha(__file__),hdf5SHA256=sha(os.environ['NUMIVIVO_HDF5_LIBRARY']),scoringStarted=False))
 def run(label,args):
  start=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
  peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr);commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-start,maximumResidentBytes=int(peak[1]) if peak else None));write(a.out/'commands.json',commands);assert r.returncode==0,(label,r.stderr[-2000:])
 for donor in ['human1','human2','human3','human4']:
  assert shutil.disk_usage(a.out).free>500000000,('capacity',donor)
  root=a.out/('scratch-'+donor);(root/'inputs').mkdir(parents=True)
  for group in ['train','query']:
   spec=freeze['sources'][donor+'/'+group];target=root/'inputs'/(group+'.h5ad')
   with gzip.open(a.inputs/spec['gzipPath'],'rb') as f,target.open('xb') as g:shutil.copyfileobj(f,g,1048576)
   assert sha(target)==spec['decodedSHA256']
  run(donor+'-fit',['singlecell-reference-fit',root/'inputs/train.h5ad','--plan',a.inputs/donor/'fit.json','--output',root/'reference'])
  run(donor+'-verify-fit',['singlecell-reference-verify',root/'reference'])
  args=['singlecell-reference-map',root/'inputs/query.h5ad','--plan',a.inputs/donor/'query.json','--reference',root/'reference','--output']
  run(donor+'-map',args+[root/'mapped']);run(donor+'-verify-map',['singlecell-reference-map-verify',root/'mapped'])
  if donor=='human1':
   run(donor+'-repeat',args+[root/'repeat']);assert sha(root/'mapped/report.json')==sha(root/'repeat/report.json')
  model=json.loads((root/'reference/model.json').read_text());logistic=model['logistic'];assert logistic['gradientMaximum']<=model['plan']['logistic']['gradientTolerance']
  write(a.out/(donor+'-optimizer.json'),logistic)
  record=pack(root,a.out/(donor+'.tar.gz'));record.update(donor=donor,path=donor+'.tar.gz',archiveVerified=True)
  opened=subprocess.run(['/usr/sbin/lsof','+D',str(root)],capture_output=True,text=True);assert opened.returncode==1 and not opened.stdout.strip() and not opened.stderr.strip()
  archives.append(record);write(a.out/'archive-progress.json',archives)
  # The original compressed inputs stay intact. Only this invocation's decoded
  # scratch is removed, after all unique objects and reconstruction paths verify.
  shutil.rmtree(root)
  print(json.dumps(dict(donor=donor,iterations=logistic['iterations'],gradientMaximum=logistic['gradientMaximum'],archiveBytes=record['bytes'],availableBytes=shutil.disk_usage(a.out).free)),flush=True)
 write(a.out/'prediction-freeze.json',dict(status='completed',allFourDonors=True,archives=archives,commands=commands,inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),binarySHA256=sha(a.binary),scoringStarted=False))
if __name__=='__main__':main()
