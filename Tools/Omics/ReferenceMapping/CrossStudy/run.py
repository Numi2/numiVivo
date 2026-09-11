#!/usr/bin/env python3
"""Execute both complete source-panel transfers and freeze before label scoring."""
import argparse,json,os,platform,re,shutil,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'Logistic'))
from prepare import sha,write
from bundles import pack

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 freeze=json.loads((a.inputs/'input-freeze.json').read_text())
 for name,h in freeze['files'].items():assert sha(a.inputs/name)==h,name
 commands=[];archives=[]
 start=dict(binarySHA256=sha(a.binary),inputFreezeSHA256=sha(a.inputs/'input-freeze.json'),runnerSHA256=sha(__file__),hdf5SHA256=sha(os.environ['NUMIVIVO_HDF5_LIBRARY']),platform=platform.platform(),hardware=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),memoryBytes=int(subprocess.check_output(['sysctl','-n','hw.memsize'],text=True)),scoringStarted=False)
 write(a.out/'execution-start.json',start)
 def run(label,args):
  t=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
  commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-t,maximumResidentBytes=int(peak[1]) if peak else None));write(a.out/'commands.json',commands)
  print(json.dumps(commands[-1]),flush=True);assert r.returncode==0,(label,r.stderr[-2500:])
 for train,query,trainFile,queryFile in [('kang','ding','kang.h5ad','ding-query.h5ad'),('ding','kang','ding-train.h5ad','kang.h5ad')]:
  direction=train+'-to-'+query;assert shutil.disk_usage(a.out).free>1500000000
  root=a.out/('scratch-'+direction);root.mkdir()
  run(direction+'-fit',['singlecell-reference-fit',a.inputs/trainFile,'--plan',a.inputs/(train+'-fit.json'),'--output',root/'reference'])
  run(direction+'-verify-fit',['singlecell-reference-verify',root/'reference'])
  args=['singlecell-reference-map',a.inputs/queryFile,'--plan',a.inputs/(query+'-query.json'),'--reference',root/'reference','--output',root/'mapped']
  run(direction+'-map',args);run(direction+'-verify-map',['singlecell-reference-map-verify',root/'mapped'])
  model=json.loads((root/'reference/model.json').read_text());write(a.out/(direction+'-optimizer.json'),model['logistic'])
  record=pack(root,a.out/(direction+'.tar.gz'));record.update(direction=direction,path=direction+'.tar.gz',archiveVerified=True);archives.append(record);write(a.out/'archive-progress.json',archives)
  opened=subprocess.run(['lsof','+D',str(root)],capture_output=True,text=True);assert opened.returncode==1 and not opened.stdout.strip(),opened.stdout
  shutil.rmtree(root)
 write(a.out/'prediction-freeze.json',dict(status='completed',bothDirections=True,queryCells=68704,archives=archives,commands=commands,execution=start,scoringStarted=False))
if __name__=='__main__':main()
