#!/usr/bin/env python3
"""Compare published whole-map PCA with windowed PCA on complete prepared cohorts."""
import argparse,hashlib,json,re,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--baseline',type=Path,required=True)
p.add_argument('--root',type=Path,required=True)
a=p.parse_args();results={}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
for cohort in ['baron','hagai']:
 prepared=a.root/(cohort+'-input');baseline=a.root/(cohort+'-baseline');baseline.mkdir(exist_ok=False)
 plan=dict(schemaVersion=1,mapping=json.loads((prepared/'mapping.json').read_text()),contrasts=[],reduction={})
 plan_path=baseline/'plan.json';plan_path.write_text(json.dumps(plan,indent=2)+'\n')
 args=[str(a.baseline),'singlecell-h5ad-pseudobulk',str(prepared/'prepared.h5ad'),'--plan',str(plan_path),'--output',str(baseline/'bundle')]
 start=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',*args],capture_output=True,text=True)
 (baseline/'run.log').write_text(r.stdout+r.stderr)
 peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr)
 baseline_command=dict(arguments=args,exitCode=r.returncode,seconds=time.monotonic()-start,maximumResidentBytes=int(peak[1]) if peak else None)
 (baseline/'command.json').write_text(json.dumps(baseline_command,indent=2)+'\n');assert r.returncode==0,r.stderr
 out=a.root/cohort
 subprocess.run([sys.executable,str(Path(__file__).with_name('check_streamed_cli.py')),'--binary',str(a.binary),'--prepared',str(prepared),'--out',str(out)],check=True)
 old=json.loads((baseline/'bundle/report.json').read_text());new=json.loads((out/'bundle/report.json').read_text())
 assert old['reduction']==new['reduction'],cohort+' changed PCA output'
 for key in ['metadata','quality','pseudobulk','canonicalNonzeros']:assert old[key]==new[key],key
 storage=new['reductionStorage'];assert storage['method']=='three-source-passes-windowed-selected-COO-v2'
 assert storage['cacheBytes']>16777216
 for key in ['cacheBytes','selectedEntries','cacheFingerprint','entryVisits','sourcePasses','options']:assert old['reductionStorage'][key]==storage[key],key
 native=json.loads((out/'checks.json').read_text())
 results[cohort]=dict(status='passed',sourceSHA256=sha(prepared/'prepared.h5ad'),cells=len(new['metadata']['cells']),features=len(new['metadata']['features']),
  allPCAFieldsExactlyUnchanged=True,allSourceStatisticsExactlyUnchanged=True,storage=storage,baseline=baseline_command,windowed=native)
 (a.root/'comparison.json').write_text(json.dumps(dict(status='in-progress',cohorts=results),indent=2)+'\n')
(a.root/'comparison.json').write_text(json.dumps(dict(status='passed',binarySHA256=sha(a.binary),baselineSHA256=sha(a.baseline),checkerSHA256=sha(Path(__file__)),
 cohorts=results,qualification='Two complete real cohorts; same-host single-run timing observations, not a general speedup or million-cell qualification'),indent=2)+'\n')
