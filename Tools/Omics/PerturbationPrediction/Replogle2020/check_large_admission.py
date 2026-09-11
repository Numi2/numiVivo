#!/usr/bin/env python3
"""Real source stays bounded by the existing default storage and explicit work limits."""
import argparse,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();r=a.root;base=json.loads((r/'prepared-complete/rna-projection-plan.json').read_text());binary=r/'runtime/numivivo-omics';rows=[]
for label,edit,reason in [('default-storage',{'maximumOutputBytes':None},'byte allowance'),('insufficient-work',{'maximumElementVisits':100000000},'element visits')]:
 plan={**base,**edit}
 if plan.get('maximumOutputBytes') is None:plan.pop('maximumOutputBytes',None)
 path=r/(label+'-plan.json');path.write_text(json.dumps(plan)+'\n');dest=r/(label+'-rejected');assert not dest.exists()
 q=subprocess.run([str(binary),'singlecell-h5ad-project',str(r/'prepared-complete/original.h5ad'),'--plan',str(path),'--output',str(dest)],capture_output=True,env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib'))
 (r/(label+'.log')).write_bytes(q.stderr);assert q.returncode==65 and reason in q.stderr.decode() and not dest.exists()
 rows.append(dict(label=label,returnCode=q.returncode,outputAbsent=True,reason=reason))
(r/'large-admission-checks.json').write_text(json.dumps(dict(status='passed',cases=rows),indent=2)+'\n');print('Default storage and insufficient work reject full source without publication')
