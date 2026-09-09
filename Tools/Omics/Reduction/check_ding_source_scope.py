#!/usr/bin/env python3
"""Verify that experiment metadata does not authorize donor correction."""
import argparse,json,shlex,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--host',required=True)
for name in ['binary','hdf5','remote-root','protocol','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
protocol=json.loads(a.protocol.read_text());options=dict(protocol['integration'],covariate='donor',ridge=1,seed=7)
plan=a.out/'donor-plan.json';plan.write_text(json.dumps(dict(schemaVersion=1,inputKind='fitted',integration=options),indent=2)+'\n');remote_plan=a.remote_root/'donor-rejection-plan.json';target=a.remote_root/'rejected-donor'
subprocess.run(['rsync',str(plan),a.host+':'+str(remote_plan)],check=True)
command=['env','NUMIVIVO_HDF5_LIBRARY='+str(a.hdf5),str(a.binary),'singlecell-pca-integrate',str(a.remote_root/'pca'),'--plan',str(remote_plan),'--output',str(target)]
result=subprocess.run(['ssh',a.host,shlex.join(command)],capture_output=True,text=True);(a.out/'rejection.log').write_text(result.stdout+result.stderr)
assert result.returncode==65 and 'known selected covariate' in result.stderr,result.stderr
check=subprocess.run(['ssh',a.host,shlex.join(['test','!','-e',str(target)])],capture_output=True,text=True);assert check.returncode==0
(a.out/'checks.json').write_text(json.dumps(dict(status='passed',command=command,exitCode=result.returncode,unknownDonorRejected=True,noDestinationPublished=True,qualification='Source experiment IDs remain experimental strata and biological-sample identifiers; no known donor identity is invented.'),indent=2)+'\n')
