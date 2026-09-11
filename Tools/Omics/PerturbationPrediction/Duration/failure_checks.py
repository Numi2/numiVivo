#!/usr/bin/env python3
"""Exercise public CLI rejection and reconstruction, retaining each failed attempt."""
import argparse,copy,hashlib,json,subprocess
from pathlib import Path
from common import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();out=a.study/'failure-checks';out.mkdir(exist_ok=False);src=a.study/'inputs/D34';reference=a.study/'native/D34/model';base=json.loads((src/'query.json').read_text());model=json.loads((reference/'model.json').read_text());records=[]
 cases=[('extrapolation',lambda p:p.update(hours=[37])),('duplicate-hours',lambda p:p.update(hours=[1,1])),('seen-donor',lambda p:p['mapping']['samples'][0].update(donorID=model['curves'][0]['donorID'],biologicalReplicateID=model['curves'][0]['donorID'])),('treated-query',lambda p:p['mapping']['samples'][0].update(condition='IFNB:2h')),('namespace',lambda p:p.update(featureNamespace='unmatched')),('units',lambda p:p['mapping'].update(countUnit='readCount')),('unknown-setting',lambda p:p.update(extrapolate=True))]
 for name,mutate in cases:
  plan=copy.deepcopy(base);mutate(plan);file=out/(name+'.json');write(file,plan);destination=out/(name+'-output')
  cmd=[str(a.binary),'singlecell-duration-predict',str(src/'query.h5ad'),'--plan',str(file),'--reference',str(reference),'--output',str(destination)]
  r=subprocess.run(cmd,capture_output=True,text=True);(out/(name+'.log')).write_text(r.stdout+r.stderr);assert r.returncode!=0 and not destination.exists() and not list(out.glob('.numivivo-*')),name
  records.append(dict(name=name,exitCode=r.returncode,outputAbsent=True,temporaryOutputAbsent=True,planSHA256=sha(file)))
 changed=out/'changed-model';subprocess.run(['/bin/cp','-cR',str(reference),str(changed)],check=True)
 model['ridgeInverse'][0][0]+=1;payload=json.dumps(model,sort_keys=True,separators=(',',':')).encode();(changed/'model.json').write_bytes(payload);receipt=json.loads((changed/'receipt.json').read_text());receipt['result']['bytes']=list(hashlib.sha256(payload).digest());write(changed/'receipt.json',receipt)
 r=subprocess.run([str(a.binary),'singlecell-duration-verify',str(changed)],capture_output=True,text=True);(out/'changed-model.log').write_text(r.stdout+r.stderr);assert r.returncode!=0 and 'model does not reconstruct' in r.stderr+r.stdout
 records.append(dict(name='modified-coefficient-with-updated-result-hash',exitCode=r.returncode,hashConsistentButReconstructionRejected=True))
 write(out/'summary.json',dict(status='passed',checks=records,binarySHA256=sha(a.binary),scriptSHA256=sha(Path(__file__))))
 print(json.dumps(dict(status='passed',checks=len(records))),flush=True)
if __name__=='__main__':main()
