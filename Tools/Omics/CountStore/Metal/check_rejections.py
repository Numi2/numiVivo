#!/usr/bin/env python3
"""Check no-publication failures and receipt binding with the real native owner."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--store',type=Path,required=True);p.add_argument('--normalized',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
receipt=(a.normalized/'receipt.json').read_bytes();results=[]
def reject(label,args,reason):
 p=subprocess.run([str(a.binary),*map(str,args)],capture_output=True)
 (a.out/(label+'.stdout')).write_bytes(p.stdout);(a.out/(label+'.log')).write_bytes(p.stderr)
 assert p.returncode==65 and reason in p.stderr.decode(),(label,p.returncode,p.stderr)
 results.append(dict(label=label,returnCode=p.returncode,expectedReason=reason))
for target in ['0','-1','nan','1000000001','0.5']:
 dest=a.out/('invalid-'+target)
 reject('target-'+target,['normalize-count-store',a.store,target,'metal-fp32',dest],
        'Metal normalization target' if target=='0.5' else 'count store normalization target')
 assert not dest.exists()
reject('overwrite',['normalize-count-store',a.store,'10000','metal-fp32',a.normalized],'destination exists')
assert (a.normalized/'receipt.json').read_bytes()==receipt
for label,change in [('device',lambda r:r['execution'].__setitem__('registryID',r['execution']['registryID']^1)),
                     ('profile',lambda r:r['execution'].__setitem__('numericalProfile','changed')),
                     ('backend',lambda r:r['execution'].__setitem__('backend','cpu-fp64'))]:
 dest=a.out/label;dest.mkdir()
 for name in ['values.bin','metadata.json','input-receipt.json']:os.link(a.normalized/name,dest/name)
 r=json.loads(receipt);change(r);(dest/'receipt.json').write_text(json.dumps(r,sort_keys=True,separators=(',',':'),ensure_ascii=False))
 reject(label,['verify-normalized-count-store',dest,a.store],
        'normalized count receipt' if label=='backend' else 'normalized count reconstruction differs')
assert (a.normalized/'receipt.json').read_bytes()==receipt
result=dict(status='passed',rejections=results,originalReceiptUnchanged=True,unavailableGPU='not exercised; both checked hosts expose physical Apple GPUs',checkerSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(status='passed',rejections=len(results))))
