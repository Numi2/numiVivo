#!/usr/bin/env python3
"""Reuse arithmetic checks only after exact output equality across the guard repair."""
import argparse,hashlib,json,shutil
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--initial',type=Path,required=True);p.add_argument('--final',type=Path,required=True);a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
old=json.loads((a.initial/'complete.json').read_text());new=json.loads((a.final/'complete.json').read_text());assert len(old)==len(new)==58
old={(r['case'],r['method']):r for r in old};results=[]
for r in new:
 key=(r['case'],r['method']);before=old[key];d=a.final/r['case']/r['method'];original=a.initial/r['case']/r['method']
 assert before['requestSHA256']==r['requestSHA256'] and before['inputSHA256']==r['inputSHA256'] and before['referenceSHA256']==r['referenceSHA256']
 for k in ['outputSHA256','logicalSHA256']:assert before[k]==r[k],(key,k)
 assert sha((d/'native.json.gz').read_bytes())==r['outputSHA256'] and sha((original/'native.json.gz').read_bytes())==r['outputSHA256']
 reused=False
 if (original/'independent-check.json').exists():
  check=json.loads((original/'independent-check.json').read_text());assert check['outputSHA256']==r['outputSHA256']
  for name in ['independent-check.json','independent-gene-errors.json.gz']:
   target=d/name
   if target.exists():assert target.read_bytes()==(original/name).read_bytes()
   else:shutil.copyfile(original/name,target)
  reused=True
 results.append(dict(case=r['case'],method=r['method'],outputSHA256=r['outputSHA256'],originalBinarySHA256=before['binarySHA256'],finalBinarySHA256=r['binarySHA256'],exactlyEqual=True,independentCheckReused=reused))
(a.final/'requalification-equality.json').write_text(json.dumps(results,sort_keys=True,indent=2)+'\n');print(json.dumps(dict(exactlyEqualArms=len(results),reusedChecks=sum(r['independentCheckReused'] for r in results))))
