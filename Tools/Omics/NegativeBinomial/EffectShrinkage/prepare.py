#!/usr/bin/env python3
"""Freeze measured treatment plans and exact baseline receipts before MAP fits."""
import argparse, gzip, hashlib, json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
p.add_argument('--baseline-home',type=Path,required=True)
a=p.parse_args();a.root.mkdir(exist_ok=True,parents=True)
records=[]
for study in ['kang','hagai']:
 for mode in ['default','active']:
  source=a.baseline_home/(f'numivivo-{study}-r-native-current-20260909' if mode=='default' else f'numivivo-{study}-nb-support-20260909/active')
  out=a.root/(study+'-'+mode);out.mkdir(exist_ok=False)
  report=(source/'report.json').read_bytes();plan=json.loads((source/'plan.json').read_text())
  assert len(plan['contrasts'])==1
  options=plan['contrasts'][0]['negativeBinomialOptions'];assert 'effectPriorStandardDeviationLog2' not in options
  options['effectPriorStandardDeviationLog2']=1.0 # Predeclared, never tuned to results.
  (out/'plan.json').write_text(json.dumps(plan,sort_keys=True)+'\n')
  (out/'baseline-report.json.gz').write_bytes(gzip.compress(report,mtime=0))
  for name in ['receipt.json','plan.json']:
   (out/('baseline-'+name)).write_bytes((source/name).read_bytes())
  receipt=json.loads((source/'receipt.json').read_text())
  records.append(dict(study=study,mode=mode,baseline=str(source),reportSHA256=hashlib.sha256(report).hexdigest(),sourceSHA256=bytes(receipt['source']['bytes']).hex(),priorSDLog2=1.0))
(a.root/'protocol.json').write_text(json.dumps(dict(prior='zero-centered normal on requested contrast',priorSDLog2=1.0,qualification='conditional numerical validation; no posterior coverage, power, or FDR qualification',records=records),sort_keys=True,indent=2)+'\n')
