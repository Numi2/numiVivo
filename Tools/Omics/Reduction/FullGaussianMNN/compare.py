#!/usr/bin/env python3
"""Summarize the fixed intervention against retained original and local64 results."""
import argparse,hashlib,json
from pathlib import Path

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser()
 for n in ('root','previous','original','out'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();rows=[];inputs={}
 for c in ('hagai','kang','ding'):
  paths={'fullGaussian':a.root/'evaluation'/c/'checks.json','local64':a.previous/'evaluation'/c/'checks.json','exactMatching':a.original/('evaluation-'+c)/'checks.json'}
  docs={k:read(v) for k,v in paths.items()};inputs.update({str(v):sha(v) for v in paths.values()});baseline=docs['fullGaussian']['baseline'];m=docs['fullGaussian']['metrics']
  keys=('cellTypeBalancedAccuracy','conditionBalancedAccuracy','programSpearman','withinStratumProgramSpearman')
  metrics={key:{method:d['metrics'].get(key) for method,d in docs.items()} for key in keys if key in m}
  losses={t:baseline['cellTypeRecall'][t]-v for t,v in m.get('cellTypeRecall',{}).items()}
  rows.append(dict(cohort=c,cells=docs['fullGaussian']['cells'],gates={k:d['gates'] for k,d in docs.items()},metrics=metrics,perTypeRecallLossVersusOriginalUncorrected=losses,programs={k:d['metrics'].get('programs') for k,d in docs.items()} if c=='ding' else None))
 out=dict(status='compared-complete-original-cohorts',cohorts=rows,inputSHA256=inputs,checkerSHA256=sha(Path(__file__)),scope='Fixed inspected-cohort development intervention. Original uncorrected biological baseline and unchanged gate margins; exact MNN and prior local64 comparisons are retained separately. No parameter tuning or independent-study claim.')
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')
 print(json.dumps({c['cohort']:c['gates']['fullGaussian'] for c in rows}))
if __name__=='__main__':main()
