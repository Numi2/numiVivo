#!/usr/bin/env python3
"""Freeze empirical-prior plans against three measured study baselines."""
import argparse,csv,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ['root','fixed_prior_root','crowell_root']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();records=[]
def load(p):return json.loads(gzip.decompress(p.read_bytes())) if p.suffix=='.gz' else json.loads(p.read_text())
cases=[(name,a.fixed_prior_root/name/'native') for name in ['kang-default','kang-active','hagai-default']]
cases += [('crowell-'+c['id'],a.crowell_root/c['id']/'native') for c in load(a.crowell_root/'protocol.json')['cases'] if c['expectedReplicationGate']=='ready']
for name,baseline in cases:
 report=load(baseline/'report.json.gz');plan=load(baseline/'plan.json.gz');archive=load(baseline/'archive.json')
 for request in plan['contrasts']:
  assert request['negativeBinomialOptions'].pop('effectPriorStandardDeviationLog2')==1
  request['negativeBinomialOptions']['effectPriorEstimation']='weightedUpperQuantile'
 d=a.root/name;d.mkdir(exist_ok=False);raw=json.dumps(plan,sort_keys=True,indent=2)+'\n';(d/'plan.json').write_text(raw)
 with (d/'prior-input.tsv').open('w') as h:
  writer=csv.writer(h,delimiter='\t');writer.writerow(['featureIndex','featureID','effectLog2','mean','trendDispersion'])
  cohort=report['contrasts'][0]
  for f,diagnostic in zip(cohort['features'],cohort['negativeBinomial']['features'],strict=True):
   if f['status']=='tested' and 'supportResolution' not in diagnostic:
    writer.writerow([f['featureIndex'],f['featureID'],f['log2FoldChange'],f['meanNormalizedCount'],diagnostic['trendDispersion']])
 records.append(dict(id=name,baselineReport=str(baseline/'report.json.gz'),baselineReportSHA256=hashlib.sha256((baseline/'report.json.gz').read_bytes()).hexdigest(),source=archive['source'],sourceSHA256=archive['sourceSHA256'],planSHA256=hashlib.sha256(raw.encode()).hexdigest(),testedFeatures=report['contrasts'][0]['testedFeatures']))
protocol=dict(records=records,baselineSourceCommit='cf2008ec240312b5baf95d75ceca73d8e5678441',method='weightedUpperQuantile',referencePolicy='All tested full-design genes; abs(log2 effect)<10; weights 1/(1/mean + trend dispersion); weighted absolute 95th percentile divided by qnorm(.975); minimum 20 genes; SD floor .01',unavailable='Crowell CPE remains replication-ineligible; no empirical fit attempted',qualification='Numerical empirical-prior and conditional MAP checks, not prior uncertainty, posterior coverage, risk, power or FDR qualification')
(a.root/'protocol.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n');print(json.dumps(dict(cases=len(records),testedFeatures=sum(r['testedFeatures'] for r in records))))
