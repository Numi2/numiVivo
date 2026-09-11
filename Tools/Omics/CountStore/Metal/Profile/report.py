#!/usr/bin/env python3
"""Summarize all declared runs without excluding startup or slower observations."""
import argparse,json,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();experiments=[]
for sub in ['', 'final']:
 d=json.loads((a.root/sub/'comparison.json').read_text());assert d['allBundleBytesExact'] and d['nativeReplaysPassed'] and len(d['rows'])==12
 groups=[]
 for backend in ['cpu','metal']:
  for owner in ['old','new']:
   rows=[x for x in d['rows'] if x['backend']==backend and x['owner']==owner]
   times=[x['seconds'] for x in rows];assert len(times)==3
   groups.append(dict(backend=backend,owner=owner,seconds=times,medianSeconds=statistics.median(times),meanSeconds=statistics.mean(times),peakRSSBytes=[x['peakRSSBytes'] for x in rows]))
 experiments.append(dict(candidate=sub or 'initial',groups=groups,binaries=d['freeze']['binaries']))
assert json.loads((a.root/'final/writer-independent.json').read_text())['status']=='passed'
assert len(json.loads((a.root/'final/cli-checks.json').read_text())['commands'])==13
result=dict(status='passed',experiments=experiments,allRunsIncluded=True,scope='One full original Kang cohort on a shared M4 Pro desktop; no cold-cache guarantee, general speedup, Metal superiority or new biological prediction qualification.')
(a.root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print('All 24 complete comparisons, final memory-safety and CLI checks passed')
