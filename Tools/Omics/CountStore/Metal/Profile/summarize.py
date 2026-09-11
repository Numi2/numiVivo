#!/usr/bin/env python3
"""Extract inclusive main-thread stack samples, not timing or GPU utilization."""
import argparse,hashlib,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();rows=[]
for label in ['cpu','metal']:
 path=a.root/(label+'-sample.txt');text=path.read_text();graph=text.split('Call graph:\n',1)[1].split('Total number in stack',1)[0]
 nodes=[];stack=[]
 for line in graph.splitlines():
  m=re.match(r'^([ +!:|]+)(\d+) (.+)$',line)
  if not m:continue
  depth=len(m[1]);name=m[3]
  while stack and stack[-1]['depth']>=depth:stack.pop()
  node=dict(depth=depth,count=int(m[2]),name=name,ancestors=[x['name'] for x in stack]);nodes.append(node);stack.append(node)
 marker='static VivoH5ADCountStore.normalize('
 roots=[n for n in nodes if n['name'].startswith(marker)]
 assert len(roots)==1
 denominator=roots[0]['count']
 def sample_count(predicate):
  return sum(n['count'] for n in nodes if any(x.startswith(marker) for x in n['ancestors']) and predicate(n))
 counts=dict(sourceReconstruction=sample_count(lambda n:n['name'].startswith(('static VivoH5ADCountStore.build(', 'specialized static VivoH5ADCountStore.build('))),
             fileSnapshotHashing=sample_count(lambda n:n['name'].startswith('static VivoOmicsFileSnapshot.fingerprint(')),
             recordWriterInclusive=sample_count(lambda n:n['name'].startswith('VivoCountRecordWriter.append(')),
             gpuWait=sample_count(lambda n:'MTLCommandBuffer waitUntilCompleted' in n['name']))
 assert all(0<=x<=denominator for x in counts.values())
 rows.append(dict(backend=label,normalizationMainThreadSamples=denominator,inclusiveSamples=counts,
                  conditionalPercent={k:v/denominator*100 for k,v in counts.items()},SHA256=hashlib.sha256(path.read_bytes()).hexdigest()))
result=dict(status='passed',profiles=rows,scope='Inclusive main-thread samples conditioned on normalization; nested categories overlap. Not wall-clock stage durations or GPU utilization.')
(a.root/'profile-summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
