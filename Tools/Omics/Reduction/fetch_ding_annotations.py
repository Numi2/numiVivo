#!/usr/bin/env python3
"""Fetch the frozen SCP424 annotation snapshot through its documented public API."""
import argparse,csv,hashlib,io,json,urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
expected={'CellType':'a211f2a98c99ae286ecd25efbf859092d289f804a57d8b4c7b3a320184733eeb','Experiment':'f5fe806825d175d23731bbb3e99fce64739d294bf6fbed59f9e16b8339f2b7ef','Method':'351c2e588e990435cdafedeae85766f9a85ca686abf75f8f1f34bb9a3cf5b28e'}
def fetch(field):
 url='https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP424/annotations/'+field+'/cell_values?annotation_type=group&annotation_scope=study'
 request=urllib.request.Request(url,headers={'Accept':'text/plain'})
 with urllib.request.urlopen(request,timeout=60) as response:
  assert response.status==200;raw=response.read(5000001);assert len(raw)<=5000000
 digest=hashlib.sha256(raw).hexdigest();(a.out/(field+'.tsv')).write_bytes(raw)
 assert digest==expected[field],(field,'Source snapshot changed',digest)
 reader=csv.DictReader(io.StringIO(raw.decode('utf-8')),delimiter='\t');assert reader.fieldnames==['NAME',field];rows=list(reader)
 names={r['NAME'] for r in rows};assert len(rows)==len(names)==31021
 return dict(field=field,url=url,sha256=digest,bytes=len(raw),rows=len(rows)),names
with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(fetch,expected))
assert all(names==results[0][1] for _,names in results)
result=dict(status='passed',sourceStudy='SCP424',apiDocumentation='https://singlecell.broadinstitute.org/single_cell/api/v1',sourceDescription='https://singlecell.broadinstitute.org/single_cell/study/SCP424/single-cell-comparison-pbmc-data',files=[r for r,_ in results],allFieldCellIdentitiesExactlyAligned=True,authenticationUsed=False,qualification='Official public visualization metadata, not the authenticated raw meta.txt file download. CellType is the source study-level individual-analysis annotation, separately described from Harmony-derived cluster annotations. Values are source annotations, not newly inferred labels.')
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
