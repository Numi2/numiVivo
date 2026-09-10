#!/usr/bin/env python3
"""Verify the streamed annotation source cap without materializing a large fixture."""
import argparse
import json
from pathlib import Path
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=False)
source=a.out/'over-one-gib.h5ad'
with source.open('wb') as f: f.truncate(1_073_741_825)
plan=a.out/'plan.json'
plan.write_text(json.dumps(dict(schemaVersion=1,source=dict(bytes=[0]*32),provenance='Software admission-bound fixture',
    edits=[dict(path='uns/check',mode='add',value=dict(string=dict(shape=[],values=['fixture'])))])))
destination=a.out/'rejected.h5ad'
command=[str(a.binary.resolve()),'singlecell-h5ad-annotate',str(source),'--plan',str(plan),'--output',str(destination)]
result=subprocess.run(command,capture_output=True,text=True)
(a.out/'rejection.log').write_text(result.stdout+result.stderr)
assert result.returncode==65 and 'omics input type or bytes' in result.stderr
assert not destination.exists() and source.stat().st_size==1_073_741_825
# This is a generated sparse invalid-input fixture, never the experimental file.
source.unlink()
(a.out/'report.json').write_text(json.dumps(dict(status='passed',sourceBytes=1_073_741_825,
    expectedReturnCode=65,actualReturnCode=result.returncode,noOutputPublished=True,sparseFixtureRemoved=True),indent=2)+'\n')
print((a.out/'report.json').read_text())
