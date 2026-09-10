#!/usr/bin/env python3
"""Check the full CLI's 2 MiB pseudobulk plan cap before reading a source matrix."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

p=argparse.ArgumentParser(description=__doc__)
for key in ('binary','source','plan','out'):p.add_argument('--'+key,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
original=a.plan.read_bytes();assert 131_072<len(original)<2_097_152
oversized=a.out/'over-two-mib.json';oversized.write_bytes(original+b' '*(2_097_153-len(original)))
destination=a.out/'unexpected-output'
result=subprocess.run([str(a.binary.resolve()),'singlecell-h5ad-pseudobulk',str(a.source),'--plan',str(oversized),'--output',str(destination)],capture_output=True,text=True)
(a.out/'rejection.log').write_text(result.stdout+result.stderr)
assert result.returncode==65 and 'Artifact-store I/O bound exceeded: over-two-mib.json' in result.stderr
assert not destination.exists()
(a.out/'report.json').write_text(json.dumps(dict(status='passed',originalPlanBytes=len(original),
    rejectedBytes=oversized.stat().st_size,returnCode=result.returncode,noOutputPublished=True,
    originalPlanSHA256=hashlib.sha256(original).hexdigest()),indent=2)+'\n')
print((a.out/'report.json').read_text())
