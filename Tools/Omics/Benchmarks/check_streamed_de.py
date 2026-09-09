#!/usr/bin/env python3
"""Compare streamed NB results to an existing native resident cohort report.

MEX can reorder cells by sample. Compare source membership by cell identity,
then require every remaining design, diagnostic and statistical field to agree.
"""
import argparse,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--resident-report',type=Path,required=True)
p.add_argument('--bundle',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
old=json.loads(a.resident_report.read_text());new=json.loads((a.bundle/'report.json').read_text())
old_cells=old['processed']['dataset']['cells'];new_cells=new['metadata']['cells']
def identity(cell):return (cell['sampleID'],cell['barcode'],cell.get('group'))
assert sorted(map(identity,old_cells))==sorted(map(identity,new_cells))
assert len(old['contrasts'])==len(new['contrasts'])>0
for left,right in zip(old['contrasts'],new['contrasts']):
    lobs=left['design']['observations'];robs=right['design']['observations']
    assert len(lobs)==len(robs)
    for l,r in zip(lobs,robs):
        li=l.pop('sourceCellIndices');ri=r.pop('sourceCellIndices')
        assert sorted(identity(old_cells[i]) for i in li)==sorted(identity(new_cells[i]) for i in ri)
    assert left==right,'numerical results or experimental design changed'
result=dict(status='passed-exact-DE-comparison-with-identity-resolved-membership',cells=len(new_cells),nonzeros=new['canonicalNonzeros'],contrasts=len(new['contrasts']),receipt=json.loads((a.bundle/'receipt.json').read_text()),
    qualification='All statistical and diagnostic fields exactly equal resident report. Source row indices differ because MEX groups source cells by sample; cell membership verified by sample/barcode/group identity.')
with a.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(result['status'])
