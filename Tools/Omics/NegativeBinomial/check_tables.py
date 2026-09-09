#!/usr/bin/env python3
"""Validate already-exported native NB tables against their complete report."""
import argparse,csv,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--report',type=Path,required=True);p.add_argument('--tables',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();report=json.loads(a.report.read_text());index=json.loads((a.tables/'index.json').read_text())
checked=0
for declaration,contrast in zip(index['contrasts'],report['contrasts'],strict=True):
    assert declaration['method']==contrast['method'] and declaration['id']==contrast['request']['id']
    if 'negativeBinomial' not in contrast:continue
    diagnostics=json.loads((a.tables/declaration['diagnostics']).read_text())
    assert diagnostics==contrast['negativeBinomial']
    with (a.tables/declaration['statistics']).open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==len(contrast['features'])
    for row,feature in zip(rows,contrast['features'],strict=True):
        assert row['feature_id']==feature['featureID'] and row['status']==feature['status']
        assert row['t']=='' and row['df']==''
        for key,column in [('zStatistic','z'),('pValue','p_value'),('adjustedPValue','BH_adjusted_p')]:
            if key in feature:assert float(row[column])==feature[key]
            else:assert row[column]==''
    checked+=1
assert checked>0
summary=dict(status='passed-NB-table-diagnostics-and-missing-values',contrasts=checked,
    reportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),indexSHA256=hashlib.sha256((a.tables/'index.json').read_bytes()).hexdigest())
with a.out.open('x') as f:json.dump(summary,f,indent=2);f.write('\n')
print(json.dumps(summary,indent=2))
