#!/usr/bin/env python3
"""Require the complete original abundance/global-scale qualification family."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path,required=True)
args = parser.parse_args()
sha = lambda raw: hashlib.sha256(raw).hexdigest()
root = args.root
protocol = json.loads((root/'family/protocol.json').read_text())
assert protocol['protocolSHA256']==sha(Path(__file__).with_name('PROTOCOL.md').read_bytes())
runs = json.loads((root/'family/complete.json').read_text())
assert len(runs)==58
assert {(r['case'],r['method']) for r in runs}=={(r['case'],r['method']) for r in protocol['family']}
controls = json.loads((root/'controls/checks.json').read_text())
assert controls['status']=='passed' and controls['cases']==105
assert controls['binarySHA256']==protocol['binarySHA256']
assert controls['protocolSHA256']==protocol['protocolSHA256']
for filename,key in [('input.json','inputSHA256'),('native.json','nativeSHA256'),('reference.json','referenceSHA256')]:
    assert sha((root/'controls'/filename).read_bytes())==controls[key]
failures=[]
for run in runs:
    directory=root/'family'/run['case']/run['method']
    assert json.loads((directory/'run.json').read_text())==run
    raw=(directory/'native.json.gz').read_bytes()
    assert sha(raw)==run['outputSHA256'] and sha(gzip.decompress(raw))==run['logicalSHA256']
    assert run['binarySHA256']==protocol['binarySHA256']
    if run['status']!='passed':
        failures.append(run);continue
    raw=(directory/'comparison.json.gz').read_bytes()
    assert sha(raw)==run['comparisonSHA256']
    table=json.loads(gzip.decompress(raw))
    assert len(table)==run['genes'] and all(row['passed'] for row in table)
    assert run['failedGenes']==0 and run['exitCode']==0
    assert run['maximumAbundanceAbsoluteError']<=2e-8
    assert run['maximumAbundanceScaledScore']<=1e-7+1e-10
    assert run['maximumFinalScaledScore']<=1e-7+1e-10
    assert run['maximumRefitMeanRelativeError']<=2e-7
    assert max(run['maximumResidualRelativeError'].values())<=2e-7
    assert run['scaleRelativeError']<=2e-7
    assert run['maximumAbundanceBracketWidth']<=2e-12
summary=dict(status='passed' if not failures else 'failed',arms=len(runs),
             geneArmFits=sum(r.get('genes',0) for r in runs),controls=controls['cases'],
             finiteControlledCases=sum(c['expected']=='finite' for c in controls['checks']),
             expectedControlledErrors=sum(c['expected']=='error' for c in controls['checks']),
             maximumControlledAbsoluteError=controls['maximumIndependentAbsoluteError'],
             referenceControlDisagreements=controls['referenceDisagreements'],
             abundanceScoreEvaluations=sum(r.get('abundanceScoreEvaluations',0) for r in runs),
             omittedGeneUpdates=sum(len(e) for r in runs for e in r.get('omittedIndices',[])),
             failures=failures,results=runs)
if not failures:
    for key in ['maximumAbundanceAbsoluteError','maximumAbundanceScaledScore','maximumFinalScaledScore',
                'maximumRefitMeanRelativeError','scaleRelativeError','maximumAbundanceScoreEvaluations',
                'maximumAbundanceBracketWidth']:
        summary[key]=max(r[key] for r in runs)
    summary['maximumResidualRelativeError']={key:max(r['maximumResidualRelativeError'][key] for r in runs)
                                           for key in ['deviance','degreesOfFreedom','quasiDispersion']}
(root/'summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['results','failures']},indent=2))
raise SystemExit(0 if summary['status']=='passed' else 1)
