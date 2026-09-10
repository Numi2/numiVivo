#!/usr/bin/env python3
"""Summarize all declared populations, normalization policies and exceptions."""
import argparse, json
from pathlib import Path
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
a = p.parse_args(); r = a.root
protocol = json.loads((r/'protocol.json').read_text())
cases = []; run_count = 0; diagnostics = []
for case in protocol['cases']:
    d = r/case['id']; row = dict(case)
    if not (d/'comparison.json').exists():
        assert case['expectedReplicationGate'] == 'insufficient'
        assert 'insufficient biological replication after selection' in (d/'publish.log').read_text()
        row['status'] = 'unavailable-insufficient-independent-animals'
    else:
        comparison = json.loads((d/'comparison.json').read_text())
        assert comparison['independentBHVerified']
        row['status'] = 'completed-descriptive-comparison'
        row['comparisons'] = comparison['results']
        runs = json.loads((d/'reference/runs.json').read_text())
        assert len(runs) == 6
        for method, run in runs.items():
            run_count += 1
            if run['status'] != 'completed' or run['warnings'] or run['messages']:
                diagnostics.append(dict(case=case['id'], method=method, **run))
        row['modelCheck'] = json.loads((d/'model-check.json').read_text())
    cases.append(row)
model = json.loads((r/'model-check.json').read_text())
prob = json.loads((r/'probability-diagnostic.json').read_text())
result = dict(status='third-real-animal-study-complete-with-retained-numerical-exceptions',
    nativeSourceCommit='3823fb121d41ccd82efb0683bb4a787bde662a90',
    rdaURL='https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/muscData/Crowell19_4vs4.Rda',
    protocol=protocol, referenceFits=run_count, referenceDiagnostics=diagnostics,
    referenceBetaNonconvergence=sum(v['referenceBetaNonconvergence'] or 0 for c in cases for v in c.get('comparisons', {}).values()),
    testedGenePopulationPairs=sum(c['testedGenes'] for c in model['cases']),
    conditionalFits=sum(c['conditionalFits'] for c in model['cases']),
    originalStrictModelCheckStatus=model['status'], originalStrictModelCheckFailedFits=model['failedFits'],
    maximumStoredWaldArithmeticError=max(c['maximumStoredArithmeticError'] for c in prob['cases']),
    changedRefitBH005Decisions=sum(c['changedBH005Decisions'] for c in prob['cases']),
    cases=cases,
    qualification='Descriptive experimental comparison. Different testing families, estimators and tests retained. No FDR calibration, posterior coverage, causal validation or speedup claim.')
(r/'summary.json').write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'cases','protocol'}}, indent=2))
