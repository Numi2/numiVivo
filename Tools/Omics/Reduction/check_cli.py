#!/usr/bin/env python3
"""Exercise real count publication, sparse reduction, verified replay and rejection."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--imported', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
commands = []
def run(*args, success=True):
    start = time.monotonic()
    r = subprocess.run([str(a.binary.resolve()), *map(str, args)], capture_output=True, text=True, timeout=600)
    number = len(commands)
    (a.out / f'{number:02d}.stdout').write_text(r.stdout)
    (a.out / f'{number:02d}.stderr').write_text(r.stderr)
    commands.append(dict(arguments=list(map(str,args)), exitCode=r.returncode, seconds=time.monotonic()-start, expectedSuccess=success))
    (a.out / 'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    assert (r.returncode == 0) == success, r.stderr[-4000:]
    if not success:
        assert r.returncode > 0, 'signal termination is not controlled rejection'
    return r
store = a.out / 'store'
counts = a.out / 'counts.json'
receipt = a.out / 'analysis.json'
plan = dict(schemaVersion=1, id='sparse-reduction-qualification', normalizationTarget=10000,
            filter=dict(minimumCounts=1,minimumDetectedFeatures=1), contrasts=[],
            reduction=dict(highlyVariableFeatures=2000,meanBins=20,components=20,maximumBasis=128,relativeResidualTolerance=1e-6,seed=7))
plan_path = a.out / 'plan.json'
plan_path.write_text(json.dumps(plan,indent=2)+'\n')
run('singlecell-run', a.imported / 'manifest.json', '--store', store, '--output', counts)
run('singlecell-analyze', counts, '--plan', plan_path, '--store', store, '--output', receipt)
run('singlecell-analysis-verify', receipt, '--store', store)
run('singlecell-analysis-export', receipt, '--store', store, '--output', a.out / 'report.json')
run('singlecell-analyze', counts, '--plan', plan_path, '--store', store, '--output', a.out / 'repeat.json')
assert json.loads(receipt.read_text()) == json.loads((a.out/'repeat.json').read_text())
plan['reduction']['maximumBasis'] = 20
short = a.out / 'insufficient-basis.json'
short.write_text(json.dumps(plan,indent=2)+'\n')
rejected = a.out / 'rejected.json'
failure = run('singlecell-analyze', counts, '--plan', short, '--store', store, '--output', rejected, success=False)
assert 'residual' in failure.stderr.lower(), failure.stderr
assert not rejected.exists()
result = dict(status='passed', checks=['count publication','reduction publication','native reconstruction','report export','deterministic receipt','controlled unconverged rejection without receipt'],
              binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(), commands=commands,
              qualification='Software reconstruction and numerical convergence; not biological qualification')
(a.out / 'checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='commands'},indent=2))
