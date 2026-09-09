#!/usr/bin/env python3
"""Verify native full-source Baron publication/replay and scientific rejections."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--prepared', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
commands = []


def run(arguments, expected=None):
    start = time.monotonic()
    result = subprocess.run([str(a.binary.resolve()), *map(str, arguments)], capture_output=True, text=True, timeout=300)
    slot = len(commands)
    (a.out/f'{slot:02}.stdout').write_text(result.stdout)
    (a.out/f'{slot:02}.stderr').write_text(result.stderr)
    commands.append(dict(arguments=list(map(str, arguments)), exitCode=result.returncode,
                         seconds=time.monotonic()-start, expectedFailure=expected))
    (a.out/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    if expected is None:
        assert result.returncode == 0, result.stderr
    else:
        assert result.returncode > 0 and expected in result.stderr.lower(), result.stderr
    return result


source = a.prepared/'prepared.h5ad'
bundle = a.out/'bundle'
run(['singlecell-h5ad-pseudobulk', source, '--plan', a.prepared/'stream-plan.json', '--output', bundle])
run(['singlecell-h5ad-pseudobulk-verify', bundle])
repeated = a.out/'repeat'
run(['singlecell-h5ad-pseudobulk', source, '--plan', a.prepared/'stream-plan.json', '--output', repeated])
assert (bundle/'receipt.json').read_bytes() == (repeated/'receipt.json').read_bytes()
resident = a.out/'resident-rejected'
run(['singlecell-h5ad-import', source, '--plan', a.prepared/'mapping.json', '--output', resident], 'oversized sparse array lengths')
assert not resident.exists()
# No cell pseudoreplication: donor4 supplies only one biological replicate.
plan = json.loads((a.prepared/'stream-plan.json').read_text())
plan['contrasts'] = [dict(id='unreplicated-T2D', model='negativeBinomial', controlCondition='nonT2D',
    treatmentCondition='T2D', cellGroup='beta', design='independentReplicates', adjustForBatch=False,
    minimumReplicatesPerCondition=3)]
invalid_plan = a.out/'unreplicated-disease-plan.json'
invalid_plan.write_text(json.dumps(plan, indent=2)+'\n')
disease = a.out/'disease-rejected'
run(['singlecell-h5ad-pseudobulk', source, '--plan', invalid_plan, '--output', disease], 'insufficient biological replication')
assert not disease.exists()
summary = dict(status='passed', binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest(),
    checks=['full-source streamed publication', 'native reconstruction', 'exact repeated receipt',
            'resident count limit rejects without output', 'single T2D donor fails replicated NB DE without output'],
    commands=commands, qualification='Native count workflow and rejection checks; not integration or disease-effect evidence')
(a.out/'checks.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({key: summary[key] for key in ['status', 'checks', 'binarySHA256']}, indent=2))
