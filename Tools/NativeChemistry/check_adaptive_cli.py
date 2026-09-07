#!/usr/bin/env python3
"""Execute the real native campaign CLI. Python supplies no production numerics."""
import copy
import json
import pathlib
import subprocess
import sys

binary = pathlib.Path(sys.argv[1]).resolve(strict=True)
root = pathlib.Path(sys.argv[2]).resolve()
root.mkdir(parents=True, exist_ok=True)
checks = []
commands = []


def run(*arguments, expected=0):
    command = [str(binary), *map(str, arguments)]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    commands.append(dict(command=command, returncode=result.returncode,
                         stdout=result.stdout, stderr=result.stderr))
    (root / 'commands.json').write_text(json.dumps(commands, indent=2) + '\n')
    if result.returncode != expected:
        raise AssertionError((command, expected, result.returncode, result.stdout, result.stderr))
    return result


for name in ('equilibrium-correction', 'electronic-crosscheck'):
    request = root / (name + '.json')
    plan = root / (name + '.plan.json')
    output = root / (name + '.result.json')
    store = root / (name + '-store')
    run('campaign-template', name, '--output', request)
    run('campaign-plan', request, '--output', plan)
    assert len(json.loads(plan.read_text())['actionPlans']) == 2
    run('campaign-run', request, '--store', store, '--output', output)
    original = json.loads(output.read_text())
    assert original['termination'] == 'criteriaSatisfied'
    assert len(original['steps']) == 2
    assert all(s['record']['metric']['nativeChecksPassed'] for s in original['steps'])
    resumed_path = root / (name + '.resumed.json')
    run('campaign-resume', original['checkpointArtifact'], '--store', store, '--output', resumed_path)
    resumed = json.loads(resumed_path.read_text())
    assert resumed['termination'] == 'criteriaSatisfied'
    assert resumed['steps'] == original['steps']
    repeated_path = root / (name + '.cached.json')
    run('campaign-run', request, '--store', store, '--output', repeated_path)
    repeated = json.loads(repeated_path.read_text())
    assert repeated['termination'] == 'criteriaSatisfied'
    assert all(not s['record']['entirelyUncachedExecution'] for s in repeated['steps'])
    assert [s['record']['metric'] for s in repeated['steps']] == [s['record']['metric'] for s in original['steps']]
    run('campaign-resume', original['checkpointArtifact'], '--store', store, '--output', resumed_path, expected=1)
    checks.append(name + ': native run, checkpoint replay, cache validation and output protection')

base = json.loads((root / 'equilibrium-correction.json').read_text())
for name, change, status, termination in (
    ('admission-exhausted', lambda r: r.update(maximumDeclaredWorkUnits=1), 2, 'declaredWorkLimit'),
    ('wrong-observable', lambda r: r['criteria'][0].update(observableFingerprint='f' * 64), 1, None),
    ('unknown-schema', lambda r: r.update(schema='unknown-required-schema'), 1, None),
):
    value = copy.deepcopy(base)
    change(value)
    request = root / (name + '.json')
    output = root / (name + '.result.json')
    request.write_text(json.dumps(value))
    run('campaign-run', request, '--store', root / 'equilibrium-correction-store', '--output', output, expected=status)
    if termination:
        result = json.loads(output.read_text())
        assert result['termination'] == termination and result['steps'] == []
        assert pathlib.Path(str(output) + '.receipt.json').exists()
    checks.append(name + ': expected native failure status')

report = dict(schema='numivivo.org/adaptive-cli-observations/v1', status='passed',
              binary=str(binary), checks=checks, commandCount=len(commands),
              scope='Real native CLI, workflow and artifact execution on finite-state/algebraic fixtures; not molecular accuracy or GPU speed validation')
(root / 'results.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
