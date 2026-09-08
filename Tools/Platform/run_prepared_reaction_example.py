#!/usr/bin/env python3
"""Run the published prepared H2 -> sampled geometry -> H3 -> encounter example.

Uses public CLI processes and existing numerical engines. The harmonic sampling
model is synthetic; its finite prefix supplies initial geometry only. The H3
rate and maintained H2 bath remain explicit, conditional model assumptions.
Run on an Apple host after coordinating GPU ownership. Keep the output folder.
"""
import argparse
import copy
import hashlib
import json
import math
import pathlib
import shutil
import subprocess
import time


def digest(value):
    return bytes(value['bytes']).hex()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=pathlib.Path, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    args = parser.parse_args()
    binary, out = args.binary.resolve(strict=True), args.out.resolve()
    example = pathlib.Path(__file__).resolve().parents[2] / 'Examples/prepared-reaction'
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; previous evidence is never deleted')
    out.mkdir(parents=True, exist_ok=True)
    store = out / 'store'
    for name in ['preparation.json', 'campaign.json']:
        shutil.copyfile(example / name, out / name)
    config = json.loads((out / 'campaign.json').read_text())
    if config['schema'] != 'numivivo.org/examples/prepared-reaction/v1':
        raise RuntimeError('Unsupported published campaign schema')
    report = {'schema': 'numivivo.org/test-evidence/prepared-reaction-cli/v1',
              'binarySHA256': sha(binary), 'inputs': {n: sha(out / n) for n in ['preparation.json', 'campaign.json']},
              'scope': 'Synthetic prepared harmonic H2 geometry seed to freshly qualified H3 exchange and conditional first-event probability; no ensemble, experimental or performance claim',
              'checks': [], 'passed': False}

    def save():
        (out / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def check(condition, label):
        report['checks'].append({'label': label, 'passed': bool(condition)})
        save()
        if not condition:
            raise RuntimeError(label)

    def run(label, *arguments, expected=0):
        print(label, flush=True)
        command = [str(binary), *map(str, arguments)]
        started = time.monotonic()
        done = subprocess.run(command, capture_output=True, text=True, timeout=1800)
        (out / (label + '.stdout')).write_text(done.stdout)
        (out / (label + '.stderr')).write_text(done.stderr)
        report['checks'].append({'label': label, 'command': command, 'exitCode': done.returncode,
                                 'expectedExit': expected, 'elapsedSeconds': time.monotonic() - started,
                                 'passed': done.returncode == expected})
        save()
        if done.returncode != expected:
            raise RuntimeError(label + ': unexpected exit; inspect retained stdout/stderr')
        return done

    def read(name):
        return json.loads((out / name).read_text())

    def write(name, value):
        path = out / name
        with path.open('x') as handle:
            handle.write(json.dumps(value, indent=2, allow_nan=False) + '\n')
        return path

    def output(report_name, name):
        return next(x['artifact'] for x in read(report_name)['exports'] if x['name'] == name)

    def recipe(identifier, operation, inputs, budget):
        return {'schema': 'numivivo.org/workflow-recipe/v1', 'identifier': identifier,
                'artifacts': [{'identifier': key, 'source': {'stored': value}} for key, value in inputs.items()],
                'nodes': [{'identifier': identifier, 'operation': operation, 'version': '1',
                           'inputs': {key: {'artifact': {'identifier': key}} for key in inputs}, 'configuration': {},
                           'resources': {'numericalBackend': 'cpu-fp64', 'budget': budget,
                                         'maximumInputBytes': budget['maximumBytes'], 'maximumOutputBytes': budget['maximumBytes']}}],
                'outputs': [{'name': 'result', 'node': identifier, 'port': 'result'}],
                'policy': {'maximumConcurrentTasks': 1, 'maximumConcurrentMetalTasks': 1,
                           'maximumReservedBytes': 2147483648, 'maximumNodes': 256, 'maximumInlineBytes': 134217728}}

    help_text = run('seed-help', 'reaction-seed-help').stdout
    check(all(x in help_text for x in ['reaction-seed-verify', 'reaction-encounter-request']), 'new commands discoverable')
    catalog = json.loads(run('catalog', 'workflow-catalog').stdout)
    check('vivo.native.conditional-encounter' in {x['identifier'] for x in catalog}, 'encounter operation registered')
    missing = read('preparation.json')
    missing['forceField']['bondParameters'] = []
    run('missing-parameter', 'molecule-prepare', write('missing-parameter.json', missing), '--store', out / 'missing-store',
        '--output', out / 'missing-result.json', expected=65)
    check(not (out / 'missing-result.json').exists(), 'missing parameters cannot publish a successful preparation')
    run('prepare', 'molecule-prepare', out / 'preparation.json', '--store', store, '--output', out / 'prepared.json')
    run('sampling-template', 'molecule-sampling-template', out / 'prepared.json', '--context', 'published synthetic harmonic H2 geometry seed',
        '--distance', '0,1', '--temperature', '300', '--cutoff', '0.5', '--output', out / 'sampling-template.json')
    sampled = read('sampling-template.json')
    cfg = config['sampling']
    sampled['replicaSeeds'] = cfg['replicaSeeds']
    initial = sampled['initialStates'][0]
    sampled['initialStates'] = [dict(copy.deepcopy(initial), sourceTimePS=t) for t in cfg['sourceTimesPS']]
    sampled.pop('minimization', None)
    sampled['md'].update(timeStepPS=cfg['timeStepPS'], targetTemperatureK=cfg['temperatureK'],
                         frictionPerPS=cfg['frictionPerPS'], cutoffNM=cfg['cutoffNM'], neighborSkinNM=0,
                         ensemble='nvt', thermostat='langevinMiddle', neighborListEnabled=False, electrostatics='cutoff')
    for key in ['equilibrationSteps', 'stepsPerBlock', 'sampleEvery', 'maximumBlocks', 'requiredConsecutivePasses', 'trajectoryChunkBytes']:
        sampled[key] = cfg[key]
    for key in ['minimumRetainedFramesPerReplica', 'minimumReplicas', 'maximumRHat', 'minimumEffectiveSamplesPerReplica', 'maximumAutocorrelationLag']:
        sampled['convergence'][key] = cfg[key]
    sampled['observables'] = [sampled['observables'][0]]
    sampled['observables'][0]['maximumMeanStandardError'] = cfg['maximumMeanStandardError']
    run('sample', 'molecule-sampling-run', write('sampling-request.json', sampled), '--store', store,
        '--output', out / 'sampling-receipt.json', expected=75)
    sampling = read('sampling-receipt.json')
    check(sampling['status'] == 'budgetExhausted' and sampling['completedBlocks'] == cfg['maximumBlocks'], 'partial accepted sampling is explicit')
    imported = json.loads(run('import-sampling-receipt', 'workflow-import', out / 'sampling-receipt.json',
                             '--kind', 'molecular-sampling-receipt', '--store', store).stdout)
    run('export-sampling', 'molecule-sampling-export', '--store', store, '--checkpoint', digest(sampling['checkpoint']),
        '--replica', 0, '--receipt', digest(imported['fingerprint']), '--verify-all-payloads', '--output', out / 'sampling-export.json')
    run('destination-template', 'reaction-template', config['reactionTemplate'], '--output', out / 'destination.json')
    seed = {'schema': 'numivivo.org/sampling-reaction-seed-request/v1', 'sourceExport': read('sampling-export.json'),
            'destination': read('destination.json'), 'modelTransferStatement': config['modelTransferStatement'],
            'assignments': [{'target': {'connectedEndpoint': {'identifier': name, 'componentIndex': component}},
                             'sourceAtomByNucleus': [0, 1]} for name, component in [('H0-H1_plus_H2', 0), ('H0_plus_H1-H2', 1)]]}
    seed_path = write('seed-request.json', seed)
    run('seed', 'reaction-seed', seed_path, '--store', store, '--output', out / 'seed-receipt.json')
    run('verify-seed', 'reaction-seed-verify', out / 'seed-receipt.json', '--store', store,
        '--verify-all-payloads', '--output', out / 'seed-verification.json')
    run('seed-repeated', 'reaction-seed', seed_path, '--store', store, '--output', out / 'seed-repeated.json')
    check(read('seed-receipt.json') == read('seed-repeated.json'), 'repeated seed publication retains exact immutable identities')
    bad_seed = copy.deepcopy(seed)
    bad_seed['assignments'][0]['sourceAtomByNucleus'] = [0, 0]
    run('duplicate-mapping', 'reaction-seed', write('duplicate-mapping.json', bad_seed), '--store', store, expected=65)
    run('seed-no-clobber', 'reaction-seed', seed_path, '--store', store, '--output', out / 'seed-receipt.json', expected=65)
    run('seed-store-alias', 'reaction-seed', seed_path, '--store', store, '--output', store / 'forbidden.json', expected=65)
    seed_receipt = read('seed-receipt.json')
    check(not seed_receipt['provenance']['source']['declaredObservableCriteriaSatisfied'], 'seed preserves unmet sampling criteria')
    transfers = seed_receipt['provenance']['transfers']
    check(len(transfers) == 2 and all(x['sourceMassesDa'] == [1, 1] and x['destinationMassesDa'] == [1.008, 1.008]
          and x['destinationThermochemistry']['temperatureK'] == 298.15 for x in transfers), 'mass and temperature transfer contexts retained')
    request_output = next(x for x in seed_receipt['outputs'] if x['name'] == 'request')
    budget = seed['destination']['calculation']['connectedReaction']['request']['saddle']['model']['budget']
    reaction_recipe = recipe('reaction', 'vivo.native.reaction-qualification',
                             {'request': {'kind': request_output['kind'], 'fingerprint': request_output['artifact']}}, budget)
    run('reaction', 'workflow-run', write('reaction-workflow.json', reaction_recipe), '--store', store, '--output', out / 'reaction-report.json')
    check(read('reaction-report.json')['allTasksSucceeded'], 'fresh connected reaction workflow succeeded')
    check(not read('reaction-report.json')['nodes'][0]['outcome']['succeeded']['reused'], 'reaction qualification was freshly executed')
    reaction_output = output('reaction-report.json', 'result')
    run('export-reaction', 'workflow-export', digest(reaction_output['artifact']), '--kind', reaction_output['kind'],
        '--store', store, '--output', out / 'reaction-result.json')
    tst = read('reaction-result.json')['connectedReaction']['result']
    check(tst['molecularity'] == 2 and tst['rateUnits'] == 'Pa^-1 s^-1', 'separated H plus H2 retains bimolecular pressure-rate units')
    conditions = {'schema': 'numivivo.org/conditional-encounter-conditions/v1', 'identifier': 'tagged-H-in-maintained-H2',
                  'reactantEndpointIdentifier': 'H0-H1_plus_H2', 'temperatureK': 298.15, 'taggedComponentIndex': 1,
                  'reservoirs': [{'componentIndex': 0, 'value': config['encounter']['reservoirPressurePa'], 'unit': 'Pa'}],
                  'observationTimesSeconds': config['encounter']['observationTimesSeconds'], 'includeReservoirSensitivities': True}
    run('encounter-request', 'reaction-encounter-request', out / 'reaction-result.json', '--conditions', write('conditions.json', conditions),
        '--output', out / 'encounter-request.json')
    for label, update in [('wrong-temperature', {'temperatureK': 300}), ('missing-reservoir', {'reservoirs': []}),
                          ('wrong-units', {'reservoirs': [{'componentIndex': 0, 'value': 1, 'unit': 'mol/L'}]})]:
        invalid = dict(copy.deepcopy(conditions), **update)
        run(label, 'reaction-encounter-request', out / 'reaction-result.json', '--conditions', write(label + '.json', invalid), expected=65)
    imported = json.loads(run('import-encounter-request', 'workflow-import', out / 'encounter-request.json',
                             '--kind', 'vivo.conditional-encounter-request', '--store', store).stdout)
    encounter_recipe = recipe('encounter', 'vivo.native.conditional-encounter',
                              {'request': {'kind': imported['kind'], 'fingerprint': imported['fingerprint']},
                               'reaction': {'kind': reaction_output['kind'], 'fingerprint': reaction_output['artifact']}}, budget)
    run('encounter', 'workflow-run', write('encounter-workflow.json', encounter_recipe), '--store', store, '--output', out / 'encounter-report.json')
    check(read('encounter-report.json')['allTasksSucceeded'], 'conditional encounter workflow succeeded')
    encounter_output = output('encounter-report.json', 'result')
    run('export-encounter', 'workflow-export', digest(encounter_output['artifact']), '--kind', encounter_output['kind'],
        '--store', store, '--output', out / 'encounter-result.json')
    result = read('encounter-result.json')['encounter']
    hazard = result['conditionalHazardPerSecond']
    check(math.isclose(hazard, tst['rateConstant'] * config['encounter']['reservoirPressurePa'], rel_tol=1e-12), 'rate and bath pressure produce inverse-second hazard')
    for point in result['kinetics']['observations']:
        expected = -math.expm1(-hazard * point['timeSeconds'])
        check(abs(point['reactedProbability'] - expected) <= 5e-10, 'first-event probability at ' + str(point['timeSeconds']) + ' seconds')
    check(result['transmission'] == tst['request']['transmission'] and result['transmission']['kind'] == 'assumed', 'assumed transmission remains unchanged')
    report['reactionArtifact'] = digest(reaction_output['artifact'])
    report['encounterArtifact'] = digest(encounter_output['artifact'])
    report['observations'] = result['kinetics']['observations']
    report['rateConstant'] = tst['rateConstant']
    report['rateUnits'] = tst['rateUnits']
    report['conditionalHazardPerSecond'] = hazard
    report['passed'] = all(x['passed'] for x in report['checks'])
    save()
    print(json.dumps({'passed': report['passed'], 'checks': len(report['checks']), 'output': str(out)}), flush=True)


if __name__ == '__main__':
    main()
