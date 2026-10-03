"""Virtual Wet Lab experiment records over NumiVivo's native donor predictor.

No prediction model is implemented here. The native owner fits, predicts and
replays source counts. This module selects eligible specimens, freezes plans,
and scores the held-out RNA measurements only after prediction is sealed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone

FORMAT = 'numivivo-virtual-wet-lab/v1'
LIMITS = ('Conditional donor-average RNA estimates only; no phenotype, survival, '
          'mechanistic or clinical prediction. Previously inspected public data: '
          'development replay, not a fresh prospective validation. Cells are '
          'subsamples; donors are biological units. Replays are not replicates.')


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def native(binary, args, log):
    result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
    write(log, {'argv': [str(binary), *map(str, args)], 'exitCode': result.returncode,
                'stdout': result.stdout, 'stderr': result.stderr})
    require(result.returncode == 0, 'Native NumiVivo failed; retained log: ' + str(log))
    return result


def inventory(root):
    files = sorted(p for p in root.rglob('*') if p.is_file())
    require(not any(p.is_symlink() for p in root.rglob('*')), 'Symlinks are not experiment artifacts')
    return {str(p.relative_to(root)): sha(p) for p in files}


def pairs(bulk, assay):
    """Expose only donor pairs with adequate source cells in both declared arms."""
    found = {}
    for i, g in enumerate(bulk['groups']):
        if g.get('cellGroup') != assay['cellGroup']:
            continue
        d, c = g.get('donorID'), g['condition']
        if not d or c not in (assay['controlCondition'], assay['treatmentCondition']):
            continue
        found.setdefault(d, {}).setdefault(c, []).append(i)
    eligible = {}
    for d, conditions in sorted(found.items()):
        if any(len(conditions.get(c, [])) != 1 for c in
               (assay['controlCondition'], assay['treatmentCondition'])):
            continue
        rows = [conditions[c][0] for c in (assay['controlCondition'], assay['treatmentCondition'])]
        gs = [bulk['groups'][r] for r in rows]
        if gs[0]['biologicalReplicateID'] != gs[1]['biologicalReplicateID']:
            continue
        if min(len(g['sourceCellIndices']) for g in gs) < assay['minimumCellsPerArm']:
            continue
        require(set(gs[0]['sourceCellIndices']).isdisjoint(gs[1]['sourceCellIndices']), 'Control/treatment cells overlap')
        eligible[d] = rows
    require(len(eligible) >= 3, 'Assay needs at least three eligible paired donors')
    units = [bulk['groups'][rows[0]]['biologicalReplicateID'] for rows in eligible.values()]
    require(len(set(units)) == len(units), 'Biological unit aliases span donors')
    return eligible


def load_assay(config):
    config = Path(config).resolve()
    assay = read(config)
    require(assay['format'] == FORMAT, 'Unsupported assay format')
    require(assay['family'] == 'cell-response', 'This release executes cell-response assays only')
    require(assay['observationHours'] > 0 and assay['minimumCellsPerArm'] >= 1, 'Invalid assay support')
    source = (config.parent / assay['sourceBundle']).resolve()
    require(sha(source / 'report.json') == assay['sourceReportSHA256'], 'Source report changed')
    receipt = read(source / 'receipt.json')
    require(bytes(receipt['report']['bytes']).hex() == assay['sourceReportSHA256'], 'Source receipt differs')
    require(read(source / 'plan.json')['mapping']['evidence'] == 'measured', 'Assay requires measured source data')
    bulk = read(source / 'report.json')['pseudobulk']
    return assay, source, bulk, pairs(bulk, assay)


def catalog(config):
    assay, _, bulk, eligible = load_assay(config)
    return {**{k: assay[k] for k in ('id', 'title', 'family', 'cellGroup', 'intervention',
             'observationHours', 'sourceCitation', 'sourceReportSHA256', 'provenance')},
            'limits': LIMITS, 'model': 'Native context ridge; fixed alpha=1',
            'featureCount': len(bulk['featureIDs']), 'specimens': [
                {'donor': d, 'controlCells': len(bulk['groups'][rows[0]]['sourceCellIndices']),
                 'treatedCells': len(bulk['groups'][rows[1]]['sourceCellIndices']),
                 'trainingDonors': len(eligible) - 1} for d, rows in eligible.items()]}


def runtime_identity(binary):
    binary = Path(binary).resolve()
    hdf5 = os.environ.get('NUMIVIVO_HDF5_LIBRARY')
    require(hdf5 and Path(hdf5).is_file(), 'Set NUMIVIVO_HDF5_LIBRARY to the native HDF5 library')
    return {'binarySHA256': sha(binary), 'hdf5SHA256': sha(hdf5),
            'adapterSHA256': sha(__file__), 'backend': 'native Swift CPU, log1p-CPM donor aggregates'}


def predict(config, binary, workspace, donor, hours, intervention):
    assay, source, bulk, eligible = load_assay(config)
    require(donor in eligible, 'Unsupported specimen')
    require(hours == assay['observationHours'], 'Unsupported observation time')
    require(intervention == assay['intervention'], 'Unsupported intervention')
    runtime = runtime_identity(binary)
    workspace = Path(workspace).resolve(); workspace.mkdir(parents=True, exist_ok=True)
    run = workspace / uuid.uuid4().hex; run.mkdir()
    training = [r for d, rows in eligible.items() if d != donor for r in rows]
    control, treated = eligible[donor]
    fold = {'id': 'held-out-donor', 'perturbationID': assay['id'],
            'controlCondition': assay['controlCondition'], 'treatmentCondition': assay['treatmentCondition'],
            'cellGroup': assay['cellGroup'], 'trainingGroupIndices': training, 'queryGroupIndices': [control]}
    plan = {'schemaVersion': 1, 'sourceReport': {'bytes': list(bytes.fromhex(assay['sourceReportSHA256']))},
            'featureNamespace': assay['featureNamespace'], 'provenance': assay['provenance'], 'folds': [fold]}
    registration = {'format': FORMAT, 'createdAt': timestamp(), 'assay': assay,
                    'donor': donor, 'trainingDonors': [d for d in eligible if d != donor],
                    'controlRow': control, 'heldOutRow': treated, 'runtime': runtime,
                    'randomness': 'none: deterministic fit; repeating a run adds no biological replicate',
                    'purpose': 'development replay', 'primaryModel': 'contextRidge',
                    'analysis': 'all-source-genes RNA-response RMSE and MAE in natural-log(1+CPM); '
                                'strict improvement over both noChange and meanResponse',
                    'stoppingRule': 'One selected held-out donor; retain every model and every gene.',
                    'limits': LIMITS}
    # Retain the exact adapter and preparation provenance alongside native source snapshots.
    shutil.copy2(__file__, run / 'wetlab-adapter.py')
    if assay.get('preparationSHA256'):
        preparation = Path(config).resolve().parent / 'preparation.json'
        require(sha(preparation) == assay['preparationSHA256'], 'Source preparation changed')
        shutil.copy2(preparation, run / 'preparation.json')
        original = preparation.parent / 'original-kang.h5ad'
        require(sha(original) == read(preparation)['originalSHA256'], 'Original source changed')
        shutil.copy2(original, run / 'original-source.h5ad')
    write(run / 'registration.json', registration)
    write(run / 'batch-plan.json', plan)
    try:
        native(binary, ['singlecell-perturbation-batch', source, '--plan', run / 'batch-plan.json',
                        '--output', run / 'prediction'], run / 'predict-log.json')
        receipt = read(run / 'prediction/receipt.json')
        require(all(f['status'] == 'completed' for f in receipt['folds']), 'Native fold failed; inspect retained receipt')
        require(read(run / 'prediction/folds/000/model.json')['trainingDonors'] == sorted(registration['trainingDonors']),
                'Native training donors differ from registration')
        # Seal prediction, raw source, model, plan, exact owner and logs BEFORE scoring.
        write(run / 'seal.json', {'format': FORMAT, 'sealedAt': timestamp(), 'files': inventory(run)})
    except Exception as error:
        write(run / 'failure.json', {'at': timestamp(), 'message': str(error)})
        raise
    return run


def check_seal(run, binary=None):
    run = Path(run)
    seal = read(run / 'seal.json')
    require(seal['format'] == FORMAT, 'Unsupported experiment seal')
    for name, digest in seal['files'].items():
        path = run / name
        require(path.resolve().is_relative_to(run.resolve()) and not path.is_symlink(), 'Invalid artifact path')
        require(sha(path) == digest, 'Sealed artifact changed: ' + name)
    registration = read(run / 'registration.json')
    if binary:
        require(runtime_identity(binary) == registration['runtime'], 'Replay needs the recorded executable, HDF5 and adapter')
    return registration


def score_values(run):
    registration = check_seal(run)
    bulk = read(run / 'prediction/source/report.json')['pseudobulk']
    report = read(run / 'prediction/folds/000/prediction.json')
    require(report['featureIDs'] == bulk['featureIDs'], 'Prediction and observation gene axes differ')
    require(len(report['predictions']) == 1, 'Expected one held-out query')
    p = report['predictions'][0]; row = registration['heldOutRow']; group = bulk['groups'][row]
    assay = registration['assay']
    require(group['donorID'] == registration['donor'] and group['condition'] == assay['treatmentCondition']
            and group['cellGroup'] == assay['cellGroup'], 'Held-out specimen identity differs')
    require(p['group'] == bulk['groups'][registration['controlRow']], 'Control identity differs')
    m = bulk['matrix']; counts = [0] * len(bulk['featureIDs'])
    for k in range(m['rowOffsets'][row], m['rowOffsets'][row + 1]):
        counts[m['featureIndices'][k]] = m['counts'][k]
    total = sum(counts); require(total > 0, 'Empty held-out library')
    observed = [math.log1p(x / total * 1e6) for x in counts]
    truth = [o - c for o, c in zip(observed, p['control'])]
    metrics = []
    for e in p['estimates']:
        require(len(e['predictedResponse']) == len(truth), 'Incomplete prediction vector')
        errors = [a - b for a, b in zip(e['predictedResponse'], truth)]
        metrics.append({'model': e['baseline'], 'rmse': math.sqrt(math.fsum(x*x for x in errors)/len(errors)),
                        'mae': math.fsum(abs(x) for x in errors)/len(errors)})
    by_model = {m['model']: m['rmse'] for m in metrics}
    advantage = all(by_model['contextRidge'] < by_model[b] for b in ('noChange', 'meanResponse'))
    return {'format': FORMAT, 'donor': registration['donor'], 'units': 'natural-log(1+CPM)',
            'features': len(truth), 'metrics': metrics, 'observedTreated': observed,
            'observedResponse': truth, 'verdict': 'beats-both-baselines' if advantage else 'does-not-beat-both-baselines',
            'limits': LIMITS}


def reveal(run, binary):
    run = Path(run); check_seal(run, binary)
    require(not (run / 'comparison.json').exists(), 'Measurements were already revealed; use verify')
    # Native reconstruction checks every original count and the donor split.
    log = run / ('reveal-verify-' + uuid.uuid4().hex + '.json')
    native(binary, ['singlecell-perturbation-batch-verify', run / 'prediction'], log)
    result = score_values(run)
    result['predictionSealSHA256'] = sha(run / 'seal.json')
    write(run / 'comparison.json', result)
    return result


def verify(run, binary):
    run = Path(run); check_seal(run, binary)
    native(binary, ['singlecell-perturbation-batch-verify', run / 'prediction'],
           run / ('replay-' + uuid.uuid4().hex + '.json'))
    if (run / 'comparison.json').exists():
        expected = score_values(run); expected['predictionSealSHA256'] = sha(run / 'seal.json')
        require(read(run / 'comparison.json') == expected, 'Held-out comparison does not reconstruct')
    return {'status': 'verified', 'run': run.name, 'revealed': (run / 'comparison.json').exists()}


def summary(run):
    run = Path(run)
    registration = check_seal(run)
    report = read(run / 'prediction/folds/000/prediction.json')
    result = {'id': run.name, 'recordDirectory': str(run.resolve()), 'registration': registration, 'featureIDs': report['featureIDs'],
              'prediction': report['predictions'][0], 'revealed': (run / 'comparison.json').exists()}
    if result['revealed']:
        comparison = read(run / 'comparison.json')
        expected = score_values(run); expected['predictionSealSHA256'] = sha(run / 'seal.json')
        require(comparison == expected, 'Comparison changed')
        result['comparison'] = comparison
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('catalog'); p.add_argument('config', type=Path)
    p = sub.add_parser('predict'); p.add_argument('config', type=Path)
    p.add_argument('--binary', type=Path, required=True); p.add_argument('--workspace', type=Path, required=True)
    p.add_argument('--donor', required=True); p.add_argument('--hours', type=float, required=True)
    p.add_argument('--intervention', required=True)
    for name in ('reveal', 'verify', 'inspect'):
        p = sub.add_parser(name); p.add_argument('run', type=Path)
        if name != 'inspect': p.add_argument('--binary', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'catalog': result = catalog(args.config)
    elif args.command == 'predict': result = {'run': str(predict(args.config, args.binary, args.workspace, args.donor, args.hours, args.intervention))}
    elif args.command == 'reveal': result = reveal(args.run, args.binary)
    elif args.command == 'verify': result = verify(args.run, args.binary)
    else: result = summary(args.run)
    print(json.dumps(result, allow_nan=False))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        print('Virtual Wet Lab: ' + str(error), file=sys.stderr)
        sys.exit(1)
