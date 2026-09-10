#!/usr/bin/env python3
"""Run every Norman held-target native fit using a bounded sequential disk workspace.

Standard-library driver; no fitting or response prediction occurs in Python.
Full training aggregates are prepared once from a freshly verified source. Each
native fitter receives only control and the other 104 single-target rows.
"""
import argparse
import copy
import gzip
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

SOURCE = 'efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc'
REPORT = 'c36a13a80371e80fc5c0eb0dee598cdaed69dc955c7df312a4174588bffcee9c'
ANNOTATIONS = 'dc2b8bcc73d4ea362b2df803b99f8601729784a0865004b27c878763ffef21e4'
COVERAGE = '9c9b7c510a356941204ed83cd4f9ed48214d116ed810b7d63aa4d1019d03835c'
CONTEXT = dict(id='Norman2019-K562-pooled', organism='NCBITaxon:9606', featureNamespace='Ensembl-gene',
               perturbationNamespace='Norman2019-guide-target', countUnit='umiCount')
NAMESPACE = 'GO-direct-BP-MF-CC-MyGene-20260906'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''): digest.update(block)
    return digest.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')


def compress(source, destination):
    with source.open('rb') as raw, destination.open('wb') as output:
        with gzip.GzipFile(filename='', fileobj=output, mode='wb', mtime=0) as compressed:
            shutil.copyfileobj(raw, compressed)


def subset(training, keep):
    original = training['matrix']; offsets, columns, values = [0], [], []
    for row in [0, *[i+1 for i in keep]]:
        start, end = original['rowOffsets'][row:row+2]
        columns.extend(original['featureIndices'][start:end]); values.extend(original['counts'][start:end]); offsets.append(len(values))
    return dict(training, selection=dict(training['selection'], targets=[training['selection']['targets'][i] for i in keep]),
                matrix=dict(original, cellCount=len(keep)+1, rowOffsets=offsets, featureIndices=columns, counts=values))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('binary', 'pseudobulk', 'annotations', 'coverage', 'out'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    assert sha(a.annotations) == ANNOTATIONS and sha(a.coverage) == COVERAGE
    assert sha(a.pseudobulk/'report.json') == REPORT
    a.out.mkdir(parents=True, exist_ok=False)
    annotations = json.loads(a.annotations.read_text()); coverage = json.loads(a.coverage.read_text())
    targets = [c['target'] for c in coverage]
    assert len(targets) == 105 and targets == sorted(targets)
    descriptor = {}
    for c in coverage:
        target = c['target']
        status = {'supported': 'available', 'no-usable-GO-annotations': 'noData', 'unresolved-source-symbol': 'unresolvedIdentity'}[c['status']]
        descriptor[target] = dict(targetID=target, status=status, terms=annotations.get(target, {}).get('terms', []))
        if 'sourceFeatureID' in c: descriptor[target]['featureID'] = c['sourceFeatureID']
    selection = dict(schemaVersion=1, context=CONTEXT, controlCondition='control',
                     targets=[dict(id=t, condition=t) for t in targets],
                     provenance='Norman 2019 full filtered release; fresh native source reconstruction. Pooled K562, no donor replication claim. Each target-kernel fold removes its held target before fitting; all paired conditions excluded.')
    write(a.out/'selection.json', selection)
    commands = []; started = time.monotonic()
    def run(name, command, expected=0):
        before = time.monotonic()
        with (a.out/(name+'.log')).open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-l', str(a.binary.resolve()), *map(str, command)], stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(name=name, returnCode=result.returncode, expected=expected, elapsedSeconds=time.monotonic()-before))
        write(a.out/'commands.json', commands)
        assert result.returncode == expected, (name, result.returncode)
    run('prepare', ['singlecell-composition-prepare', a.pseudobulk, '--plan', a.out/'selection.json', '--output', a.out/'training.json'])
    training = json.loads((a.out/'training.json').read_text())
    assert training['selection'] == selection and training['evidence'] == 'measured'
    assert bytes(training['source']['bytes']).hex() == SOURCE and bytes(training['sourceReport']['bytes']).hex() == REPORT
    assert training['matrix']['cellCount'] == 106 and training['matrix']['featureCount'] == 33694
    compress(a.out/'training.json', a.out/'training.json.gz')
    records = []
    for held, target in enumerate(targets):
        fold = a.out/'folds'/target; fold.mkdir(parents=True)
        scratch = a.out/'scratch'; scratch.mkdir()
        keep = [i for i in range(105) if i != held]
        selected = subset(training, keep)
        assert target not in [t['id'] for t in selected['selection']['targets']]
        # Selection is insensitive to the held row; replace all its entries and select again.
        mutated = dict(training, matrix=dict(training['matrix'], counts=training['matrix']['counts'].copy()))
        start, end = mutated['matrix']['rowOffsets'][held+1:held+3]
        mutated['matrix']['counts'][start:end] = [1]*(end-start)
        assert subset(mutated, keep) == selected
        del mutated
        write(scratch/'training.json', selected); del selected
        plan = dict(schemaVersion=1, context=CONTEXT, descriptorNamespace=NAMESPACE, source=dict(bytes=list(bytes.fromhex(ANNOTATIONS))),
                    provenance='Pinned direct GO terms from exact original Ensembl identities; MyGene build 20260906. Fixed lambda=1, no outer-outcome selection.',
                    targets=[descriptor[targets[i]] for i in keep], regularization=1, maximumWork=200000000)
        query = dict(schemaVersion=1, context=CONTEXT, descriptorNamespace=NAMESPACE, source=plan['source'],
                     queries=[dict(id=target, descriptor=descriptor[target])], maximumWork=200000000)
        write(fold/'plan.json', plan); write(fold/'query-input.json', query)
        run(target+'-fit', ['singlecell-target-kernel-fit', scratch/'training.json', '--plan', fold/'plan.json', '--output', scratch/'model'])
        run(target+'-verify', ['singlecell-target-kernel-verify', scratch/'model'])
        run(target+'-predict', ['singlecell-target-kernel-predict', scratch/'model', '--plan', fold/'query-input.json', '--output', scratch/'prediction'])
        run(target+'-prediction-verify', ['singlecell-target-kernel-prediction-verify', scratch/'prediction'])
        files = []
        for name in ('training.json', 'plan.json', 'model.json', 'receipt.json'):
            path = scratch/'model'/name
            assert sha(path) == sha(scratch/'prediction/reference'/name)
            files.append(dict(path='model/'+name, bytes=path.stat().st_size, sha256=sha(path)))
        for name in ('query.json', 'report.json', 'receipt.json'):
            path = scratch/'prediction'/name
            files.append(dict(path='prediction/'+name, bytes=path.stat().st_size, sha256=sha(path)))
            if name == 'report.json': compress(path, fold/(name+'.gz'))
            else: shutil.copy2(path, fold/name)
        shutil.copy2(scratch/'model/receipt.json', fold/'model-receipt.json')
        result = json.loads((scratch/'prediction/report.json').read_text())
        record = dict(target=target, status=result['queries'][0]['status'], trainingTargets=[targets[i] for i in keep],
                      supportedTrainingTargets=result['supportedTrainingTargets'], files=files, heldOutcomeMutationSelectedCountsExact=True)
        write(fold/'manifest.json', record); records.append(record)
        if held == 0:
            # Retain one complete reference bundle; all other models are reconstructible from full training plus fold selection.
            retained = a.out/'first-model'; retained.mkdir()
            for path in (scratch/'model').iterdir(): compress(path, retained/(path.name+'.gz'))
            run('first-repeat', ['singlecell-target-kernel-predict', scratch/'model', '--plan', fold/'query-input.json', '--output', scratch/'repeat'])
            assert sha(scratch/'repeat/report.json') == sha(scratch/'prediction/report.json')
            assert sha(scratch/'repeat/receipt.json') == sha(scratch/'prediction/receipt.json')
            reordered = subset(training, list(reversed(keep)))
            write(scratch/'reordered.json', reordered); del reordered
            run('reordered-fit', ['singlecell-target-kernel-fit', scratch/'reordered.json', '--plan', fold/'plan.json', '--output', scratch/'reordered'])
            assert sha(scratch/'reordered/model.json') == sha(scratch/'model/model.json')
        shutil.rmtree(scratch)
        write(a.out/'folds.json', records)
        print(f'{held+1}/105 {target}: {record["status"]}', flush=True)
    assert sum(r['status']=='predicted' for r in records) == 101
    assert not list(a.out.rglob('.numivivo-target-kernel-*'))
    write(a.out/'receipt.json', dict(status='passed', targets=105, supported=101, features=33694, commands=len(commands),
          binarySHA256=sha(a.binary), driverSHA256=sha(Path(__file__)), annotationsSHA256=ANNOTATIONS, coverageSHA256=COVERAGE,
          sourceSHA256=SOURCE, sourceReportSHA256=REPORT, trainingSHA256=sha(a.out/'training.json'),
          elapsedSeconds=time.monotonic()-started, firstFoldReplayExact=True, firstFoldReorderedTrainingModelExact=True,
          retention='All prediction reports, queries, receipts and model hashes; one full model. Other models reconstruct from full training and the per-fold selected target list.'))

if __name__ == '__main__': main()
