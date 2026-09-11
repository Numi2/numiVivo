#!/usr/bin/env python3
"""Freeze all 150 inputs, then execute the unchanged native predictor and replay.

The driver selects counts and orchestrates native commands; it does not fit,
predict or score in Python. Original gemgroup IDs remain unnamed technical groups.
"""
import argparse
import datetime
import gzip
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

from prepare_training import BINARY, encoded, sha, subset, write


def compress(source, destination):
    with source.open('rb') as raw, destination.open('wb') as stream:
        with gzip.GzipFile(filename='', fileobj=stream, mode='wb', mtime=0) as zipped:
            shutil.copyfileobj(raw, zipped)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--descriptors', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root, descriptors, out = args.root.resolve(), args.descriptors.resolve(), args.out.resolve()
    preparation = root/'training-preparation'
    binary = root/'runtime/numivivo-omics'
    assert sha(binary) == BINARY
    checks = json.loads((preparation/'checks.json').read_text())
    assert checks['status'] == 'passed-count-preparation-only'
    frozen_counts = json.loads((preparation/'count-folds.json').read_text())
    receipt = json.loads((descriptors/'receipt.json').read_text())
    assert receipt['status'] == 'captured-and-frozen' and receipt['targets'] == 30
    for item in receipt['files']:
        assert sha(descriptors/item['path']) == item['SHA256']
    annotation_sha = sha(descriptors/'annotations.json')
    annotations = json.loads((descriptors/'annotations.json').read_text())
    coverage = {x['guide']: x for x in json.loads((descriptors/'coverage.json').read_text())}
    assert len(coverage) == 30
    namespace = 'GO-direct-BP-MF-CC-MyGene-'+str(receipt['buildVersion'])
    descriptor = {t: dict(targetID=t, featureID=c['sourceFeatureID'],
        status='available' if t in annotations else 'noData', terms=annotations.get(t, {}).get('terms', []))
        for t, c in coverage.items()}
    out.mkdir(parents=True, exist_ok=False)
    freeze, training_by_context = [], {}
    for context in frozen_counts['contexts']:
        gem = context['gemgroup']
        path = preparation/(gem+'-training.json')
        assert sha(path) == context['trainingSHA256']
        training = json.loads(path.read_text())
        assert training['matrix']['cellCount'] == 31 and training['matrix']['featureCount'] == 33694
        assert set(x['id'] for x in training['selection']['targets']) == set(coverage)
        training_by_context[gem] = training
        compress(path, out/(gem+'-training.json.gz'))
    for count_fold in frozen_counts['folds']:
        gem, target = count_fold['gemgroup'], count_fold['heldGuide']
        training = training_by_context[gem]
        keep = [i for i, t in enumerate(training['selection']['targets']) if t['id'] != target]
        selected = subset(training, keep)
        assert hashlib.sha256(encoded(selected)).hexdigest() == count_fold['selectedTrainingSHA256']
        targets = [t['id'] for t in selected['selection']['targets']]
        assert targets == count_fold['trainingGuides'] and len(targets) == 29
        fold = out/'folds'/gem/target
        fold.mkdir(parents=True)
        plan = dict(schemaVersion=1, context=training['selection']['context'], descriptorNamespace=namespace,
            source=dict(bytes=list(bytes.fromhex(annotation_sha))), regularization=1, maximumWork=200000000,
            provenance='Fixed Replogle protocol: original guide labels; FBA reproduced guide/gene panel and exact reference core support intended gene identities. Direct GO, lambda1, no outcome tuning. Technical gemgroup only.',
            targets=[descriptor[t] for t in targets])
        query = dict(schemaVersion=1, context=plan['context'], descriptorNamespace=namespace, source=plan['source'],
            maximumWork=200000000, queries=[dict(id=target, descriptor=descriptor[target])])
        write(fold/'plan.json', plan)
        write(fold/'query-input.json', query)
        freeze.append(dict(gemgroup=gem, target=target, keep=keep,
            selectedTrainingSHA256=count_fold['selectedTrainingSHA256'], planSHA256=sha(fold/'plan.json'),
            querySHA256=sha(fold/'query-input.json')))
    assert len(freeze) == 150
    write(out/'input-freeze.json', dict(frozenUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        binarySHA256=BINARY, driverSHA256=sha(Path(__file__)), countSelectorSHA256=sha(Path(__file__).with_name('prepare_training.py')),
        countPreparationChecksSHA256=sha(preparation/'checks.json'), countFoldsSHA256=sha(preparation/'count-folds.json'),
        descriptorReceiptSHA256=sha(descriptors/'receipt.json'), annotationsSHA256=annotation_sha,
        sourceIdentitySHA256=receipt['identitiesSHA256'], folds=freeze, scoringStarted=False, fittingStarted=False))
    commands, results = [], []

    def run(label, command):
        start = time.monotonic()
        with (out/(label+'.log')).open('wb') as log:
            result = subprocess.run(['/usr/bin/time', '-l', str(binary), *map(str, command)], stdout=log, stderr=subprocess.STDOUT)
        commands.append(dict(label=label, command=list(map(str, command)), returnCode=result.returncode, seconds=time.monotonic()-start))
        write(out/'commands.json', commands)
        assert result.returncode == 0, label

    for i, item in enumerate(freeze):
        gem, target = item['gemgroup'], item['target']
        fold, scratch = out/'folds'/gem/target, out/'scratch'
        scratch.mkdir()
        selected = subset(training_by_context[gem], item['keep'])
        write(scratch/'training.json', selected)
        assert sha(scratch/'training.json') == item['selectedTrainingSHA256']
        assert sha(fold/'plan.json') == item['planSHA256'] and sha(fold/'query-input.json') == item['querySHA256']
        stem = gem+'-'+target
        for label, command in [
            ('fit', ['singlecell-target-kernel-fit', scratch/'training.json', '--plan', fold/'plan.json', '--output', scratch/'model']),
            ('model-verify', ['singlecell-target-kernel-verify', scratch/'model']),
            ('predict', ['singlecell-target-kernel-predict', scratch/'model', '--plan', fold/'query-input.json', '--output', scratch/'prediction']),
            ('prediction-verify', ['singlecell-target-kernel-prediction-verify', scratch/'prediction'])]:
            run(stem+'-'+label, command)
        files = []
        for directory, names in [('model', ['training.json', 'plan.json', 'model.json', 'receipt.json']),
                                 ('prediction', ['query.json', 'report.json', 'receipt.json'])]:
            for name in names:
                path = scratch/directory/name
                files.append(dict(path=directory+'/'+name, bytes=path.stat().st_size, SHA256=sha(path)))
                if directory == 'model':
                    assert sha(path) == sha(scratch/'prediction/reference'/name)
                elif name == 'report.json':
                    compress(path, fold/'report.json.gz')
                else:
                    shutil.copy2(path, fold/name)
        shutil.copy2(scratch/'model/receipt.json', fold/'model-receipt.json')
        report = json.loads((scratch/'prediction/report.json').read_text())
        expected = 'predicted' if target in annotations else 'noData'
        assert report['queries'][0]['status'] == expected
        record = dict(gemgroup=gem, target=target, status=expected, files=files,
            trainingTargets=[t['id'] for t in selected['selection']['targets']],
            supportedTrainingTargets=report['supportedTrainingTargets'], countsExcludedBeforeNativeFit=True,
            inputFreezeSHA256=sha(out/'input-freeze.json'))
        write(fold/'manifest.json', record)
        results.append(record)
        if i == 0:
            retained = out/'first-model'
            retained.mkdir()
            for path in (scratch/'model').iterdir():
                compress(path, retained/(path.name+'.gz'))
            run('first-repeat', ['singlecell-target-kernel-predict', scratch/'model', '--plan', fold/'query-input.json', '--output', scratch/'repeat'])
            assert sha(scratch/'repeat/report.json') == sha(scratch/'prediction/report.json')
            assert sha(scratch/'repeat/receipt.json') == sha(scratch/'prediction/receipt.json')
        # Only this driver's newly created, verified per-fold scratch is removed.
        # Predictions/receipts are retained; other models reconstruct from the frozen inputs.
        shutil.rmtree(scratch)
        write(out/'folds.json', results)
        print(f'{i+1}/150 {gem} {target}: {expected}', flush=True)
    write(out/'prediction-freeze.json', dict(status='all-native-predictions-frozen',
        frozenUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), folds=150,
        predictedFolds=sum(x['status']=='predicted' for x in results), binarySHA256=BINARY,
        inputFreezeSHA256=sha(out/'input-freeze.json'), foldManifestSHA256=sha(out/'folds.json'),
        commands=len(commands), allNativeReplaysPassed=True, firstRepeatedPredictionExact=True, scoringStarted=False))


if __name__ == '__main__':
    main()
