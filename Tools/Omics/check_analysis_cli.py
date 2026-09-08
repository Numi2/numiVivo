#!/usr/bin/env python3
"""Exercise the real native cohort CLI; Python supplies no statistical fit."""
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--out', required=True, type=Path)
    args = p.parse_args()
    binary, out = args.binary.resolve(strict=True), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    store, example = out / 'store', out / 'example'
    commands, checks = [], []
    passed = False

    def save(path, value):
        path.write_text(json.dumps(value, indent=2) + '\n')
        return path

    def load(path):
        return json.loads(path.read_text())

    def run(*values, ok=True):
        command = [str(binary), *map(str, values)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=180)
        index = len(commands)
        (out / f'{index:03d}.stdout').write_text(result.stdout)
        (out / f'{index:03d}.stderr').write_text(result.stderr)
        commands.append({'arguments': command, 'exitCode': result.returncode, 'expectedSuccess': ok})
        if (result.returncode == 0) != ok:
            raise RuntimeError(f'Unexpected exit: {command}\n{result.stderr[-8000:]}')
        return result.stdout

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def digest(value):
        return bytes(value['bytes']).hex()

    try:
        run('singlecell-example', '--output', example)
        counts, analysis, report_file = out / 'counts.json', out / 'analysis.json', out / 'report.json'
        run('singlecell-run', example / 'manifest.json', '--store', store, '--output', counts)
        run('singlecell-analyze', counts, '--plan', example / 'analysis.json', '--store', store, '--output', analysis)
        run('singlecell-analysis-verify', analysis, '--store', store)
        run('singlecell-analysis-export', analysis, '--store', store, '--output', report_file)
        report = load(report_file)
        processed, contrast = report['processed'], report['contrasts'][0]
        check('QC retains all decisions and selected source indices', len(processed['decisions']) == 36 and len(processed['sourceCellIndices']) == 24)
        check('paired design uses six donors, not 24 independent cells', contrast['design']['controlReplicates'] == 6 and contrast['design']['residualDegreesOfFreedom'] == 5)
        check('all fixture genes tested with explicit empirical prior', contrast['testedFeatures'] == 32 and contrast.get('variancePrior') is not None)
        check('declared synthetic signal recovered', all(x['log2FoldChange'] > 1.5 for x in contrast['features'][:3]))
        check('source evidence stays synthetic', processed['dataset']['evidence'] == 'synthetic')
        tables, mex = out / 'tables', out / 'mex'
        run('singlecell-analysis-tables', analysis, '--store', store, '--output', tables)
        check('tables include design and missing-value convention', (tables / 'contrast-0000-design.tsv').exists() and 'empty TSV fields' in load(tables / 'index.json')['missingValues'])
        run('singlecell-analysis-mex', analysis, '--store', store, '--output', mex)
        reimport, exported = out / 'reimport.json', out / 'reimport-report.json'
        run('singlecell-run', mex / 'manifest.json', '--store', store, '--output', reimport)
        run('singlecell-export', reimport, '--store', store, '--output', exported)
        check('MEX round trip preserves selected raw dataset', load(exported)['dataset'] == processed['dataset'])
        run('singlecell-analysis-tables', analysis, '--store', store, '--output', tables, ok=False)
        manifest = load(example / 'manifest.json')
        paths = {lib[key] for lib in manifest['libraries'] for key in ['matrix', 'features', 'barcodes']}
        for path in paths:
            (example / (path + '.gz')).write_bytes(gzip.compress((example / path).read_bytes(), mtime=0))
        for lib in manifest['libraries']:
            for key in ['matrix', 'features', 'barcodes']:
                lib[key] += '.gz'
        compressed = save(example / 'compressed.json', manifest)
        gzip_receipt, gzip_report, raw_report = out / 'gzip.json', out / 'gzip-report.json', out / 'raw-report.json'
        run('singlecell-run', compressed, '--store', store, '--output', gzip_receipt)
        run('singlecell-export', gzip_receipt, '--store', store, '--output', gzip_report)
        run('singlecell-export', counts, '--store', store, '--output', raw_report)
        check('compressed sources reconstruct identical raw counts', load(gzip_report)['dataset'] == load(raw_report)['dataset'])
        check('compressed source identity is not erased', load(gzip_receipt)['input'] != load(counts)['input'])
        bad = load(example / 'analysis.json')
        bad['contrasts'][0]['design'] = 'independentReplicates'
        bad_output = out / 'must-not-publish.json'
        run('singlecell-analyze', counts, '--plan', save(out / 'bad-plan.json', bad), '--store', store, '--output', bad_output, ok=False)
        check('invalid replication publishes no success receipt', not bad_output.exists())
        # Exercise the existing workflow scheduler using stored source/plan bytes.
        descriptor = json.loads(run('artifact-put', example / 'analysis.json', '--kind', 'vivo.singlecell-analysis-plan-v1', '--media-type', 'application/json', '--store', store))
        recipe = {
            'schema': 'numivivo.org/workflow-recipe/v1', 'identifier': 'native-singlecell-analysis',
            'artifacts': [
                {'identifier': 'counts', 'source': {'stored': {'kind': 'vivo.singlecell-input-bundle-v1', 'fingerprint': load(counts)['input']}}},
                {'identifier': 'plan', 'source': {'stored': {'kind': 'vivo.singlecell-analysis-plan-v1', 'fingerprint': descriptor['fingerprint']}}}],
            'nodes': [{'identifier': 'analyze', 'operation': 'vivo.platform.singlecell-analyze', 'version': '1', 'configuration': {},
                'inputs': {key: {'artifact': {'identifier': value}} for key, value in [('input', 'counts'), ('plan', 'plan')]},
                'resources': {'budget': {'maximumBasisFunctions': 64, 'maximumBytes': 1073741824, 'maximumDeterminants': 512, 'maximumOperatorApplications': 100000000},
                              'maximumInputBytes': 134217728, 'maximumOutputBytes': 536870912, 'numericalBackend': 'cpu-fp64'}}],
            'outputs': [{'name': 'analysis', 'node': 'analyze', 'port': 'report'}],
            'policy': {'maximumConcurrentMetalTasks': 0, 'maximumConcurrentTasks': 1, 'maximumInlineBytes': 134217728, 'maximumNodes': 256, 'maximumReservedBytes': 2147483648}}
        recipe_path, workflow = save(out / 'workflow.json', recipe), out / 'workflow-run.json'
        run('workflow-run', recipe_path, '--store', store, '--output', workflow)
        check('analysis is connected to common DAG', load(workflow)['allTasksSucceeded'])
        workflow_receipt = load(Path(str(workflow) + '.receipt.json'))
        dag_report = out / 'workflow-report.json'
        run('workflow-export', digest(workflow_receipt['reportArtifact']), '--name', 'analysis', '--store', store, '--output', dag_report)
        check('DAG and dedicated CLI share one analysis authority', load(dag_report) == report)
        receipt = load(analysis)
        result_descriptor = json.loads(run('artifact-show', digest(receipt['result']), '--store', store))
        false_result = load(store / result_descriptor['objectPath'])
        false_result['report']['processed']['decisions'][0]['accepted'] = False
        false_descriptor = json.loads(run('artifact-put', save(out / 'false-result.json', false_result), '--kind', 'vivo.singlecell-analysis-result-v1', '--media-type', 'application/json', '--store', store))
        wrong_receipt = copy.deepcopy(receipt); wrong_receipt['result'] = false_descriptor['fingerprint']
        run('singlecell-analysis-verify', save(out / 'false-receipt.json', wrong_receipt), '--store', store, ok=False)
        checks.append('correctly hashed false analysis is rejected by reconstruction')
        passed = True
    finally:
        save(out / 'observations.json', {'schema': 'numivivo.org/singlecell-analysis-cli-checks/v1', 'passed': passed,
            'checks': checks, 'commands': commands, 'binarySHA256': hashlib.sha256(binary.read_bytes()).hexdigest(),
            'scope': 'Synthetic native execution and interchange checks; not external biological calibration.'})
    print(f'Native cohort CLI checks passed: {len(checks)} assertions, {len(commands)} commands')


if __name__ == '__main__':
    main()
