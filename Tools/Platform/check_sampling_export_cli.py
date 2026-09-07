#!/usr/bin/env python3
"""Exercise sampling export through real CLI processes using retained native state.

The input is the successful MolecularSamplingNativeBridgeTests receipt. The
checker copies its small store before export/corruption tests and never changes
the native fixture. It does not start MD: terminal resume is already at its block
budget, and the downstream recipe contains only snapshot/electronic operations.
This qualifies a connected finite routing fixture, not a sampled ensemble or rate.
"""
import argparse
import copy
import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    if isinstance(value, dict):
        value = bytes(value['bytes']).hex()
    if not isinstance(value, str) or not re.fullmatch(r'[0-9a-f]{64}', value):
        raise RuntimeError('Expected a SHA-256 artifact identity')
    return value


def inventory(root):
    paths = sorted(x for x in root.rglob('*') if x.is_file())
    if len(paths) > 10000 or sum(x.stat().st_size for x in paths) > 256 * 1024 * 1024:
        raise RuntimeError('Native routing fixture exceeds this bounded checker campaign')
    if any(x.is_symlink() for x in root.rglob('*')):
        raise RuntimeError('The retained native fixture must not contain symlinked entries')
    return {str(x.relative_to(root)): sha(x) for x in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=pathlib.Path, required=True)
    parser.add_argument('--fixture', type=pathlib.Path, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    args = parser.parse_args()
    binary, fixture_path, out = args.binary.resolve(strict=True), args.fixture.resolve(strict=True), args.out.resolve()
    fixture = json.loads(fixture_path.read_text())
    if (fixture.get('schema') != 'numivivo.org/test-evidence/sampling-prefix/v1'
            or fixture.get('status') != 'success'
            or fixture.get('nativeStatus') != 'accepted-metal-sampling-and-exact-prefix-resume'
            or not re.fullmatch(r'[0-9a-f]{40}', fixture.get('sourceCommit', ''))
            or fixture.get('exactResumedCursor') is not True
            or fixture.get('sourceTermination') != 'budgetExhausted'
            or fixture.get('declaredObservableCriteriaSatisfied') is not False
            or fixture.get('completedBlocks') != 4 or len(fixture.get('replicas', [])) != 2):
        raise RuntimeError('Expected a successful retained native B4 fixture with exact B1 continuation')
    source_store = pathlib.Path(fixture['store']).resolve(strict=True)
    request_source = pathlib.Path(fixture['requestFile']).resolve(strict=True)
    if out == source_store or source_store in out.parents or out == fixture_path.parent or fixture_path.parent in out.parents:
        raise RuntimeError('Keep checker outputs outside the retained native fixture')
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; existing results are not deleted')
    original_inventory = inventory(source_store)
    out.mkdir(parents=True, exist_ok=True)
    store = out / 'working-store'
    shutil.copytree(source_store, store)
    request_path = out / 'request.json'
    shutil.copyfile(request_source, request_path)
    if sha(request_path) != fingerprint(fixture['request']):
        raise RuntimeError('Retained request bytes disagree with their source artifact')
    checks = []
    report = {'schema': 'numivivo.org/test-evidence/sampling-export-cli/v1',
              'fixtureSHA256': sha(fixture_path), 'fixtureSourceCommit': fixture['sourceCommit'],
              'numericalContract': fixture['numericalContract'], 'binarySHA256': sha(binary),
              'sourceStore': str(source_store), 'workingStore': str(store), 'checks': checks,
              'scope': 'Accepted finite synthetic harmonic prefix to exact snapshot/electronic workflow; no ensemble, rate or performance claim',
              'passed': False}

    def save():
        (out / 'checks.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

    def check(condition, label):
        if not condition:
            report['passed'] = False
        checks.append({'label': label, 'passed': bool(condition)})
        save()
        if not condition:
            raise RuntimeError(label)

    def run(label, *arguments, expected=0):
        command = [str(binary), *map(str, arguments)]
        result = subprocess.run(command, capture_output=True, text=True, timeout=900)
        (out / (label + '.stdout')).write_text(result.stdout)
        (out / (label + '.stderr')).write_text(result.stderr)
        checks.append({'label': label + ' exit', 'command': command, 'expectedExit': expected,
                       'exitCode': result.returncode, 'passed': result.returncode == expected})
        save()
        if result.returncode != expected:
            raise RuntimeError(label + ': unexpected exit; inspect retained logs')
        return result

    def write(name, value):
        path = out / name
        with path.open('x') as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write('\n')
        return path

    def read(name):
        return json.loads((out / name).read_text())

    def object_path(identity):
        digest = fingerprint(identity)
        return store / 'objects/sha256' / digest[:2] / digest[2:4] / digest

    def object_json(identity):
        path = object_path(identity)
        if sha(path) != fingerprint(identity):
            raise RuntimeError('Fixture object hash mismatch')
        return json.loads(path.read_text())

    cursor, source_receipt = fingerprint(fixture['checkpoint']), fingerprint(fixture['sourceReceipt'])
    base = ['molecule-sampling-export', '--store', store, '--checkpoint', cursor, '--replica', '0']
    bound = [*base, '--receipt', source_receipt]
    molecular_help = run('molecular-help', 'molecule-help').stdout
    check(all(x in molecular_help for x in ['molecule-sampling-export', 'molecule-sampling-export-verify',
          '--maximum-export-bytes', '--read-limits', '--require-converged', 'no-clobber']), 'export options discoverable')
    check('molecule-sampling-export-verify' in run('general-help', '--help').stdout, 'global export discovery')

    try:
        for label, command in [
                ('missing-replica', base[:-2]), ('missing-source', ['molecule-sampling-export', '--store', store, '--replica', '0']),
                ('both-sources', [*base, '--reference', fixture['reference']]), ('wrong-replica', [*base[:-1], '2']),
                ('negative-replica', [*base[:-1], '-1']), ('bad-hash', [*base[:4], 'xyz', '--replica', '0']),
                ('duplicate-switch', [*base, '--require-converged', '--require-converged']),
                ('unknown-option', [*base, '--unknown', '1']), ('bad-export-budget', [*base, '--maximum-export-bytes', '0']),
                ('tiny-export-budget', [*base, '--maximum-export-bytes', '1'])]:
            before = inventory(store)
            run(label, *command, expected=65)
            check(inventory(store) == before, label + ' publishes nothing')
        unknown = write('unknown-limits.json', {'schema': 'numivivo.org/molecular-sampling-read-limits/v1', 'maximumReadByte': 99})
        malformed = write('malformed-limits.json', {'schema': 'numivivo.org/molecular-sampling-read-limits/v1', 'maximumReadBytes': None})
        tight = write('tight-limits.json', {'schema': 'numivivo.org/molecular-sampling-read-limits/v1', 'maximumAcceptedSteps': 35})
        cursor_limit = write('cursor-limits.json', {'schema': 'numivivo.org/molecular-sampling-read-limits/v1', 'maximumCursorBytes': 1})
        sufficient = write('sufficient-limits.json', {'schema': 'numivivo.org/molecular-sampling-read-limits/v1', 'maximumAcceptedSteps': 36})
        for name, limits in [('unknown-limits', unknown), ('null-limits', malformed), ('tight-limits', tight)]:
            run(name, *base, '--read-limits', limits, expected=65)
        run('reference-cursor-limit', 'molecule-sampling-export', '--store', store, '--reference', fixture['reference'],
            '--replica', '0', '--read-limits', cursor_limit, expected=65)
        run('require-converged', *bound, '--require-converged', '--output', out / 'not-converged.json', expected=65)
        check(not (out / 'not-converged.json').exists(), 'unmet policy never publishes a success file')
        run('wrong-source-receipt', *base, '--receipt', fixture['request'], expected=65)
        zero = None
        for path in (store / 'descriptors/sha256').rglob('*.json'):
            descriptor = json.loads(path.read_text())
            if descriptor['kind'] == 'molecular-sampling-checkpoint':
                value = object_json(descriptor['fingerprint'])
                if value['completedBlocks'] == 0:
                    zero = fingerprint(descriptor['fingerprint'])
        check(zero is not None, 'actual native block-zero cursor retained')
        run('zero-block', 'molecule-sampling-export', '--store', store, '--checkpoint', zero, '--replica', '0', expected=65)

        run('selected-zero', *bound, '--read-limits', sufficient, '--output', out / 'selected.json',
            '--checkpoint-output', out / 'checkpoint.json', '--snapshot-output', out / 'snapshot.json',
            '--mapping-output', out / 'mapping.json')
        selected = read('selected.json')
        provenance = selected['provenance']
        replica = fixture['replicas'][0]
        checkpoint = object_json(replica['checkpoint'])
        check((out / 'checkpoint.json').read_bytes() == object_path(replica['checkpoint']).read_bytes(), 'original accepted checkpoint bytes retained')
        check(provenance['replicaIndex'] == 0 and provenance['replicaSeed'] == replica['seed']
              and fingerprint(provenance['mdCheckpointFingerprint']) == fingerprint(replica['checkpoint'])
              and provenance['acceptedStep'] == replica['acceptedStep'] and provenance['timePS'] == replica['timePS'],
              'selected replica seed, step and original clock bind native checkpoint')
        check(provenance['sourceTermination'] == 'budgetExhausted' and not provenance['declaredObservableCriteriaSatisfied']
              and provenance['completedBlocks'] == 4 and provenance['validationScope'] == 'restart',
              'successful export retains partial sampling status and recorded payload scope')
        snapshot = read('snapshot.json')
        check(snapshot['stepIndex'] == checkpoint['acceptedStep'] and all(snapshot.get(k) == checkpoint.get(k)
              for k in ['systemFingerprint', 'configurationFingerprint', 'timePS', 'positionsNM', 'velocitiesNMPerPS', 'periodicCell']),
              'full-particle snapshot preserves geometry, velocity, cell and accepted clock')
        check(read('mapping.json')['atomToParticle'] == [0, 1], 'explicit atom-to-particle mapping retained')
        ports = {x['name']: x for x in selected['outputs']}
        check(ports['checkpoint']['kind'] == 'vivo.md-checkpoint'
              and fingerprint(ports['checkpoint']['artifact']) != fingerprint(replica['checkpoint']),
              'workflow checkpoint envelope avoids immutable raw-kind collision')

        before = inventory(store)
        run('verify', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store, '--output', out / 'verified.json')
        run('verify-stronger', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store,
            '--verify-all-payloads', '--output', out / 'verified-all.json')
        stronger = read('verified-all.json')
        check(stronger['recordedValidationScope'] == 'restart' and stronger['validationScope'] == 'allPayloads'
              and stronger['verifiedPayloads'] == provenance['trajectoryChunkCount'], 'stronger fresh verification preserves original export scope')
        check(inventory(store) == before, 'fresh export verification is read-only')
        run('verify-unmet-policy', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store, '--require-converged', expected=65)
        run('reference-reuse', 'molecule-sampling-export', '--store', store, '--reference', fixture['reference'], '--replica', '0',
            '--receipt', source_receipt, '--output', out / 'reference-selected.json')
        check(read('reference-selected.json') == selected, 'one resolved reference selects identical immutable output')
        run('selected-one', 'molecule-sampling-export', '--store', store, '--checkpoint', cursor, '--replica', '1',
            '--verify-all-payloads', '--output', out / 'replica-one.json')
        second = read('replica-one.json')['provenance']
        check(second['replicaIndex'] == 1 and second['sourceTermination'] == 'notRecorded'
              and second['timePS'] == fixture['replicas'][1]['timePS'] and second['validationScope'] == 'allPayloads',
              'second replica and absent termination receipt remain explicit')

        linked_store = out / 'store-link'; linked_store.symlink_to(store, target_is_directory=True)
        alias = out / 'export-hardlink.json'; os.link(out / 'selected.json', alias)
        dangling = out / 'dangling.json'; dangling.symlink_to(store / 'not-created.json')
        for label, extra in [
                ('no-clobber', ['--output', out / 'selected.json']),
                ('store-output', ['--output', store / 'absent' / 'new.json']),
                ('symlink-store-output', ['--output', linked_store / 'absent' / 'new.json']),
                ('dangling-output', ['--output', dangling]),
                ('overlapping-files', ['--output', out / 'overlap', '--checkpoint-output', out / 'overlap' / 'child.json']),
                ('identical-files', ['--output', out / 'same.json', '--checkpoint-output', out / 'same.json']),
                ('existing-optional', ['--output', out / 'must-not-exist.json', '--checkpoint-output', out / 'checkpoint.json']),
                ('limits-input-alias', ['--read-limits', sufficient, '--output', sufficient])]:
            before = inventory(store)
            run(label, *bound, *extra, expected=65)
            check(inventory(store) == before, label + ' rejects before artifact publication')
        run('verify-input-alias', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store, '--output', alias, expected=65)
        check(not (out / 'must-not-exist.json').exists() and not (out / 'overlap').exists(), 'collective output preflight avoids partial file publication')
        check(read('selected.json') == selected, 'no-clobber failures preserve the original export file')

        # Corrupt only the oldest payload in the copied fixture. A restart-scope
        # verification is deliberately narrower than a fresh all-payload check.
        manifest = object_json(replica['trajectory'])
        link = object_json(manifest['tail'])
        while link.get('previous') is not None:
            link = object_json(link['previous'])
        payload = object_path(link['payload']); original = payload.read_bytes()
        try:
            payload.write_bytes(bytes([original[0] ^ 1]) + original[1:])
            run('restart-scope-old-payload', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store)
            run('strong-scope-old-payload', 'molecule-sampling-export-verify', out / 'selected.json', '--store', store, '--verify-all-payloads', expected=65)
        finally:
            payload.write_bytes(original)
        forged = copy.deepcopy(selected); forged['provenance']['acceptedStep'] += 1
        run('forged-export', 'molecule-sampling-export-verify', write('forged-export.json', forged), '--store', store, expected=65)
        run('terminal-resume-tight', 'molecule-sampling-run', request_path, '--store', store, '--resume', cursor,
            '--read-limits', tight, expected=65)
        terminal = json.loads(run('terminal-resume', 'molecule-sampling-run', request_path, '--store', store,
            '--resume', cursor, '--read-limits', sufficient, expected=75).stdout)
        check(fingerprint(terminal['checkpoint']) == cursor and terminal['status'] == 'budgetExhausted',
              'real CLI resume forwards reader limits and preserves terminal accepted prefix without MD')

        # Existing template supplies explicit H2 STO-3G and one alpha/one beta
        # electron. Replace only the MD producer with verified exported inputs.
        recipe = json.loads(run('electronic-template', 'workflow-template', 'md-electronic-analysis').stdout)
        source_request = json.loads(request_path.read_text())
        check(len(source_request['structure']['atoms']) == 2
              and all(x['element']['atomicNumber'] == 1 for x in source_request['structure']['atoms']),
              'downstream H2 electronic settings match this finite routing fixture')
        recipe['identifier'] = 'sampling-selected-checkpoint-electronic-route'
        recipe['artifacts'] = [x for x in recipe['artifacts'] if x['identifier'] == 'basis']
        for identifier, output in [('structure', 'source-structure'), ('system', 'system'), ('accepted-checkpoint', 'checkpoint')]:
            port = ports[output]
            recipe['artifacts'].append({'identifier': identifier, 'source': {'stored': {'kind': port['kind'], 'fingerprint': port['artifact']}}})
        recipe['nodes'] = [x for x in recipe['nodes'] if x['identifier'] not in ['segment-1', 'segment-2']]
        snapshot_node = next(x for x in recipe['nodes'] if x['identifier'] == 'snapshot')
        snapshot_node['inputs']['checkpoint'] = {'artifact': {'identifier': 'accepted-checkpoint'}}
        recipe['outputs'] = [x for x in recipe['outputs'] if x['node'] not in ['segment-1', 'segment-2']]
        check(len(recipe['nodes']) == 7 and not any(x['operation'] in ['vivo.platform.md-start', 'vivo.platform.md-continue'] for x in recipe['nodes']),
              'downstream recipe contains snapshot and electronics only, with no MD engine')
        recipe_path = write('selected-electronic-recipe.json', recipe)
        run('electronic-plan', 'workflow-plan', recipe_path)
        run('electronic-route', 'workflow-run', recipe_path, '--store', store, '--output', out / 'electronic.json')
        electronic = read('electronic.json')
        check(electronic['allTasksSucceeded'] and len(electronic['nodes']) == 7, 'exact accepted checkpoint reaches all existing electronic stages')
        mapping = next(x['artifact'] for x in electronic['exports'] if x['name'] == 'snapshot-mapping')
        mapped = json.loads(run('electronic-mapping', 'workflow-export', fingerprint(mapping['artifact']), '--kind', mapping['kind'], '--store', store).stdout)
        check(mapped == read('mapping.json'), 'downstream mapper consumes the identical original accepted checkpoint')
        run('electronic-reuse', 'workflow-run', recipe_path, '--store', store, '--output', out / 'electronic-reused.json')
        check(all(x['outcome']['succeeded']['reused'] for x in read('electronic-reused.json')['nodes']), 'unchanged connected electronic route reuses verified tasks')
        report_id = fingerprint(read('electronic.json.receipt.json')['reportArtifact'])
        run('electronic-verify', 'workflow-verify', report_id, '--store', store, '--output', out / 'electronic-verified.json')
        check(read('electronic-verified.json')['allTasksSucceeded'], 'workflow and export share the same executable implementation identity')
        report['passed'] = True
    finally:
        check(inventory(source_store) == original_inventory, 'original native fixture remains byte-identical')
        save()
    print('PASS', len(checks), 'real sampling export/verify and connected workflow checks; evidence in', out)


if __name__ == '__main__':
    main()
