#!/usr/bin/env python3
"""Exercise real general-workflow commands; no numerical engine is mocked."""
import argparse
import copy
import json
import pathlib
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('binary', type=pathlib.Path)
    parser.add_argument('out', type=pathlib.Path)
    args = parser.parse_args()
    binary, out = args.binary.resolve(), args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; existing data are not deleted')
    out.mkdir(parents=True, exist_ok=True)
    store = out / 'store'
    checks = []

    def check(condition, label):
        checks.append({'label': label, 'passed': bool(condition)})
        (out / 'checks.json').write_text(json.dumps({'checks': checks, 'passed': all(x['passed'] for x in checks)}, indent=2))
        if not condition:
            raise RuntimeError(label)

    def run(label, *arguments, expected=0):
        result = subprocess.run([str(binary), *map(str, arguments)], capture_output=True, text=True, timeout=900)
        (out / (label + '.log')).write_text(result.stdout + '\n' + result.stderr)
        check(result.returncode == expected, label + ' exit status')
        return result

    def write(name, value):
        path = out / name
        path.write_text(json.dumps(value, sort_keys=True, allow_nan=False) + '\n')
        return path

    def read(name):
        return json.loads((out / name).read_text())

    def digest(value):
        return bytes(value['bytes']).hex()

    def outcome(report, node):
        return next(x['outcome'] for x in report['nodes'] if x['identifier'] == node)

    catalog = json.loads(run('catalog', 'workflow-catalog').stdout)
    registered = {x['identifier'] for x in catalog}
    check(len(registered) == len(catalog), 'registered operation identifiers unique')
    check({'vivo.platform.md-start', 'vivo.platform.md-continue', 'vivo.platform.target-reference',
           'vivo.platform.structure-electronic-system', 'vivo.native.advanced-many-body',
           'vivo.platform.qmmm-transmission-analyze', 'vivo.platform.qmmm-apply-transmission',
           'vivo.platform.qmmm-chemical-qualification', 'vivo.platform.qmmm-chemical-state-network',
           'vivo.platform.qmmm-chemical-exchange-network'}.issubset(registered),
          'catalog exposes cross-domain, transmission and chemical-qualification operation families')
    qmmm_help = run('qmmm-help', 'qmmm-free-energy-help').stdout
    check('qmmm-transmission-analyze' in qmmm_help and 'qmmm-transmission-apply' in qmmm_help and
          'qmmm-chemical-qualify' in qmmm_help and 'qmmm-chemical-state-network' in qmmm_help and
          'qmmm-chemical-exchange-network' in qmmm_help,
          'real CLI routes transmission, qualification, rapid-equilibrium and transient-exchange commands')
    recipe_path = out / 'recipe.json'
    run('template', 'workflow-template', 'molecular-analysis', '--output', recipe_path)
    recipe = read('recipe.json')
    plan = json.loads(run('plan', 'workflow-plan', recipe_path).stdout)
    check(len(plan['nodes']) == 6 and not store.exists(), 'planning six typed stages performs no store writes')
    run('first', 'workflow-run', recipe_path, '--store', store, '--output', out / 'first.json')
    first = read('first.json')
    check(first['allTasksSucceeded'] and len(first['nodes']) == 6, 'complete molecular pipeline executes')
    check(first['maximumAdmittedTasks'] == 2, 'independent many-body tasks admitted together')
    check(all('succeeded' in n['outcome'] for n in first['nodes']), 'all stages have successful receipts')
    check(not any(n['outcome']['succeeded']['reused'] for n in first['nodes']), 'fresh pipeline actually executes')
    run('second', 'workflow-run', recipe_path, '--store', store, '--output', out / 'second.json')
    second = read('second.json')
    check(all(n['outcome']['succeeded']['reused'] for n in second['nodes']), 'every unchanged stage reuses verified outputs')
    check(first['exports'] == second['exports'], 'exports retain identical artifact identities on resume')
    report_id = digest(read('first.json.receipt.json')['reportArtifact'])
    run('verify', 'workflow-verify', report_id, '--store', store, '--output', out / 'verified.json')
    check(read('verified.json')['allTasksSucceeded'], 'stored report verifies against task reconstruction')
    fci = next(x['artifact'] for x in first['exports'] if x['name'] == 'fci')
    fci_id = digest(fci['artifact'])
    run('export', 'workflow-export', fci_id, '--kind', fci['kind'], '--store', store, '--output', out / 'fci.json')
    check('configurationInteraction' in read('fci.json'), 'export contains the actual selected numerical method')
    run('wrong-export-kind', 'workflow-export', fci_id, '--kind', 'wrong.kind', '--store', store, expected=65)
    run('non-ascii-digest', 'workflow-export', 'é' * 32, '--kind', fci['kind'], '--store', store, expected=65)
    original = recipe_path.read_bytes()
    run('input-alias', 'workflow-run', recipe_path, '--output', recipe_path, '--force', '--store', store, expected=65)
    check(recipe_path.read_bytes() == original, 'input bytes survive rejected overwrite')
    link = out / 'recipe-hardlink.json'
    link.hardlink_to(recipe_path)
    run('hardlink-alias', 'workflow-run', recipe_path, '--output', link, '--force', '--store', store, expected=65)
    run('store-alias', 'workflow-run', recipe_path, '--output', store / 'uncreated' / 'output.json', '--store', store, expected=65)
    store_link = out / 'linked-store'
    store_link.symlink_to(store, target_is_directory=True)
    run('store-symlink-alias', 'workflow-run', recipe_path, '--output', store_link / 'new-leaf.json', '--store', store, expected=65)
    run('no-clobber', 'workflow-run', recipe_path, '--output', out / 'first.json', '--store', store, expected=65)
    malformed = copy.deepcopy(recipe)
    malformed['nodes'][0]['operation'] = 'not.registered'
    path = write('unknown.json', malformed)
    run('unknown-preflight', 'workflow-run', path, '--store', out / 'never-created', expected=65)
    check(not (out / 'never-created').exists(), 'invalid operation rejected before artifact-store creation')
    malformed = copy.deepcopy(recipe)
    malformed['nodes'][0]['version'] = '999'
    run('version-preflight', 'workflow-plan', write('version.json', malformed), expected=65)
    failed = copy.deepcopy(recipe)
    reference = next(n for n in failed['nodes'] if n['identifier'] == 'reference')
    independent = copy.deepcopy(reference)
    independent['identifier'] = 'independent-reference'
    failed['nodes'].append(independent)
    reference['configuration']['maximumIterations'] = 1
    run('partial', 'workflow-run', write('partial-recipe.json', failed), '--store', store, '--output', out / 'partial.json', expected=1)
    partial = read('partial.json')
    check(not partial['allTasksSucceeded'], 'partial report never becomes a completed workflow')
    check('failed' in outcome(partial, 'reference') and 'blocked' in outcome(partial, 'hamiltonian'),
          'numerical failure and blocked dependent are distinct')
    check('succeeded' in outcome(partial, 'independent-reference'), 'unrelated valid calculation retained')
    check(all(x.get('artifact') is None for x in partial['exports']), 'failed branch exports remain absent')
    forged = copy.deepcopy(first)
    forged['allTasksSucceeded'] = False
    imported = json.loads(run('import-forged', 'workflow-import', write('forged.json', forged), '--kind', 'vivo.workflow-run', '--store', store).stdout)
    run('reject-forged', 'workflow-verify', digest(imported['fingerprint']), '--store', store, expected=65)
    object_path = store / 'objects' / 'sha256' / fci_id[:2] / fci_id[2:4] / fci_id
    original = object_path.read_bytes()
    try:
        object_path.write_bytes(original + b'\n')
        run('corrupt-cached-output', 'workflow-run', recipe_path, '--store', store, '--output', out / 'corrupt.json', expected=1)
        check('failed' in outcome(read('corrupt.json'), 'fci'), 'corrupt cache is retained as a failed node, not silently recomputed')
    finally:
        object_path.write_bytes(original)
    run('restored-cache', 'workflow-run', recipe_path, '--store', store, '--output', out / 'restored.json')
    sampled = out / 'sampled-recipe.json'
    run('sampled-template', 'workflow-template', 'md-electronic-analysis', '--output', sampled)
    run('sampled-run', 'workflow-run', sampled, '--store', store, '--output', out / 'sampled.json')
    report = read('sampled.json')
    check(report['allTasksSucceeded'] and len(report['nodes']) == 9, 'Metal-to-electronic workflow executes all nine typed stages')
    item = next(x['artifact'] for x in report['exports'] if x['name'] == 'snapshot-mapping')
    run('sampled-mapping', 'workflow-export', digest(item['artifact']), '--kind', item['kind'], '--store', store, '--output', out / 'mapping.json')
    mapping = read('mapping.json')
    check(mapping['atomToParticle'] == [0, 1] and mapping['step'] == 20, 'production snapshot preserves particle mapping and accepted step')
    check(mapping['sourceStructure'] != mapping['snapshotStructure'], 'new geometry has a distinct source-bound structure identity')
    run('sampled-resume', 'workflow-run', sampled, '--store', store, '--output', out / 'sampled-cached.json')
    check(all(x['outcome']['succeeded']['reused'] for x in read('sampled-cached.json')['nodes']), 'all MD and electronic stages resume from verified artifacts')
    print('PASS', len(checks), 'real CLI, planner, cache, export and failure-preservation checks')


if __name__ == '__main__':
    main()
