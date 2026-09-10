#!/usr/bin/env python3
"""Qualify immutable graph-reference reuse on complete existing Baron bundles."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_text())


def replace(path, raw):
    # Scratch trees share immutable inputs by hard link. Never write through one.
    temporary = path.with_name(path.name+'.replacement')
    temporary.write_bytes(raw)
    os.replace(temporary, path)


def update(path, edit):
    value = read(path)
    edit(value)
    replace(path, (json.dumps(value, sort_keys=True)+'\n').encode())


def bind(path, key, source):
    update(path, lambda value: value.update({key: {'bytes': list(bytes.fromhex(sha(source)))}}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    original = root/'clustering/baron'
    current = root/'clustering-storage-adaptive/baron'
    reference = root/'clustering/baron-reference'
    protocol = root/'clustering/baron-protocol.json'
    graph_check = root/'graph/baron-independent-check-final/checks.json'
    checker = Path(__file__).with_name('reference_clustering.py')
    assert read(protocol)['cells'] == 8569
    assert sha(original/'input/receipt.json') != sha(current/'input/receipt.json')
    args.out.mkdir(parents=True, exist_ok=False)
    inputs = [p for folder in [original, current, reference] for p in folder.rglob('*') if p.is_file()]
    inputs += [protocol, graph_check, checker, Path(__file__)]
    before = {str(p): sha(p) for p in inputs}
    records = []

    def run(name, bundle, reuse, expected, ref=reference, reason=None):
        command = [sys.executable, str(checker), 'check', '--bundle', str(bundle),
                   '--protocol', str(protocol), '--graph-check', str(graph_check),
                   '--reference', str(ref), '--out', str(args.out/name)]
        if reuse:
            command += ['--reference-graph', str(original/'input')]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (args.out/(name+'.log')).write_text(result.stdout)
        assert (result.returncode == 0) == expected, (name, result.stdout)
        if reason:
            assert reason in result.stdout, (name, result.stdout)
        record = dict(name=name, returnCode=result.returncode, expectedSuccess=expected,
                      command=command, logSHA256=sha(args.out/(name+'.log')))
        if expected:
            checked = read(args.out/name/'checks.json')
            assert checked['status'] == 'passed' and checked['cells'] == 8569
            assert checked['inputBindingsStableThroughoutCheck']
            record['checksSHA256'] = sha(args.out/name/'checks.json')
        else:
            assert not (args.out/name/'checks.json').exists()
        records.append(record)

    run('original-direct', original, False, True)
    run('current-explicit-reuse', current, True, True)
    old = read(args.out/'original-direct/checks.json')
    new = read(args.out/'current-explicit-reuse/checks.json')
    for key in ['resultSHA256', 'clusters', 'independentModularity', 'references',
                'referencePairwiseARI', 'disconnectedCommunities', 'edgeVisits']:
        assert old[key] == new[key], key
    assert len(new['referenceReuse']['equivalentPayloads']) == 13
    run('current-without-equivalence-rejected', current, False, False)

    # Every altered artifact is internally rehashed. Equality must still reject
    # it rather than accepting a new receipt as the original numerical evidence.
    mutations = ['graph-payload', 'pca-payload', 'source-identity', 'pca-implementation',
                 'unfitted-graph', 'reference-partitions']
    for name in mutations:
        with tempfile.TemporaryDirectory(prefix='numivivo-reference-reuse-') as temporary:
            scratch = Path(temporary)
            bundle = scratch/'bundle'
            shutil.copytree(current, bundle, copy_function=os.link)
            graph = bundle/'input'
            pca = graph/'input'
            ref = reference
            reason = None
            if name == 'graph-payload':
                replace(graph/'edges.bin', (graph/'edges.bin').read_bytes()[:-1]+b'\x00')
                assert sha(graph/'edges.bin') != sha(current/'input/edges.bin')
                bind(graph/'graph.json', 'edges', graph/'edges.bin')
                bind(graph/'receipt.json', 'graph', graph/'graph.json')
                reason = 'Reference reuse'
            elif name == 'pca-payload':
                replace(pca/'model.json', (pca/'model.json').read_bytes()+b'\n')
                bind(pca/'receipt.json', 'model', pca/'model.json')
                reason = 'Reference reuse'
            elif name == 'source-identity':
                update(pca/'receipt.json', lambda value: value['source']['bytes'].__setitem__(0, value['source']['bytes'][0]^1))
                reason = 'Reference reuse source identity'
            elif name == 'pca-implementation':
                update(pca/'receipt.json', lambda value: value['implementation'].update({'qualificationTamper': True}))
                reason = 'Current PCA/graph implementation mismatch'
            elif name == 'unfitted-graph':
                update(graph/'plan.json', lambda value: value.update({'inputKind': 'query'}))
                bind(graph/'receipt.json', 'plan', graph/'plan.json')
            else:
                ref = scratch/'reference'
                shutil.copytree(reference, ref, copy_function=os.link)
                replace(ref/'partitions.npz', (ref/'partitions.npz').read_bytes()+b'changed')
            bind(graph/'receipt.json', 'input', pca/'receipt.json')
            bind(bundle/'receipt.json', 'input', graph/'receipt.json')
            run(name+'-rejected', bundle, True, False, ref=ref, reason=reason)

    assert all(sha(Path(path)) == digest for path, digest in before.items()), 'Original inputs changed'
    result = dict(status='passed', cells=8569, commands=records, sourceIdentities=before,
                  originalArtifactsUnchanged=True, allThirteenPayloadsCompared=True,
                  scope='Complete Baron reference applicability, numerical partition checks and rehashed corruption rejection. No completed HIRISA clustering or biological qualification is claimed.')
    (args.out/'qualification.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(status='passed', commands=len(records), expectedRejections=7,
                         cells=8569, qualificationSHA256=sha(args.out/'qualification.json'))))


if __name__ == '__main__':
    main()
