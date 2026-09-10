#!/usr/bin/env python3
"""Recompute saved full-graph igraph partition objectives with NumPy blocks."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
from reference_clustering import read, sha, source, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['graph', 'protocol', 'graph-check', 'reference', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    protocol = read(args.protocol)
    record = read(args.reference/'reference.json')
    bindings = {name: sha(path) for name, path in {
        'graphReceiptSHA256': args.graph/'receipt.json', 'graphSHA256': args.graph/'graph.json',
        'protocolSHA256': args.protocol, 'graphCheckSHA256': args.graph_check,
        'partitionsSHA256': args.reference/'partitions.npz'}.items()}
    assert record['status'] == 'passed' and all(record[k] == v for k, v in bindings.items())
    reference_sha = sha(args.reference/'reference.json')
    report, edges = source(args.graph, protocol, read(args.graph_check))
    n = report['cells']
    assert record['cells'] == n and record['undirectedEdges']*2 == len(edges)
    assert record['versions']['igraph'] == protocol['referenceVersion']
    labels = {}
    with np.load(args.reference/'partitions.npz', allow_pickle=False) as archive:
        assert set(archive.files) == {'seed'+str(seed) for seed in protocol['referenceSeeds']}
        for seed in protocol['referenceSeeds']:
            values = archive['seed'+str(seed)]
            assert values.shape == (n,) and np.issubdtype(values.dtype, np.integer)
            assert values.min() == 0 and values.max() < n
            labels[seed] = values
    degrees = np.zeros(n)
    internal = dict.fromkeys(labels, 0.0)
    for begin in range(0, len(edges), 1048576):
        block = edges[begin:begin+1048576]
        degrees += np.bincount(block['row'], weights=block['value'], minlength=n)
        for seed, values in labels.items():
            within = values[block['row']] == values[block['column']]
            internal[seed] += float(block['value'][within].sum())
    total = float(degrees.sum())
    assert total > 0
    checked = []
    assert [v['seed'] for v in record['references']] == protocol['referenceSeeds']
    for old in record['references']:
        values = labels[old['seed']]
        k = old['clusters']
        assert np.array_equal(np.unique(values), np.arange(k))
        volumes = np.bincount(values, weights=degrees, minlength=k)
        objective = internal[old['seed']]/total-protocol['resolution']*float(np.sum((volumes/total)**2))
        np.testing.assert_allclose(objective, old['modularity'], rtol=1e-10, atol=1e-10)
        checked.append(dict(seed=old['seed'], clusters=k, savedModularity=old['modularity'],
                            independentModularity=objective, absoluteError=abs(objective-old['modularity'])))
    # Recheck the graph and reference identities after traversing every edge.
    _, after_edges = source(args.graph, protocol, read(args.graph_check))
    del after_edges
    assert sha(args.reference/'reference.json') == reference_sha
    assert sha(args.reference/'partitions.npz') == bindings['partitionsSHA256']
    assert sha(args.protocol) == bindings['protocolSHA256'] and sha(args.graph_check) == bindings['graphCheckSHA256']
    result = dict(status='passed', cells=n, directedEdges=len(edges), references=checked,
                  bindings=bindings, referenceSHA256=reference_sha, checkerSHA256=sha(Path(__file__)),
                  seconds=time.monotonic()-start,
                  scope='All saved reference labels and every graph edge enter independent modularity reconstruction. This checks numerical reference integrity; it neither chooses new partitions nor qualifies native clustering or biology.')
    write(args.out/'checks.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
