#!/usr/bin/env python3
"""Stream original CSR counts into source-bound measured expression programs."""
import argparse
import json
from pathlib import Path
import platform
import resource
import time
import h5py
import numpy as np
from scipy import sparse
from check_integration_response import sha


def strings(group):
    assert set(group) == {'values', 'mask'} and not group['mask'][:].any()
    return group['values'].asstr()[:].tolist()


def resolve(ids, names, definitions):
    assert len(ids) == len(names) == len(set(ids))
    positions = {}
    for i, name in enumerate(names):
        positions.setdefault(name, []).append(i)
    weights = np.zeros((len(ids), len(definitions)))
    resolved = []
    for j, definition in enumerate(definitions):
        assert definition['featureNamespace'] == 'HUMAN_GENE_SYMBOL'
        members = definition['members']
        assert len({m['featureID'] for m in members}) == len(members)
        assert all(np.isfinite(m['weight']) and m['weight'] != 0 for m in members)
        requested = {m['featureID']: float(m['weight']) for m in members}
        assert all(len(positions.get(name, [])) <= 1 for name in requested), 'Ambiguous original symbol'
        matched = [i for i, name in enumerate(names) if name in requested]
        missing = sorted(set(requested)-set(names))
        absolute = sum(abs(requested[names[i]]) for i in matched)
        coverage = absolute/sum(abs(v) for v in requested.values())
        assert absolute > 0 and coverage >= definition['minimumWeightCoverage']
        weights[matched, j] = [requested[names[i]]/absolute for i in matched]
        resolved.append(dict(definition=definition, missingSymbols=missing, weightCoverage=coverage,
                             featureIndices=matched, featureIDs=[ids[i] for i in matched],
                             originalSymbols=[names[i] for i in matched], effectiveWeights=weights[matched, j].tolist()))
    return weights, resolved


def score_block(counts, indices, offsets, weights, target):
    n = len(offsets)-1
    assert offsets[0] == 0 and offsets[-1] == len(counts) == len(indices)
    assert (np.diff(offsets) >= 0).all() and (indices >= 0).all() and (indices < len(weights)).all()
    assert np.issubdtype(counts.dtype, np.unsignedinteger)
    # The frozen source uses uint16 counts; this bound prevents uint64 sum overflow.
    assert len(counts)*int(counts.max(initial=0)) < 2**64
    total = np.zeros(n, dtype=np.uint64)
    nonempty = np.flatnonzero(np.diff(offsets))
    if len(nonempty):
        boundaries = np.zeros(len(indices), dtype=bool)
        boundaries[offsets[nonempty]] = True
        assert ((np.diff(indices) > 0) | boundaries[1:]).all(), 'Source CSR is not canonical'
        total[nonempty] = np.add.reduceat(counts.astype(np.uint64), offsets[nonempty])
    assert total.max(initial=0) <= 2**53
    present = np.flatnonzero(np.any(weights != 0, axis=1)[indices] & (counts > 0))
    rows = np.searchsorted(offsets[1:], present, side='right')
    genes = indices[present]
    log_values = np.log1p(counts[present].astype(float)/total[rows]*target)
    scores = np.zeros((n, weights.shape[1]))
    detected = np.zeros_like(scores, dtype=np.uint16)
    for j in range(weights.shape[1]):
        values = log_values*weights[genes, j]
        scores[:, j] = np.bincount(rows, weights=values, minlength=n)
        use = weights[genes, j] != 0
        detection = np.bincount(rows[use], minlength=n)
        assert detection.max(initial=0) <= np.iinfo(np.uint16).max
        detected[:, j] = detection
    # A separate sparse matrix product checks all program score values in the block.
    matrix = sparse.csr_matrix((log_values, (rows, genes)), shape=(n, len(weights)))
    oracle = np.asarray(matrix @ weights)
    np.testing.assert_allclose(scores, oracle, rtol=1e-12, atol=1e-12)
    error = float(np.max(np.abs(scores-oracle), initial=0))
    scores[total == 0] = np.nan
    return total, scores, detected, error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'protocol', 'definitions', 'qc-reference', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    began = time.monotonic()
    protocol = json.loads(args.protocol.read_text())
    bindings = {name: sha(getattr(args, name)) for name in ['source', 'protocol', 'definitions']}
    assert bindings['source'] == protocol['sourceSHA256']
    assert bindings['definitions'] == protocol['definitionsSHA256']
    assert sha(args.qc_reference/'receipt.json') == protocol['qualityReferenceSHA256']
    quality = json.loads((args.qc_reference/'receipt.json').read_text())
    assert quality['sourceSHA256'] == bindings['source'] and quality['cells'] == protocol['cells']
    expected_total = np.zeros(protocol['cells'], dtype=np.uint64)
    end = 0
    for library in quality['libraries']:
        assert library['start'] == end
        path = args.qc_reference/(library['accession']+'.npz')
        assert sha(path) == library['referenceSHA256']
        with np.load(path, allow_pickle=False) as values:
            expected_total[library['start']:library['end']] = values['totalCounts']
        end = library['end']
    assert end == protocol['cells']
    definitions = json.loads(args.definitions.read_text())['definitions']
    assert [d['id'] for d in definitions] == protocol['programs']
    args.out.mkdir(parents=True, exist_ok=False)
    with h5py.File(args.source, 'r') as h:
        ids = strings(h['var/id']); names = strings(h['var/name'])
        assert ids == strings(h['var/_index'])
        weights, resolved = resolve(ids, names, definitions)
        matrix = h['X']
        assert matrix.attrs['encoding-type'] == 'csr_matrix'
        assert tuple(matrix.attrs['shape']) == (protocol['cells'], len(ids))
        assert str(matrix['data'].dtype) == protocol['sourceCountDtype']
        ptr = matrix['indptr'][:]
        assert ptr.shape == (protocol['cells']+1,) and ptr[0] == 0 and (np.diff(ptr) >= 0).all()
        assert ptr[-1] == len(matrix['data']) == len(matrix['indices'])
        scores = np.lib.format.open_memmap(args.out/'scores.npy', mode='w+', dtype='<f8', shape=(protocol['cells'], len(definitions)))
        detected = np.lib.format.open_memmap(args.out/'detected.npy', mode='w+', dtype='<u2', shape=scores.shape)
        total = np.lib.format.open_memmap(args.out/'total-counts.npy', mode='w+', dtype='<u8', shape=(protocol['cells'],))
        start = 0; error = 0.0; blocks = 0; max_rows = 0; max_entries = 0
        while start < protocol['cells']:
            stop = min(start+protocol['maximumSourceRows'], protocol['cells'],
                       int(np.searchsorted(ptr, ptr[start]+protocol['maximumSourceEntries'], side='right')-1))
            assert stop > start, 'One cell exceeds the declared sparse block budget'
            first, last = int(ptr[start]), int(ptr[stop])
            values = matrix['data'][first:last]; indices = matrix['indices'][first:last]
            counts, observed, detection, discrepancy = score_block(values, indices, ptr[start:stop+1]-first, weights, protocol['normalizationTarget'])
            np.testing.assert_array_equal(counts, expected_total[start:stop])
            total[start:stop] = counts; scores[start:stop] = observed; detected[start:stop] = detection
            error = max(error, discrepancy); blocks += 1
            max_rows = max(max_rows, stop-start); max_entries = max(max_entries, last-first)
            start = stop
        scores.flush(); detected.flush(); total.flush()
        empty = int(np.count_nonzero(total == 0))
        umi_total = int(total.sum(dtype=np.uint64))
    for name, digest in bindings.items():
        assert sha(getattr(args, name)) == digest, 'Source input changed'
    payloads = {name: dict(bytes=(args.out/name).stat().st_size, SHA256=sha(args.out/name)) for name in ['scores.npy', 'detected.npy', 'total-counts.npy']}
    result = dict(status='passed', cells=protocol['cells'], features=len(ids), entries=int(ptr[-1]),
                  originalUMIs=umi_total, emptyLibraries=empty, programs=resolved, bindings=bindings,
                  payloads=payloads, qualityReferenceSHA256=protocol['qualityReferenceSHA256'],
                  everyCellLibraryTotalMatchesIndependentSourceQC=True, maximumSparseProductError=error,
                  blocks=blocks, maximumBlockRows=max_rows, maximumBlockEntries=max_entries,
                  preparerSHA256=sha(Path(__file__)), seconds=time.monotonic()-began,
                  maximumResidentBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system() == 'Darwin' else 1024),
                  scope='Independent measured program reference over all original counts, using exact original symbol mapping and declared missing coverage. Empty-library scores remain NaN in arrays. No native million-cell program publication, pathway activity, cell label, or integration-preservation claim.')
    with (args.out/'reference.json').open('x') as f:
        json.dump(result, f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')
    print(json.dumps({k: result[k] for k in ['status', 'cells', 'entries', 'originalUMIs', 'emptyLibraries', 'maximumSparseProductError', 'seconds', 'maximumResidentBytes']}))


if __name__ == '__main__':
    main()
