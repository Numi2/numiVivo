#!/usr/bin/env python3
"""Run the full product on all Norman cells; compare every aggregate and cell QC value."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('binary', 'source', 'prepared', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    checks = []

    def run(name, command, success=True):
        with (args.out / (name + '.log')).open('w') as log:
            result = subprocess.run(['/usr/bin/time', '-l', str(args.binary.resolve()), *map(str, command)],
                                    stdout=log, stderr=subprocess.STDOUT)
        assert result.returncode == (0 if success else 65), (name, result.returncode)
        checks.append(name)

    bundle = args.out / 'bundle'
    run('aggregate', ['singlecell-h5ad-pseudobulk', args.source, '--plan', args.prepared / 'plan.json', '--output', bundle])
    run('verify', ['singlecell-h5ad-pseudobulk-verify', bundle])
    run('reject-overwrite', ['singlecell-h5ad-pseudobulk', args.source, '--plan', args.prepared / 'plan.json', '--output', bundle], False)
    audit = json.loads((args.prepared / 'audit.json').read_text())
    # NpzFile lazily decompresses on every key lookup. Materialize each small
    # reference array once, not once per cell in the exact-value comparisons.
    with np.load(args.prepared / 'reference.npz', allow_pickle=False) as archive:
        reference = {key: archive[key] for key in archive.files}
    report = json.loads((bundle / 'report.json').read_text())
    assert fingerprint(bundle / 'original.h5ad') == audit['source']['sha256']
    assert report['canonicalNonzeros'] == audit['storedEntries']
    assert report['contrasts'] == []
    metadata, bulk = report['metadata'], report['pseudobulk']
    assert bulk['featureIDs'] == reference['feature_ids'].tolist()
    assert [f['id'] for f in metadata['features']] == reference['feature_ids'].tolist()
    assert [f['name'] for f in metadata['features']] == reference['gene_symbols'].tolist()
    assert len(metadata['cells']) == audit['shape'][0] == len(report['quality'])
    for i, (cell, quality) in enumerate(zip(metadata['cells'], report['quality'])):
        assert cell['barcode'] == quality['barcode'] == reference['barcodes'][i]
        assert cell['sampleID'] == quality['sampleID'] == reference['cell_conditions'][i]
        assert cell['group'] == 'K562'
        for key, expected in [('totalCounts', 'totals'), ('detectedFeatures', 'detected'), ('mitochondrialCounts', 'mitochondrial')]:
            assert quality[key] == reference[expected][i], (i, key)
        assert quality['mitochondrialFeatureCount'] == audit['mitochondrialFeatureCount']
        expected_fraction = float(reference['mitochondrial'][i]) / float(reference['totals'][i]) if reference['totals'][i] else None
        assert quality.get('mitochondrialFraction') == expected_fraction
    assert [g['condition'] for g in bulk['groups']] == reference['conditions'].tolist()
    matrix = bulk['matrix']
    assert matrix['cellCount'] == 237 and matrix['featureCount'] == 33694
    assert len(matrix['counts']) == audit['aggregateNonzeros']
    for i, group in enumerate(bulk['groups']):
        assert group['biologicalReplicateID'] == 'K562-pooled-replication-unresolved'
        assert group.get('donorID') is None
        assert group['cellGroup'] == 'K562' and group['organism'] == 'NCBITaxon:9606'
        assert group['sampleIDs'] == [group['condition']]
        assert group['batchIDs'] == ['pooled-eight-gemgroups']
        assert group['sourceCellIndices'] == np.flatnonzero(reference['cell_conditions'] == group['condition']).tolist()
        start, end = matrix['rowOffsets'][i:i + 2]
        expected = reference['counts'][i]
        indices = np.flatnonzero(expected)
        assert matrix['featureIndices'][start:end] == indices.tolist()
        assert matrix['counts'][start:end] == expected[indices].tolist()
    assert not list(args.out.glob('.numivivo-stream-*'))
    result = dict(status='passed', checks=checks, source=audit['source'],
                  binarySHA256=fingerprint(args.binary), reportSHA256=fingerprint(bundle / 'report.json'),
                  exactAggregateValues=audit['aggregateNonzeros'], exactCellQCRecords=audit['shape'][0],
                  sourceEntries=audit['storedEntries'], totalCounts=audit['totalCounts'],
                  limitations='One complete filtered release; no donor replication, DE, prediction, million-cell or GPU qualification.')
    (args.out / 'checks.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
