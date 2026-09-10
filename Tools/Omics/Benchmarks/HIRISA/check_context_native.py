#!/usr/bin/env python3
"""Check every selected native context aggregate against the sparse reference."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import resource
import numpy as np


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('bundle', 'frozen', 'reference', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--lineage', required=True)
    args = parser.parse_args()
    start = time.monotonic()
    frozen = json.loads((args.frozen/'freeze.json').read_text())
    reference = json.loads((args.reference/'check.json').read_text())
    assert sha(args.frozen/'freeze.json') == reference['freezeSHA256']
    assert sha(args.reference/'counts.npz') == reference['countsSHA256']
    assert sha(args.reference/'identities.json') == reference['identitiesSHA256']
    name = args.lineage+'-plan.json'
    assert sha(args.frozen/name) == frozen['files'][name]['sha256']
    plan = json.loads((args.frozen/name).read_text())
    assert json.loads((args.bundle/'plan.json').read_text()) == plan
    report = json.loads((args.bundle/'report.json').read_text())
    rows = plan['cellSelection']['observationIndices']
    assert report['sourceCellCount'] == frozen['sourceCells']
    assert report['sourceObservationIndices'] == rows
    assert len(report['metadata']['cells']) == len(report['quality']) == len(rows)
    identities = json.loads((args.reference/'identities.json').read_text())
    aggregates = json.loads((args.frozen/'aggregates.json').read_text())
    index = {k: i for i, k in enumerate(identities['groups'])}
    values = np.load(args.reference/'counts.npz', allow_pickle=False)
    bulk = report['pseudobulk']
    assert bulk['featureIDs'] == identities['features']
    matrix = bulk['matrix']
    found, members = set(), []
    total = 0
    for i, group in enumerate(bulk['groups']):
        donor, preparation = group['biologicalReplicateID'].rsplit('-', 1)
        role = {'none': 'control', 'IFNa': 'treated'}[group['condition']]
        key = f'{args.lineage}:{preparation}:{donor}:{role}'
        assert key not in found
        found.add(key)
        aggregate = aggregates[key]
        assert group['donorID'] == donor and group['cellGroup'] == aggregate['authorLabel']
        assert group['sampleIDs'] == aggregate['accessions'] and group['batchIDs'] == aggregate['batchPools']
        local_rows = group['sourceCellIndices']
        assert len(local_rows) == aggregate['cells']
        original_rows = np.asarray([rows[j] for j in local_rows], dtype='<u4')
        assert hashlib.sha256(original_rows.tobytes()).hexdigest() == aggregate['originalRowsSHA256']
        for j in local_rows:
            cell = report['metadata']['cells'][j]
            assert cell['sampleID'] in aggregate['accessions'] and cell['group'] == aggregate['authorLabel']
        members.extend(local_rows)
        first, last = matrix['rowOffsets'][i:i+2]
        columns = matrix['featureIndices'][first:last]
        assert columns == sorted(set(columns))
        actual = np.zeros(len(identities['features']), dtype=np.int64)
        actual[columns] = matrix['counts'][first:last]
        np.testing.assert_array_equal(actual, values['counts'][index[key]])
        quality_total = sum(report['quality'][j]['totalCounts'] for j in local_rows)
        assert quality_total == int(values['totals'][index[key]])
        total += quality_total
    assert found == {k for k in aggregates if aggregates[k]['lineage'] == args.lineage}
    assert sorted(members) == list(range(len(rows)))
    result = dict(status='passed', lineage=args.lineage, cells=len(rows), groups=len(found),
                  features=len(identities['features']), totalUMIs=total,
                  everyNativeCountAndOriginalMembershipExact=True,
                  referenceCheckSHA256=sha(args.reference/'check.json'), reportSHA256=sha(args.bundle/'report.json'),
                  scriptSHA256=sha(Path(__file__)), seconds=time.monotonic()-start,
                  maximumRSSBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  predictionFittingPerformed=False)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
