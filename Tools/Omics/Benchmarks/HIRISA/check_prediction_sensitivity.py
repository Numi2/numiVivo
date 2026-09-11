#!/usr/bin/env python3
"""Independently check every sensitivity value with NumPy source-score arithmetic."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--result', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        raise ValueError('Output must be new')
    raw = a.result.read_bytes()
    report = json.loads(raw)
    root = Path(__file__).resolve().parent
    lookup = {}
    for identity in report['inputs']:
        directory = root / 'evidence' / identity['archive']
        manifest = json.loads((directory / 'manifest.json').read_bytes())
        parts = sorted((r for r in manifest['records'] if r['sourcePath'] == identity['sourcePath']),
                       key=lambda r: r.get('sourceOffset', 0))
        payload = b''.join(gzip.decompress((directory / r['storedPath']).read_bytes())
                           if r['gzipEncoded'] else (directory / r['storedPath']).read_bytes()
                           for r in parts)
        assert hashlib.sha256(payload).hexdigest() == identity['SHA256']
        study = 'known-treatment' if 'prediction-scores' in identity['sourcePath'] else 'preparation-transfer'
        for fold in json.loads(payload)['folds']:
            contrast = ((fold['population'], fold['treatment']) if study == 'known-treatment'
                        else (fold['kind'], fold['queryPreparation'], fold['lineage']))
            for score in fold['scores']:
                key = study, contrast, score['family'], score['method'], fold['heldOutDonor']
                assert key not in lookup
                lookup[key] = score['responseRMSE']
    checked = 0
    for row in report['comparisons'] + report['transferPenalties']:
        if 'study' in row:
            c = b = row['study'], tuple(row['contrast']), row['family']
            candidate, baseline = row['candidate'], row['baseline']
        else:
            c = 'preparation-transfer', ('cross', row['queryPreparation'], row['lineage']), 'all-source-genes'
            b = 'preparation-transfer', ('within', row['queryPreparation'], row['lineage']), 'all-source-genes'
            candidate = row['candidate'].removeprefix('cross-')
            baseline = row['baseline'].removeprefix('within-')
        donors = sorted(d for s, k, f, m, d in lookup if (s, k, f) == c and m == candidate)
        assert len(donors) == row['donorCount'] and set(donors) == set(row['pairedDifferences'])
        x = np.asarray([[lookup[(*c, candidate, d)], lookup[(*b, baseline, d)]] for d in donors])
        differences = x[:, 0] - x[:, 1]
        # Separate means of the candidate/baseline after boolean donor removal.
        means = np.array([(x[np.arange(len(donors)) != i, 0].mean() -
                           x[np.arange(len(donors)) != i, 1].mean()) for i in range(len(donors))])
        full = x[:, 0].mean() - x[:, 1].mean()
        np.testing.assert_allclose(differences, [row['pairedDifferences'][d] for d in donors], atol=1e-12, rtol=0)
        np.testing.assert_allclose(means, [row['omittedDonorMeans'][d] for d in donors], atol=1e-12, rtol=0)
        np.testing.assert_allclose([full, means.min(), means.max()],
                                   [row['meanDifference'], row['omittedMinimum'], row['omittedMaximum']],
                                   atol=1e-12, rtol=0)
        expected = ('stable-improvement' if full < 0 and np.all(means < 0) else
                    'stable-worse' if full > 0 and np.all(means > 0) else
                    'all-tied' if full == 0 and np.all(means == 0) else 'donor-sensitive')
        assert row['status'] == expected
        assert row['donorWins'] == int((differences < 0).sum())
        assert row['donorTies'] == int((differences == 0).sum())
        checked += 1
    assert len(lookup) == 199 * 8 and checked == 368
    # Hand-calculated examples exercise the interpretation, including a full-mean
    # improvement that reverses when one influential donor score is omitted.
    from prediction_sensitivity import comparison
    for values, expected in [([-3., 1., 1.], 'donor-sensitive'),
                             ([-4., -2., -1.], 'stable-improvement'),
                             ([4., 2., 1.], 'stable-worse'), ([0., 0., 0.], 'all-tied')]:
        assert comparison(dict(zip('abc', values)), 'c', 'b')['status'] == expected
    result = dict(status='passed', sourceScores=len(lookup), comparisons=checked,
                  handCalculatedCases=4, tolerance=1e-12, numpyVersion=np.__version__,
                  diagnosticSHA256=hashlib.sha256(raw).hexdigest(),
                  checkerSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    a.out.write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
