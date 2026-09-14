"""Check frozen counts, identities and conditional NB fit arithmetic; no truth scoring."""
import collections
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(sys.argv[1])
reports = []
for receipt in json.loads((root / 'native-freeze.json').read_text())['records']:
    site = receipt['site']
    raw = gzip.decompress((root / 'native' / (site + '.json.gz')).read_bytes())
    assert hashlib.sha256(raw).hexdigest() == receipt['outputSHA256']
    obj = json.loads(raw)
    assert obj.get('error') is None and obj.get('result') is not None
    result = obj['result']; design = result['design']
    x = np.asarray(design['rows']); contrast = np.asarray(design['contrast'])
    sf = np.asarray(design['sizeFactorValues'])
    counts = pd.read_csv(root / 'counts' / site / 'counts.tsv.gz', sep='\t').to_numpy(dtype=np.uint64)
    counts = counts[:, design['sourcePseudobulkIndices']]
    assert counts.sum(axis=0).tolist() == design['libraryCounts']
    assert len(result['features']) == len(result['negativeBinomial']['features']) == 25794
    maxima = collections.defaultdict(float); fits = collections.Counter()
    for i, (feature, diagnostic) in enumerate(zip(result['features'], result['negativeBinomial']['features'])):
        assert feature['featureIndex'] == diagnostic['featureIndex'] == i
        assert feature['featureID'] == 'SEQC:RefSeq:row:' + str(i)
        assert feature['totalCounts'] == int(counts[i].sum())
        y = counts[i].astype(float)
        for key in ('finalFit', 'effectShrinkageFit'):
            fit = diagnostic.get(key)
            if fit is None:
                continue
            beta = np.asarray(fit['coefficients']); mu = np.asarray(fit['means']); alpha = fit['dispersion']
            np.testing.assert_allclose(mu, np.exp(x @ beta + np.log(sf)), rtol=2e-12)
            if key == 'finalFit':
                weights = mu / (1 + alpha * mu)
                information = x.T @ (x * weights[:, None])
                se = fit.get('standardError')
            else:
                precision = 1 / fit['priorStandardDeviation'] ** 2
                weights = mu * (1 + alpha * y) / (1 + alpha * mu) ** 2
                information = x.T @ (x * weights[:, None]) + precision * np.outer(contrast, contrast)
                se = fit.get('posteriorStandardDeviation')
            if se is not None:
                expected = np.sqrt(contrast @ np.linalg.solve(information, contrast))
                error = abs(expected - se) / max(1, abs(se))
                assert error < 2e-9, (site, i, key, error)
                maxima[key + 'StandardError'] = max(maxima[key + 'StandardError'], error)
            fits[key] += 1
    report = dict(site=site, features=25794, fits=dict(fits), maximumErrors=dict(maxima),
                  statuses=dict(collections.Counter(f['status'] for f in result['features'])))
    reports.append(report)
    print(json.dumps(report), flush=True)
(root / 'native-checks.json').write_text(json.dumps(reports, indent=2) + '\n')
