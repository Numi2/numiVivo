#!/usr/bin/env python3
"""Audit the original Muraro numeric table without rounding or selecting cells."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
sha = hashlib.sha256(a.source.read_bytes()).hexdigest()
assert sha == '2f253ffb1f6d54f6bb259e862195f5c20ea3f13e00471b9864ff92773658f513'
features = nonzeros = fractional = invalid = 0
examples = []
with gzip.open(a.source, 'rt') as stream:
    cells = next(csv.reader([stream.readline()], delimiter='\t'))
    assert len(cells) == len(set(cells))
    for line in stream:
        feature, data = line.rstrip('\n').split('\t', 1)
        values = np.fromstring(data, sep='\t')
        assert len(values) == len(cells)
        features += 1
        nonzeros += int(np.count_nonzero(values))
        invalid += int(np.sum(~np.isfinite(values) | (values < 0)))
        mask = np.isfinite(values) & (values != np.floor(values))
        fractional += int(mask.sum())
        if mask.any() and len(examples) < 3:
            index = np.flatnonzero(mask)[0]
            examples.append(dict(feature=feature.strip('"'), cell=cells[index], value=float(values[index])))
assert (len(cells), features, nonzeros, fractional, invalid) == (3072, 19140, 12442034, 12442034, 0)
result = dict(status='ineligible-for-integer-count-analysis', sourceSHA256=sha,
    sourceURL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85241/suppl/GSE85241_cellsystems_dataset_4donors_updated.csv.gz',
    cells=len(cells), features=features, nonzeros=nonzeros, fractionalEntries=fractional, invalidEntries=invalid,
    cellPrefixCountsNotDonorVerification=dict(Counter(cell.split('-')[0] for cell in cells)), fractionalExamples=examples,
    decision='No rounding, inverse correction, cell/feature subset or synthetic fallback. This release cannot qualify the native integer-count path; retaining continuous expression requires a separately typed assay path.')
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
