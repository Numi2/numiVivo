#!/usr/bin/env python3
"""Descriptive target-gene activation, without alias inference or replicate tests."""
import argparse
import json
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with np.load(args.prepared / 'reference.npz', allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    conditions = data['conditions'].tolist()
    symbols = data['gene_symbols'].tolist()
    assert len(set(symbols)) == len(symbols)
    counts = data['counts']
    totals = counts.sum(axis=1)
    assert np.all(totals > 0)
    expression = np.log1p(counts.astype(float) / totals[:, None] * 1e6)
    control = conditions.index('control')
    rows = []
    for name in conditions:
        if name == 'control' or '_' in name:
            continue
        if name not in symbols:
            rows.append(dict(target=name, status='missing-feature'))
            continue
        i, j = conditions.index(name), symbols.index(name)
        rows.append(dict(target=name, featureID=str(data['feature_ids'][j]), status='observed',
                         controlCounts=int(counts[control, j]), perturbedCounts=int(counts[i, j]),
                         targetLog1pCPMChange=float(expression[i, j] - expression[control, j])))
    observed = [row for row in rows if row['status'] == 'observed']
    result = dict(scope='Descriptive target-gene activation in pooled single-target conditions; '
                  'no independent replicate inference or causal validation. '
                  'Full-gene library denominators, log1p CPM at 1e6.',
                  observedTargets=len(observed),
                  missingTargets=[row['target'] for row in rows if row['status'] == 'missing-feature'],
                  positiveChanges=sum(row['targetLog1pCPMChange'] > 0 for row in observed),
                  nonpositiveChanges=[row['target'] for row in observed if row['targetLog1pCPMChange'] <= 0],
                  rows=rows)
    with args.output.open('x') as stream:
        stream.write(json.dumps(result, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'rows'}))


if __name__ == '__main__':
    main()
