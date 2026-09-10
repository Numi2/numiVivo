#!/usr/bin/env python3
"""Reconstruct the unresolved Adamson guide inventory without expression data.

This is a metadata audit, not a control caller or prediction launcher. It reads
the preserved GEO and cohort archives, verifies their member hashes, and retains
candidate assignments as unverified even when the guide name resembles a gene.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def audit():
    inputs = []

    def read(archive, name):
        root = ROOT / 'evidence' / archive
        manifest_bytes = (root / 'manifest.json').read_bytes()
        manifest = json.loads(manifest_bytes)
        record = next(r for r in manifest.get('entries', manifest.get('files', []))
                      if r['logicalPath'] == name)
        raw = (root / record['path']).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == record['sha256'], name
        assert len(raw) == record['bytes'], name
        if record['path'] != name and record['path'].endswith('.gz'):
            raw = gzip.decompress(raw)
        assert hashlib.sha256(raw).hexdigest() == record['logicalSHA256'], name
        assert len(raw) == record['logicalBytes'], name
        inputs.append(dict(archive=archive, path=name,
                           sha256=record['logicalSHA256'],
                           manifestSHA256=hashlib.sha256(manifest_bytes).hexdigest()))
        return raw

    geo = '2026-09-10-geo-restoration'
    cohort_archive = '2026-09-10-author-cohort'
    barcodes = gzip.decompress(read(geo, 'source/GSM2406681_10X010_barcodes.tsv.gz')).decode().splitlines()
    reader = csv.DictReader(io.StringIO(gzip.decompress(read(
        geo, 'source/GSM2406681_10X010_cell_identities.csv.gz')).decode()))
    index = reader.fieldnames[0]
    records = list(reader)
    by_barcode = {r[index]: r for r in records}
    assert len(by_barcode) == len(records), 'Duplicate original guide record'
    assert len(set(barcodes)) == len(barcodes), 'Duplicate original barcode'
    cohort = json.loads(read(cohort_archive, 'cohort.json'))
    prior = json.loads(read(cohort_archive, 'audit.json'))
    candidates = json.loads(read(geo, 'descriptors/candidates.json'))['candidates']
    coverage = {r['target']: r for r in json.loads(read(geo, 'descriptors/coverage.json'))}
    source, selected = Counter(), Counter()
    selected_rows, selected_barcodes, excluded_rows = [], [], []
    missing = 0
    for row, barcode in enumerate(barcodes):
        record = by_barcode.get(barcode)
        if record is None:
            missing += 1
            excluded_rows.append(row)
            continue
        guide = record['guide identity']
        source[guide] += 1
        assert record['good coverage'] in ('TRUE', 'FALSE')
        keep = (int(record['number of cells']) == 1
                and record['good coverage'] == 'TRUE' and guide != '*')
        if keep:
            selected[guide] += 1
            selected_rows.append(row)
            selected_barcodes.append(barcode)
        else:
            excluded_rows.append(row)
    assert selected_rows == cohort['selectedRows']
    assert selected_barcodes == cohort['barcodes']
    assert excluded_rows == [r['sourceRow'] for r in cohort['excluded']]
    assert sorted(selected) == cohort['groups']
    assert len(selected_rows) == prior['selectedCells']
    assert len(excluded_rows) == prior['excludedCells']
    assert all(source[r['guide']] == r['sourceCells']
               and selected[r['guide']] == r['selectedCells'] for r in prior['guides'])
    target_by_guide = {}
    for candidate in candidates:
        assert candidate['assignment'] == 'unverified-guide-prefix'
        for guide in candidate['sourceGuides']:
            assert guide not in target_by_guide
            target_by_guide[guide] = candidate['target']
    retained = sorted({target_by_guide[g] for g in selected if g in target_by_guide})
    assert retained == prior['retainedGeneLikePrefixes']
    non_gene = sorted(set(selected) - set(target_by_guide))
    guides = []
    for guide in sorted(source):
        target = target_by_guide.get(guide)
        guides.append(dict(guide=guide, sourceCells=source[guide],
                           selectedCells=selected[guide],
                           excludedCells=source[guide] - selected[guide],
                           candidateTarget=target,
                           candidateGOStatus=coverage[target]['status'] if target else None,
                           experimentalRole='unverified'))
    return dict(schemaVersion=1, metadataReconstruction='passed',
                inputs=inputs, sourceCells=len(barcodes), missingGuideRecords=missing,
                selectedCells=len(selected_rows), excludedCells=len(excluded_rows),
                selectedGuideGroups=len(selected), selectedCandidateTargets=retained,
                supportedCandidateTargets=[t for t in retained if coverage[t]['status'] == 'supported'],
                selectedNonGeneLabels=non_gene,
                selectedGeneLikeGuideGroups=sum(g in target_by_guide for g in selected),
                guides=guides, controlsVerified=False, experimentalTargetsVerified=False,
                readyForFrozenPrediction=False, expressionOutcomesRead=False,
                fittingPerformed=False,
                pending=['Primary experimental control-to-guide mapping',
                         'Primary guide-to-target roster and reconciliation with the paper',
                         'Freeze verified roles and folds before fitting; freeze predictions before scoring'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    result = audit()
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
    print(json.dumps({k: result[k] for k in ('metadataReconstruction', 'sourceCells',
                     'selectedCells', 'selectedGuideGroups', 'selectedGeneLikeGuideGroups',
                     'selectedNonGeneLabels', 'readyForFrozenPrediction')}))


if __name__ == '__main__':
    main()
