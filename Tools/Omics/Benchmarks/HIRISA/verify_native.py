#!/usr/bin/env python3
"""Stream-check every native HIRISA cell and link all inference input counts."""
import argparse
import gzip
import json
import math
from importlib.metadata import version
from pathlib import Path

from anndata.io import read_elem
import h5py
import ijson
import numpy as np
from acquire import digest


def items(path, prefix):
    with path.open('rb') as stream:
        yield from ijson.items(stream, prefix)


def compact_fields(path):
    """Parse metadata dictionaries and aggregates while skipping cell/QC objects."""
    wanted = {'schemaVersion', 'canonicalNonzeros', 'metadata.id', 'metadata.countUnit',
              'metadata.evidence', 'metadata.sourceDescription', 'metadata.samples',
              'metadata.features', 'pseudobulk', 'sourceObservationIndices', 'sourceCellCount'}
    result = {};active = None;builder = None
    with path.open('rb') as stream:
        for prefix, event, value in ijson.parse(stream):
            if active is None:
                if prefix not in wanted or event == 'map_key':
                    continue
                assert prefix not in result, ('duplicate report field', prefix)
                if prefix == 'sourceObservationIndices':
                    assert event == 'null', 'This audit requires unselected whole-source ingestion'
                if event in ('start_map', 'start_array'):
                    active = prefix;builder = ijson.common.ObjectBuilder()
                else:
                    result[prefix] = value
                    continue
            builder.event(event, value)
            if prefix == active and event in ('end_map', 'end_array'):
                result[active] = builder.value;active = builder = None
    assert active is None
    return result


def fingerprint(value):
    return bytes(value['bytes']).hex()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    p.add_argument('--bundle', required=True, type=Path,
                   help='Completed native report/plan/receipt; replay separately verifies its retained source')
    a = p.parse_args();root = a.root;bundle = a.bundle
    report = bundle / 'report.json'
    receipt = json.loads((bundle / 'receipt.json').read_text())
    prepared = json.loads((root / 'prepared-receipt.json').read_text())
    qc = json.loads((root / 'native-qc-reference/receipt.json').read_text())
    expected_plan = json.loads((root / 'stream-plan.json').read_text())
    assert json.loads((bundle / 'plan.json').read_text()) == expected_plan
    assert receipt['schemaVersion'] == 1
    assert fingerprint(receipt['source']) == prepared['sourceSHA256'] == qc['sourceSHA256']
    assert digest(root / 'hirisa.h5ad') == prepared['sourceSHA256']
    assert fingerprint(receipt['plan']) == digest(bundle / 'plan.json')
    assert fingerprint(receipt['report']) == digest(report)
    assert qc['status'] == 'passed' and qc['cells'] == prepared['cells']
    fields = compact_fields(report)
    assert fields['schemaVersion'] == 1
    assert fields['canonicalNonzeros'] == prepared['storedEntries']
    assert fields.get('sourceCellCount') in (None, prepared['cells'])
    mapping = expected_plan['mapping']
    for key in ['id', 'countUnit', 'evidence', 'sourceDescription']:
        assert fields['metadata.' + key] == mapping[key]
    assert fields['metadata.samples'] == mapping['samples']
    samples = {s['id']: s for s in mapping['samples']}
    with h5py.File(root / 'hirisa.h5ad', 'r') as h:
        var = read_elem(h['var'])
        features = [{'id': str(g), 'name': str(n), 'mitochondrial': g in mapping['mitochondrialFeatureIDs']}
                    for g, n in zip(var.index, var['name'])]
        assert fields['metadata.features'] == features
        cells = items(report, 'metadata.cells.item')
        qualities = items(report, 'quality.item')
        cell_count = total_umis = detected_total = 0
        for lib in qc['libraries']:
            accession = lib['accession'];left, right = lib['start'], lib['end']
            assert left == cell_count
            uuids = h['obs/_index'].asstr()[left:right]
            assert set(h['obs/geo_accession'].asstr()[left:right]) == {accession}
            ref = root / 'native-qc-reference' / (accession + '.npz')
            assert digest(ref) == lib['referenceSHA256']
            with np.load(ref) as expected:
                for row, uuid in enumerate(uuids):
                    assert next(cells) == {'barcode': str(uuid), 'sampleID': accession, 'group': accession}
                    got = next(qualities)
                    assert got['barcode'] == uuid and got['sampleID'] == accession
                    for field in ['totalCounts', 'detectedFeatures', 'mitochondrialCounts']:
                        assert got[field] == int(expected[field][row]), (accession, row, field)
                    assert got['mitochondrialFeatureCount'] == qc['mitochondrialFeatureCount']
                    fraction = got.get('mitochondrialFraction')
                    if got['totalCounts'] and qc['mitochondrialFeatureCount']:
                        assert fraction is not None and math.isclose(float(fraction),
                            got['mitochondrialCounts'] / got['totalCounts'], rel_tol=1e-14, abs_tol=1e-16)
                    else:
                        assert fraction is None
                    total_umis += got['totalCounts'];detected_total += got['detectedFeatures']
            cell_count = right
            print(accession, 'native cells and QC verified', cell_count, flush=True)
        assert cell_count == prepared['cells'] and next(cells, None) is None and next(qualities, None) is None
    assert detected_total == prepared['storedEntries']
    # Only library-by-gene aggregates are resident, never a cell-by-gene matrix.
    bulk = fields['pseudobulk']
    feature_ids = [f['id'] for f in features]
    assert bulk['featureIDs'] == feature_ids and bulk['countUnit'] == 'umiCount'
    matrix = bulk['matrix'];groups = bulk['groups']
    assert matrix['cellCount'] == len(groups) == len(qc['libraries']) == 131
    assert matrix['featureCount'] == len(feature_ids)
    offsets = matrix['rowOffsets'];indices = matrix['featureIndices'];counts = matrix['counts']
    assert len(offsets) == len(groups) + 1 and offsets[0] == 0 and offsets[-1] == len(indices) == len(counts)
    ranges = {lib['accession']: lib for lib in qc['libraries']}
    inference_receipt = json.loads((root / 'inference-inputs/receipt.json').read_text())
    counts_path = root / 'inference-inputs/counts.npz'
    assert digest(counts_path) == inference_receipt['countsSHA256']
    assert inference_receipt['sourceH5ADSHA256'] == prepared['sourceSHA256']
    seen = set();aggregate_total = 0
    with np.load(counts_path) as inference:
        assert inference['featureIDs'].tolist() == feature_ids
        input_ids = inference['sampleIDs'].tolist()
        for row, group in enumerate(groups):
            assert len(group['sampleIDs']) == 1
            accession = group['sampleIDs'][0]
            assert accession in ranges and accession not in seen;seen.add(accession)
            sample = samples[accession];lib = ranges[accession]
            for key in ['biologicalReplicateID', 'donorID', 'condition', 'organism']:
                assert group[key] == sample[key]
            assert group['cellGroup'] == accession and group['batchIDs'] == [sample['batchID']]
            assert group['sourceCellIndices'] == list(range(lib['start'], lib['end']))
            start, end = offsets[row:row + 2];assert 0 <= start <= end <= len(counts)
            cols = np.array(indices[start:end], dtype=np.int64)
            assert np.all((cols >= 0) & (cols < len(features))) and np.all(np.diff(cols) > 0)
            values = counts[start:end];assert all(type(v) is int and v > 0 for v in values)
            dense = np.zeros(len(features), dtype=np.uint64);dense[cols] = values
            audit = json.loads((root / 'source-audit' / (accession + '.json')).read_text())
            original = root / 'source-audit' / audit['referencePath']
            assert digest(original) == audit['referenceSHA256'] == lib['sourceAuditSHA256']
            with np.load(original) as expected:
                assert np.array_equal(dense, expected['geneCounts'])
            assert np.array_equal(dense, inference['counts'][input_ids.index(accession)])
            aggregate_total += sum(values)
    assert seen == set(ranges) and aggregate_total == total_umis
    # Link the materialized cohort inputs, not just their shared count archive.
    manifest_path = root / 'inference-inputs/requests/manifest.json'
    assert digest(manifest_path) == inference_receipt['requestManifestSHA256']
    manifest = json.loads(manifest_path.read_text())
    assert len(manifest) == 48
    with h5py.File(root / 'hirisa.h5ad', 'r') as h, np.load(counts_path) as inference:
        for cohort in sorted({case['cohort'] for case in manifest}):
            cases = [case for case in manifest if case['cohort'] == cohort]
            assert sorted(case['method'] for case in cases) == ['lrt', 'ql', 'wald']
            meta = json.loads((root / 'inference-inputs/reference-inputs' / cohort / 'input.json').read_text())
            assert meta['sourceH5ADSHA256'] == prepared['sourceSHA256']
            assert meta['sourceCountsSHA256'] == inference_receipt['countsSHA256']
            assert meta['protocolSHA256'] == inference_receipt['protocolSHA256']
            source_report = root / 'inference-inputs' / cases[0]['sourceReport']
            assert digest(source_report) == meta['sourceReportSHA256']
            with gzip.open(source_report, 'rt') as stream:
                cohort_report = json.load(stream)
            assert cohort_report['metadata']['evidence'] == 'measured'
            assert cohort_report['metadata']['countUnit'] == 'umiCount'
            assert cohort_report['metadata']['features'] == features
            assert cohort_report['metadata']['samples'] == [samples[s] for s in meta['sampleIDs']]
            cb = cohort_report['pseudobulk'];cm = cb['matrix']
            assert cb['countUnit'] == 'umiCount'
            assert cb['featureIDs'] == feature_ids
            assert len(cb['groups']) == len(meta['sampleIDs']) == cm['cellCount']
            assert cm['featureCount'] == len(features) and cm['rowOffsets'][0] == 0
            assert cm['rowOffsets'][-1] == len(cm['counts']) == len(cm['featureIndices'])
            local_offset = 0
            for row, (accession, declared_range, group) in enumerate(zip(
                    meta['sampleIDs'], meta['sourceObservationRanges'], cb['groups'], strict=True)):
                lib = ranges[accession];n = lib['end'] - lib['start']
                assert declared_range == {'sampleID': accession, 'start': lib['start'], 'end': lib['end'],
                                          'localStart': local_offset, 'localEnd': local_offset + n}
                uuids = h['obs/_index'].asstr()[lib['start']:lib['end']]
                expected_cells = [{'barcode': str(u), 'sampleID': accession, 'group': meta['population']} for u in uuids]
                assert cohort_report['metadata']['cells'][local_offset:local_offset + n] == expected_cells
                expected_group = {key: samples[accession][key] for key in
                                  ['biologicalReplicateID', 'donorID', 'condition', 'organism']}
                expected_group.update(cellGroup=meta['population'], sampleIDs=[accession],
                                      batchIDs=[samples[accession]['batchID']],
                                      sourceCellIndices=list(range(local_offset, local_offset + n)))
                assert group == expected_group
                start, end = cm['rowOffsets'][row:row + 2]
                assert 0 <= start <= end <= len(cm['counts'])
                cols = np.asarray(cm['featureIndices'][start:end], dtype=np.int64)
                assert np.all((cols >= 0) & (cols < len(features))) and np.all(np.diff(cols) > 0)
                assert all(type(v) is int and v > 0 for v in cm['counts'][start:end])
                dense = np.zeros(len(features), dtype=np.uint64);dense[cols] = cm['counts'][start:end]
                assert np.array_equal(dense, inference['counts'][input_ids.index(accession)])
                local_offset += n
            assert local_offset == len(cohort_report['metadata']['cells'])
            for case in cases:
                assert case['sourceReportSHA256'] == meta['sourceReportSHA256']
                assert case['sourceReport'] == cases[0]['sourceReport']
                assert digest(root / 'inference-inputs/requests' / case['path']) == case['requestSHA256']
    assert len({case['cohort'] for case in manifest}) == 16
    assert digest(report) == fingerprint(receipt['report']), 'Report changed during independent verification'
    result = {'status': 'passed', 'verifierSHA256': digest(Path(__file__)), 'ijsonVersion': version('ijson'),
              'sourceSHA256': prepared['sourceSHA256'], 'reportSHA256': digest(report),
              'receiptSHA256': digest(bundle / 'receipt.json'), 'qcReferenceReceiptSHA256': digest(root / 'native-qc-reference/receipt.json'),
              'inferenceCountsSHA256': digest(counts_path), 'cells': cell_count,
              'features': len(features), 'canonicalNonzeros': detected_total, 'UMIs': total_umis,
              'sourceLibrariesVerified': len(seen), 'inferenceCohortsVerified': 16,
              'inferenceRequestsVerified': 48, 'allCellIdentityAndQCVerified': True,
              'nativeRawIngestionLinked': True, 'nativeReplayVerified': False,
              'scope': 'Independent complete-source native ingestion and exact inference-count linkage; replay is a separate native gate. No biological calibration or prediction claim.'}
    (root / 'native-independent-verification.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
