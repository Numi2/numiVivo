#!/usr/bin/env python3
"""Prepare all frozen paired experimental comparisons from verified source aggregates."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import h5py
import numpy as np
from anndata.io import read_elem
from acquire import digest, one


def write_json(path, data):
    path.write_text(json.dumps(data, sort_keys=True, indent=2) + '\n')

def table(path, header, rows):
    with path.open('w') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    tools = Path(__file__).parent
    design = json.loads((root / 'design.json').read_text())
    assert design['protocolSHA256'] == digest(tools / 'PROTOCOL.md')
    verification = json.loads((root / 'anndata-verification.json').read_text())
    assert verification['status'] == 'passed'
    assert digest(root / 'hirisa.h5ad') == verification['sourceSHA256']
    source = json.loads((root / 'source-receipt.json').read_text())
    mapping = json.loads((root / 'stream-plan.json').read_text())['mapping']
    by_sample = {s['accession']: s for s in design['samples']}
    native_samples = {s['id']: s for s in mapping['samples']}
    directory = root / 'inference-inputs'
    directory.mkdir(exist_ok=False)
    cohorts = directory / 'cohorts';cohorts.mkdir()
    requests = directory / 'requests';requests.mkdir()
    references = directory / 'reference-inputs';references.mkdir()
    ids = [s['accession'] for s in source['files']]
    ranges = {};offset = 0;all_counts = []
    for sample_id in ids:
        audit = json.loads((root / 'source-audit' / (sample_id + '.json')).read_text())
        reference = root / 'source-audit' / audit['referencePath']
        assert digest(reference) == audit['referenceSHA256']
        with np.load(reference) as counts:
            all_counts.append(counts['geneCounts'].copy())
        ranges[sample_id] = (offset, offset + audit['cells'])
        offset += audit['cells']
    assert offset == verification['cells']
    all_counts = np.array(all_counts, dtype=np.uint64)
    cases = []
    with h5py.File(root / 'hirisa.h5ad', 'r') as h:
        var = read_elem(h['var'])
        features = [{'id': str(g), 'name': str(n), 'mitochondrial': g in mapping['mitochondrialFeatureIDs']}
                    for g, n in zip(var.index, var['name'])]
        feature_ids = [f['id'] for f in features]
        assert all_counts.shape == (131, len(features))
        np.savez_compressed(directory / 'counts.npz', counts=all_counts,
                            sampleIDs=np.array(ids), featureIDs=np.array(feature_ids))
        for index, comparison in enumerate(design['comparisons']):
            cohort_id = f'{index + 1:02d}'
            population, treatment = comparison['population'], comparison['treatment']
            selected = []
            for pair in sorted(comparison['pairs'], key=lambda p: p['donor']):
                selected += [pair['control'], pair['treated']]
            assert len(set(selected)) == len(selected)
            assert len(selected) == (8 if population == 'Bcell' and treatment == 'IFNg' else 10)
            cells = [];groups = [];source_ranges = []
            for sample_id in selected:
                sample = by_sample[sample_id]
                left, right = ranges[sample_id]
                uuids = h['obs/_index'].asstr()[left:right]
                assert set(h['obs/geo_accession'].asstr()[left:right]) == {sample_id}
                local_start = len(cells)
                cells.extend({'barcode': str(u), 'sampleID': sample_id, 'group': population} for u in uuids)
                groups.append({'biologicalReplicateID': one(sample, 'subject id'),
                               'donorID': one(sample, 'subject id'), 'condition': one(sample, 'treatment'),
                               'organism': 'NCBITaxon:9606', 'cellGroup': population,
                               'sampleIDs': [sample_id], 'batchIDs': [one(sample, 'batch pool')],
                               'sourceCellIndices': list(range(local_start, len(cells)))})
                source_ranges.append({'sampleID': sample_id, 'start': left, 'end': right,
                                      'localStart': local_start, 'localEnd': len(cells)})
            y = all_counts[[ids.index(s) for s in selected]]
            row_offsets = [0];columns = [];values = []
            for row in y:
                indices = np.flatnonzero(row)
                columns.extend(indices.tolist());values.extend(row[indices].tolist());row_offsets.append(len(values))
            report = {'metadata': {'id': 'hirisa-inference-' + cohort_id, 'evidence': 'measured',
                        'sourceDescription': 'Python-prepared exact source-audit counts for frozen donor/batch pairs; every selected library cell retained; full-source native ingestion pending.',
                        'countUnit': 'umiCount', 'samples': [native_samples[s] for s in selected],
                        'features': features, 'cells': cells},
                      'pseudobulk': {'method': 'independent-source-CSC-library-aggregation-v1',
                        'countUnit': 'umiCount', 'groups': groups, 'featureIDs': feature_ids,
                        'matrix': {'cellCount': len(groups), 'featureCount': len(features),
                                   'rowOffsets': row_offsets, 'featureIndices': columns, 'counts': values}}}
            report_path = cohorts / (cohort_id + '.json.gz')
            report_path.write_bytes(gzip.compress(json.dumps(report, sort_keys=True, separators=(',', ':')).encode(), mtime=0))
            donors = sorted({g['donorID'] for g in groups})
            names = ['intercept', 'treatment-minus-control'] + ['donor:' + d for d in donors[1:]]
            x = np.array([[1, int(g['condition'] == treatment)] + [int(g['donorID'] == d) for d in donors[1:]] for g in groups])
            assert np.linalg.matrix_rank(x) == x.shape[1] and len(groups) - x.shape[1] == len(donors) - 1
            reference_genes = np.flatnonzero(np.all(y > 0, axis=0))
            assert len(reference_genes) >= 10
            logged = np.log(y[:, reference_genes].astype(float))
            ratios = np.exp(logged - logged.mean(axis=0, keepdims=True))
            log_factors = np.log(np.median(ratios, axis=1));factors = np.exp(log_factors - log_factors.mean())
            libraries = y.sum(axis=1);eligible = (y.sum(axis=0) >= 10) & ((y > 0).sum(axis=0) >= 3)
            contrast = [0, 1] + [0] * (len(names) - 2)
            ref = references / cohort_id;ref.mkdir()
            table(ref / 'counts.tsv', ['featureID'] + selected, ([f] + row.tolist() for f, row in zip(feature_ids, y.T)))
            table(ref / 'design.tsv', ['sampleID'] + names, ([s] + row.tolist() for s, row in zip(selected, x)))
            table(ref / 'samples.tsv', ['sampleID', 'donor', 'condition', 'batch', 'libraryCounts', 'sizeFactor'],
                  ([s, g['donorID'], g['condition'], g['batchIDs'][0], int(l), float(f)] for s, g, l, f in zip(selected, groups, libraries, factors)))
            meta = {'cohort': cohort_id, 'population': population, 'treatment': treatment,
                    'sourceCountsSHA256': digest(directory / 'counts.npz'), 'sourceReportSHA256': digest(report_path),
                    'protocolSHA256': design['protocolSHA256'], 'executionProtocolSHA256': digest(tools / 'INFERENCE_EXECUTION.md'),
                    'sourceH5ADSHA256': verification['sourceSHA256'], 'minimumFeatureCounts': 10,
                    'minimumExpressingPseudobulks': 3, 'eligibleFeatures': int(eligible.sum()), 'attemptedFeatures': len(features),
                    'contrast': contrast, 'columnNames': names, 'design': x.tolist(), 'sourcePseudobulkIndices': list(range(len(groups))),
                    'sampleIDs': selected, 'donorIDs': [g['donorID'] for g in groups], 'sizeFactors': factors.tolist(),
                    'libraryCounts': libraries.tolist(), 'referenceFeatureIndices': reference_genes.tolist(),
                    'sourceObservationRanges': source_ranges,
                    'normalization': 'Independent reconstruction of native all-positive median-ratio factors; compare against native outputs',
                    'nativeRawAggregationQualified': False}
            meta['files'] = {name: digest(ref / name) for name in ['counts.tsv', 'design.tsv', 'samples.tsv']}
            write_json(ref / 'input.json', meta)
            for method, code in [('wald', None), ('lrt', 'likelihoodRatio'), ('ql', 'quasiLikelihoodAdjusted')]:
                options = {'trend': 'gammaParametric', 'minimumTrendGenes': 20, 'minimumPriorVariance': .25, 'outlierStandardDeviations': 2}
                if code:options['testMethod'] = code
                request = {'id': 'hirisa-' + cohort_id + '-' + method, 'model': 'negativeBinomial',
                           'controlCondition': 'none', 'treatmentCondition': treatment, 'cellGroup': population,
                           'design': 'pairedDonors', 'sizeFactors': 'medianRatio', 'adjustForBatch': False,
                           'minimumCellsPerPseudobulk': 10, 'minimumReplicatesPerCondition': 3,
                           'minimumFeatureCounts': 10, 'minimumExpressingPseudobulks': 3,
                           'negativeBinomialOptions': options, 'includedDonorIDs': donors}
                request_path = requests / (cohort_id + '-' + method + '.json');write_json(request_path, request)
                cases.append({'cohort': cohort_id, 'method': method, 'path': request_path.name,
                              'requestSHA256': digest(request_path), 'sourceReport': str(report_path.relative_to(directory)),
                              'sourceReportSHA256': digest(report_path)})
            print(cohort_id, population, treatment, 'donors', len(donors), 'cells', len(cells), 'eligible', int(eligible.sum()), flush=True)
    cases.sort(key=lambda c: (['wald', 'lrt', 'ql'].index(c['method']), c['cohort']))
    write_json(requests / 'manifest.json', cases)
    write_json(directory / 'receipt.json', {'cases': 48, 'cohorts': 16, 'sourceH5ADSHA256': verification['sourceSHA256'],
               'protocolSHA256': design['protocolSHA256'], 'executionProtocolSHA256': digest(tools / 'INFERENCE_EXECUTION.md'),
               'preparerSHA256': digest(Path(__file__)), 'countsSHA256': digest(directory / 'counts.npz'),
               'requestManifestSHA256': digest(requests / 'manifest.json'), 'nativeRawAggregationQualified': False})

if __name__ == '__main__':
    main()
