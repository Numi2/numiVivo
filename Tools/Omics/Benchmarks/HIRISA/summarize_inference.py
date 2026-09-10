#!/usr/bin/env python3
"""Report full and jointly tested experimental families without calibration claims."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from acquire import digest


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind='stable')
    adjusted = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values));result[order] = np.minimum(1, adjusted)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs', required=True, type=Path)
    p.add_argument('--available', action='store_true')
    p.add_argument('--native-verification', type=Path, help='Completed independent full-source native linkage receipt')
    a = p.parse_args();root = a.inputs
    linkage = None
    if a.native_verification is not None:
        linkage = json.loads(a.native_verification.read_text())
        prepared = json.loads((root / 'receipt.json').read_text())
        assert linkage['status'] == 'passed' and linkage['nativeRawIngestionLinked']
        assert linkage['inferenceCountsSHA256'] == prepared['countsSHA256'] == digest(root / 'counts.npz')
        assert linkage['sourceSHA256'] == prepared['sourceH5ADSHA256']
        assert linkage['inferenceCohortsVerified'] == 16 and linkage['inferenceRequestsVerified'] == 48
    method_rows = [];comparisons = [];missing = [];checks = [];joint_families = []
    for number in range(1, 17):
        case = f'{number:02d}';meta = json.loads((root / 'reference-inputs' / case / 'input.json').read_text())
        tables = {};base = {'cohort': case, 'population': meta['population'], 'treatment': meta['treatment'],
                          'donors': len(set(meta['donorIDs'])), 'eligibleFeatures': meta['eligibleFeatures']}
        for method in ['wald', 'lrt', 'ql']:
            directory = root / 'native' / (case + '-' + method)
            if not (directory / 'receipt.json').exists():
                missing.append(case + '-' + method);continue
            receipt = json.loads((directory / 'receipt.json').read_text())
            assert digest(directory / 'output.json.gz') == receipt['outputSHA256']
            out = json.loads(gzip.decompress((directory / 'output.json.gz').read_bytes()))
            result = out.get('result')
            if result is None:
                method_rows.append({**base, 'method': 'native-' + method, 'status': 'failed', 'error': out.get('error')});continue
            frame = pd.DataFrame(result['features']).set_index('featureID')
            eligible = frame[frame.status.eq('tested')]
            assert eligible.pValue.notna().all()
            tables['native-' + method] = eligible
            record = {**base, 'method': 'native-' + method, 'status': 'completed',
                      'tested': len(eligible), 'BH05Calls': int(eligible.adjustedPValue.le(.05).sum()),
                      'featureStatuses': dict(Counter(frame.status)), 'seconds': receipt['seconds']}
            method_rows.append(record)
            if (directory / 'model-check.json').exists():
                check = json.loads((directory / 'model-check.json').read_text())
                checks.append({'cohort': case, 'method': method, 'status': check['status'],
                               'failedChecks': [c for c in check['checks'] if not c['passed']]})
        directory = root / 'reference' / case
        if not (root / 'reference' / (case + '-receipt.json')).exists():
            missing.append(case + '-reference');continue
        reference_receipt = json.loads((root / 'reference' / (case + '-receipt.json')).read_text())
        for name, sha in reference_receipt['files'].items():assert digest(directory / name) == sha
        runs = json.loads((directory / 'runs.json').read_text())
        for suffix in ['edgeR-QL', 'limma-voom', 'DESeq2']:
            name = 'native_size_factors-' + suffix;run = runs[name]
            if run['status'] != 'completed':
                method_rows.append({**base, 'method': suffix, 'status': 'failed', 'error': run['error']});continue
            frame = pd.read_csv(directory / (name + '.tsv'), sep='\t').set_index('featureID')
            eligible = frame[frame.pValue.notna()]
            tables[suffix] = eligible
            method_rows.append({**base, 'method': suffix, 'status': 'completed', 'tested': len(eligible),
                               'BH05Calls': int(eligible.adjustedPValue.le(.05).sum()),
                               'nonfinitePValues': int(frame.pValue.isna().sum()), 'warnings': run['warnings'], 'seconds': run['seconds']})
        if len(tables) == 6:
            common = sorted(set.intersection(*(set(t.index) for t in tables.values())))
            joint_families.append({**base, 'jointlyTestedAcrossSixMethods': len(common),
                                   'BH05Calls': {name: int(np.count_nonzero(bh(t.loc[common, 'pValue']) <= .05))
                                                 for name, t in tables.items()}})
        for native in ['native-wald', 'native-lrt', 'native-ql']:
            if native not in tables:continue
            for reference in ['edgeR-QL', 'limma-voom', 'DESeq2']:
                if reference not in tables:continue
                ids = sorted(set(tables[native].index) & set(tables[reference].index))
                n, r = tables[native].loc[ids], tables[reference].loc[ids]
                finite = np.isfinite(n.log2FoldChange) & np.isfinite(r.log2FoldChange)
                correlation = float(spearmanr(n.loc[finite, 'log2FoldChange'], r.loc[finite, 'log2FoldChange']).statistic)
                calls_n = set(n.index[bh(n.pValue) <= .05]);calls_r = set(r.index[bh(r.pValue) <= .05])
                comparisons.append({**base, 'native': native, 'reference': reference, 'jointlyTested': len(ids),
                                    'effectSpearman': correlation if np.isfinite(correlation) else None,
                                    'nativeJointBH05': len(calls_n), 'referenceJointBH05': len(calls_r),
                                    'jointBH05Intersection': len(calls_n & calls_r),
                                    'jointBH05Union': len(calls_n | calls_r)})
    result = {'allCasesTerminal': not missing, 'missingCases': missing, 'methods': method_rows,
              'comparisons': comparisons, 'sixMethodJointFamilies': joint_families, 'numericalChecks': checks, 'summarizerSHA256': digest(Path(__file__)),
              'nativeRawIngestionLinked': linkage is not None,
              'nativeVerificationSHA256': digest(a.native_verification) if linkage is not None else None,
              'biologicalCalibrationClaim': False,
              'scope': 'Frozen experimental counts/designs; full method families and pairwise jointly tested BH families. Agreement is not ground truth, FDR, power or coverage.'}
    destination = root / ('inference-summary-available.json' if a.available else 'inference-summary.json')
    destination.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n')
    assert a.available or not missing, 'Some frozen cases are not terminal'
    print(json.dumps({'terminalMethodRows': len(method_rows), 'pairwiseComparisons': len(comparisons), 'missingCases': len(missing)}))

if __name__ == '__main__':
    main()
