#!/usr/bin/env python3
"""Freeze all prespecified donor-held-out folds without reading expression values."""
import argparse
import json
from pathlib import Path
from acquire import digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', required=True, type=Path)
    a = p.parse_args();root = a.root;here = Path(__file__).resolve().parent
    repo = here.parents[3]
    design = json.loads((root / 'design.json').read_text())
    prepared = json.loads((root / 'prepared-receipt.json').read_text())
    assert design['protocolSHA256'] == digest(here / 'PROTOCOL.md') and prepared['complete']
    folds = []
    for number, comparison in enumerate(design['comparisons'], 1):
        pairs = sorted(comparison['pairs'], key=lambda x: x['donor'])
        assert len(pairs) == (4 if comparison['population'] == 'Bcell' and comparison['treatment'] == 'IFNg' else 5)
        for held_out in pairs:
            training = [pair for pair in pairs if pair['donor'] != held_out['donor']]
            assert len(training) == len(pairs) - 1
            training_ids = [pair[k] for pair in training for k in ['control', 'treated']]
            assert held_out['control'] not in training_ids and held_out['treated'] not in training_ids
            folds.append({'id': f'{number:02d}-' + held_out['donor'], 'cohort': f'{number:02d}',
                          'population': comparison['population'], 'treatment': comparison['treatment'],
                          'trainingPairs': training, 'queryControlAccession': held_out['control'],
                          'scoringTargetAccession': held_out['treated'], 'heldOutDonor': held_out['donor']})
    assert len(folds) == 79 and len({f['id'] for f in folds}) == 79
    owners = ['Sources/NumiVivoKit/Omics/VivoPerturbation.swift',
              'Sources/NumiVivoKit/Omics/VivoPerturbationIO.swift']
    result = {'schemaVersion': 1, 'folds': folds, 'foldCount': len(folds),
              'designSHA256': digest(root / 'design.json'), 'protocolSHA256': digest(here / 'PROTOCOL.md'),
              'predictionProtocolSHA256': digest(here / 'PREDICTION_EXECUTION.md'),
              'sourceH5ADSHA256': prepared['sourceSHA256'], 'freezerSHA256': digest(Path(__file__)),
              'ownerSourceSHA256': {name: digest(repo / name) for name in owners},
              'method': 'paired-donor-log1p-CPM-response-baselines-alpha1-v1',
              'responseModelFittingStarted': False, 'nativeEndToEndQualified': False}
    with (root / 'prediction-folds.json').open('x') as out:
        out.write(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'folds': len(folds), 'manifestSHA256': digest(root / 'prediction-folds.json')}))


if __name__ == '__main__':
    main()
