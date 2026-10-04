#!/usr/bin/env python3
"""Inspect training underfit without refitting or changing evaluation criteria."""
import argparse
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file
from wetlab import read, write, require
from intervention_design import bound
from train_intervention_design import native


def run(campaign, output):
    campaign, inputs, reg = bound(campaign)
    require(not output.exists(), 'Diagnostic already exists')
    output.mkdir()
    folder = campaign / 'target-descriptor'
    step = read(folder / 'selection.json')['selected']['step']
    native(reg['binary'], 'predict', folder / 'plan.json', folder / 'training.safetensors',
           output / 'native', folder / 'training' / f'weights-{step}.safetensors')
    q = load_file(str(folder / 'training.safetensors'))
    p = load_file(str(output / 'native/prediction.safetensors'))
    rows = read(inputs / 'training-rows.json')
    observed, predicted = [], []
    keys = sorted({(r['source'], r['context'], r['target']) for r in rows})
    for key in keys:
        indices = [i for i, r in enumerate(rows) if (r['source'], r['context'], r['target']) == key]
        observed.append((q['observed'][indices] - q['context'][indices]).mean(0))
        predicted.append((p['mean'][indices] - q['context'][indices]).mean(0))
    y, delta = np.asarray(observed), np.asarray(predicted)
    result = {'campaign': str(campaign), 'sourceContextTargetGroups': len(keys),
              'observedGroupResponseSD': float(y.std(0).mean()),
              'predictedGroupResponseSD': float(delta.std(0).mean()),
              'groupResponseRMSE': float(np.sqrt(np.mean((y - delta)**2))),
              'noChangeGroupResponseRMSE': float(np.sqrt(np.mean(y**2))),
              'scope': 'Training bags averaged by source/context/target; not independent biological replicates. Post-reveal diagnostic, no refit or new validation.',
              'biologicalPromotion': False}
    write(output / 'diagnostic.json', result)
    print(result)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--campaign', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.campaign, a.output)
