#!/usr/bin/env python3
"""Verify portable native artifacts without altering historical record paths."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file


def verify(root, output):
    manifest = json.loads((root / 'manifest.json').read_text())
    for name, digest in manifest['files'].items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError('Changed package artifact: ' + name)
    output.mkdir(exist_ok=False)
    checks = {}
    for i, model in enumerate(json.loads((root / 'native-models.json').read_text())):
        dest = output / str(i)
        subprocess.run([str(root / 'runtime/numivivo'), 'spatial-response', 'predict',
                        str(root / model['plan']), str(root / model['query']), str(dest),
                        str(root / model['weights'])], check=True)
        a, b = load_file(str(root / model['prediction'])), load_file(str(dest / 'prediction.safetensors'))
        if a.keys() != b.keys() or not all(np.array_equal(v, b[k]) for k, v in a.items()):
            raise ValueError('Native replay mismatch: ' + model['name'])
        checks[model['name']] = 'bit-exact mean and variance'
    result = {'status': 'verified', 'models': checks, 'biologicalPromotion': False,
              'scope': 'same-host portable native replay; not biological or Safari qualification'}
    (output / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    verify(Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve())
