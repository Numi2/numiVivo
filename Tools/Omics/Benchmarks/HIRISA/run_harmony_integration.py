#!/usr/bin/env python3
"""Run pinned Harmony on complete latent coordinates and retain exact output records."""
import argparse
import json
import platform
import resource
import time
from importlib.metadata import version
from pathlib import Path

import harmonypy
import ijson
import numpy as np

from check_integration_response import RECORD, blocks, sha


def run(scores, metadata, protocol_path, seed, out):
    began = time.time()
    assert version('harmonypy') == '2.0.0' and not out.exists()
    protocol = json.loads(protocol_path.read_text())
    assert seed in protocol['seeds']
    n, d = protocol['cells'], protocol['components']
    assert 1 <= n <= 2000000 and 1 <= d <= 64
    bindings = {key: dict(bytes=p.stat().st_size, SHA256=sha(p))
                for key, p in [('scores', scores), ('metadata', metadata), ('protocol', protocol_path)]}
    assert bindings['scores']['SHA256'] == protocol['scoresSHA256']
    assert bindings['metadata']['SHA256'] == protocol['metadataSHA256']
    with metadata.open('rb') as f: samples = list(ijson.items(f, 'samples.item'))
    lookup = {s['id']: s['donorID'] for s in samples}
    assert len(lookup) == len(samples) and all(lookup.values())
    levels = sorted(set(lookup.values())); assert 2 <= len(levels) <= 128
    level_index = {name: i for i, name in enumerate(levels)}
    codes = np.empty(n, dtype=np.int64)
    count = 0
    with metadata.open('rb') as f:
        for count, cell in enumerate(ijson.items(f, 'cells.item'), 1):
            assert count <= n
            codes[count-1] = level_index[lookup[cell['sampleID']]]
    assert count == n and np.all(np.bincount(codes, minlength=len(levels)) > 0)
    # Harmony itself is a resident latent-matrix reference. This layout avoids
    # a separate resident structured-record array; it is not out-of-core Harmony.
    x = np.empty((d, n), dtype=np.float64)
    for start, end, values in blocks(scores, n, d, 8192): x[:, start:end] = values.T
    out.mkdir(parents=True, exist_ok=False)
    start = time.time()
    result = harmonypy.run_harmony(x, {'donor': codes}, ['donor'],
                                  **protocol['options'], random_state=seed, verbose=True)
    fitting_seconds = time.time()-start
    y = result.Z_corr
    assert y.shape == (n, d)
    target = out/'scores.bin'
    with target.open('xb') as f:
        for start in range(0, n, 8192):
            end = min(n, start+8192)
            assert np.isfinite(y[start:end]).all()
            record = np.empty((end-start, d), dtype=RECORD)
            record['row'] = np.arange(start, end, dtype=np.uint32)[:, None]
            record['column'] = np.arange(d, dtype=np.uint32)
            record['value'] = y[start:end]
            f.write(record.tobytes())
    for start, end, values in blocks(target, n, d, 8192):
        assert np.array_equal(values.view(np.uint64), np.ascontiguousarray(y[start:end]).view(np.uint64)), 'Output roundtrip'
    for key, path in [('scores', scores), ('metadata', metadata), ('protocol', protocol_path)]:
        assert path.stat().st_size == bindings[key]['bytes'] and sha(path) == bindings[key]['SHA256']
    package = Path(harmonypy.__file__).parent
    report = dict(status='measured', cells=n, components=d, seed=seed, levels=levels,
                  libraryCellLevelsSHA256=__import__('hashlib').sha256(codes.astype('<i8').tobytes()).hexdigest(),
                  bindings=bindings, output=dict(bytes=target.stat().st_size, SHA256=sha(target)),
                  objectives=list(result.objective_harmony), kmeansRounds=list(result.kmeans_rounds),
                  fittingSeconds=fitting_seconds, seconds=time.time()-began,
                  maximumResidentBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system() == 'Darwin' else 1024),
                  versions=dict(harmonypy=version('harmonypy'), numpy=version('numpy'), ijson=version('ijson')),
                  packageFiles={str(p.relative_to(package)): sha(p) for p in sorted(package.rglob('*')) if p.suffix in ['.py', '.so']},
                  runnerSHA256=sha(Path(__file__)), readerSHA256=sha(Path(__file__).with_name('check_integration_response.py')),
                  platform=platform.platform(), scope=protocol['scope'])
    with (out/'report.json').open('x') as f: json.dump(report, f, indent=2, sort_keys=True, allow_nan=False); f.write('\n')
    print(json.dumps({k: report[k] for k in ['status', 'cells', 'seed', 'fittingSeconds', 'maximumResidentBytes']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['scores', 'metadata', 'protocol', 'out']: p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--seed', type=int, required=True)
    a = p.parse_args(); run(a.scores, a.metadata, a.protocol, a.seed, a.out)
