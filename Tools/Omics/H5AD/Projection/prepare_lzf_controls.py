#!/usr/bin/env python3
"""Create independent compressed chunks and retain the actual error owner for ASan."""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--source-root', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=False)
owner = a.source_root / 'Sources/NumiVivoKit/Omics/VivoSparseCounts.swift'
source = owner.read_text()
(a.out / 'ErrorOwner.swift').write_text(source[:source.index('struct VivoOmicsJSONKey')])
rows = []
with h5py.File(a.out / 'lzf-chunks.h5', 'x') as f:
    for i in range(24):
        x = np.tile(np.random.default_rng(i).integers(0, 256, size=32+i*3, dtype=np.uint8), 256)
        d = f.create_dataset(str(i), data=x, chunks=x.shape, compression='lzf'); f.flush()
        mask, compressed = d.id.read_direct_chunk((0,)); assert mask == 0
        name = f'lzf-chunk-{i}'
        (a.out / (name + '.lzf')).write_bytes(compressed)
        (a.out / (name + '.raw')).write_bytes(x.tobytes())
        rows.append({'compressed': name + '.lzf', 'expected': name + '.raw'})
(a.out / 'chunks.json').write_text(json.dumps(rows))
(a.out / 'source.json').write_text(json.dumps({'errorOwnerSHA256': hashlib.sha256(owner.read_bytes()).hexdigest(),
                                               'h5py': h5py.__version__, 'numpy': np.__version__}, indent=2))
