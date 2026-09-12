#!/usr/bin/env python3
"""Native ragged projection against AnnData/Awkward, including non-packed roots."""
import argparse
import copy
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import awkward as ak
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=False)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
fingerprint = lambda p: {'bytes': list(bytes.fromhex(sha(p)))}
C, I = ak.contents, ak.index
numbers = C.NumpyArray(np.array([2**53+1, -7, 0, 2**62], dtype=np.int64))
layouts = {
    'numpy': numbers,
    'boolean': C.NumpyArray(np.array([True, False, True, False])),
    'regular': C.RegularArray(C.NumpyArray(np.arange(8, dtype=np.float32)), 2),
    'list': C.ListArray(I.Index64(np.array([2, 0, 3, 1])), I.Index64(np.array([4, 0, 4, 3])), numbers),
    'offsets': C.ListOffsetArray(I.Index64(np.array([0, 2, 2, 3, 4])), numbers),
    'indexed': C.IndexedArray(I.Index32(np.array([3, 1, 3, 0], dtype=np.int32)), numbers),
    'option': C.IndexedOptionArray(I.Index64(np.array([3, -1, 1, 0])), numbers),
    'byte': C.ByteMaskedArray(I.Index8(np.array([1, 0, 1, 1], dtype=np.int8)), numbers, valid_when=True),
    'bits_lsb': C.BitMaskedArray(I.IndexU8(np.array([13], dtype=np.uint8)), numbers, valid_when=True, length=4, lsb_order=True),
    'bits_msb': C.BitMaskedArray(I.IndexU8(np.array([176], dtype=np.uint8)), numbers, valid_when=False, length=4, lsb_order=False),
    'unmasked': C.UnmaskedArray(numbers),
    'union': C.UnionArray(I.Index8(np.array([0, 1, 1, 0], dtype=np.int8)), I.Index64(np.array([0, 0, 1, 1])),
                         [C.NumpyArray(np.array([1.5, 2.5])), ak.Array(['α', 'β']).layout]),
    'record': ak.Array([{'chain': ['A', None], 'score': 1}, {'chain': [], 'score': None},
                        {'chain': ['β'], 'score': -7}, {'chain': ['Z', 'Q'], 'score': 2}]).layout,
}
def write_layout(parent, name, layout):
    form, length, buffers = ak.to_buffers(layout)
    g = parent.create_group(name)
    g.attrs.update({'encoding-type': 'awkward-array', 'encoding-version': '0.1.0',
                    'form': form.to_json(), 'length': length})
    for key, values in buffers.items(): ad.io.write_elem(g, key, values)

source = a.out / 'source.h5ad'
obj = ad.AnnData(X=sparse.csr_matrix(np.arange(12).reshape(4, 3)),
                 obs=pd.DataFrame(index=['c0', 'c1', 'c2', 'c3']), var=pd.DataFrame(index=['g0', 'g1', 'g2']))
obj.raw = obj.copy(); obj.write_h5ad(source)
with h5py.File(source, 'r+') as f:
    for key, layout in layouts.items(): write_layout(f['obsm'], key, layout)
    write_layout(f['varm'], 'ragged', ak.Array([[1, 2], [], [3]]).layout)
    write_layout(f['raw/varm'], 'ragged', ak.Array([['α'], ['b', None], []]).layout)
original = ad.read_h5ad(source)
checks = []
def project(tag, src, rows, cols, success=True, **extras):
    plan = dict(schemaVersion=1, source=fingerprint(src), provenance='Software conformance only',
                observationIndices=rows, featureIndices=cols, **extras)
    pp = a.out / (tag + '.json'); pp.write_text(json.dumps(plan)); dest = a.out / tag
    before = sha(src)
    r = subprocess.run([str(a.binary), 'singlecell-h5ad-project', str(src), '--plan', str(pp), '--output', str(dest)], capture_output=True, text=True)
    (a.out / (tag + '.log')).write_text(r.stdout + r.stderr)
    assert r.returncode == (0 if success else 65), (tag, r.returncode, r.stderr)
    assert sha(src) == before
    if not success:
        assert not dest.exists(); checks.append(tag); return
    got = ad.read_h5ad(dest / 'projected.h5ad'); prev = ad.read_h5ad(src)
    rr = list(range(prev.n_obs)) if rows is None else rows
    cc = list(range(prev.n_vars)) if cols is None else cols
    for slot, indices in [('obsm', rr), ('varm', cc)]:
        for key in getattr(prev, slot):
            expected = getattr(prev, slot)[key][np.array(indices, dtype=np.int64)]
            value = getattr(got, slot)[key]
            assert ak.to_list(value) == ak.to_list(expected), (tag, slot, key)
            assert str(ak.type(value)) == str(ak.type(expected)), (tag, slot, key, ak.type(value), ak.type(expected))
    assert ak.to_list(got.raw.varm['ragged']) == ak.to_list(prev.raw.varm['ragged'])
    np.testing.assert_array_equal(got.X.toarray(), prev.X[rr][:, cc].toarray())
    # Every original nested buffer remains byte-for-byte identical, including
    # integer width and entries that become unreachable after selection.
    with h5py.File(src) as f, h5py.File(dest / 'projected.h5ad') as g:
        for slot in ['obsm', 'varm', 'raw/varm']:
            for key in f[slot]:
                for buffer in f[slot][key]:
                    x, y = f[slot][key][buffer][()], g[slot][key][buffer][()]
                    assert x.dtype == y.dtype and x.tobytes() == y.tobytes()
    replay = subprocess.run([str(a.binary), 'singlecell-h5ad-project-verify', str(dest)], capture_output=True, text=True)
    (a.out / (tag + '-replay.log')).write_text(replay.stdout + replay.stderr)
    assert replay.returncode == 0, replay.stderr
    checks.append(tag); return dest / 'projected.h5ad'

project('identity', source, None, None)
repeat = project('repeated', source, [3, 1, 3, 0], [2, 0, 2])
project('chained', repeat, [2, 0, 2], [1, 0])
project('empty', source, [], [])
zero = a.out / 'zero-source.h5ad'
obj[:0].copy().write_h5ad(zero)
with h5py.File(zero, 'r+') as f:
    write_layout(f['obsm'], 'empty-root', C.EmptyArray())
    write_layout(f['raw/varm'], 'ragged', ak.Array([['α'], ['b', None], []]).layout)
project('empty-root', zero, None, None)
for tag in ['length', 'version', 'missing-index', 'negative-index', 'wrong-mask-type']:
    bad = a.out / (tag + '-source.h5ad'); bad.write_bytes(source.read_bytes())
    with h5py.File(bad, 'r+') as f:
        g = f['obsm/indexed']
        if tag == 'length': g.attrs['length'] = 5
        elif tag == 'version': g.attrs['encoding-version'] = '99.0.0'
        elif tag == 'missing-index': del g['node0-index']
        elif tag == 'negative-index': g['node0-index'][0] = -1
        else:
            g = f['obsm/byte']; del g['node0-mask']; ad.io.write_elem(g, 'node0-mask', np.ones(4, dtype=np.float64))
    project('reject-' + tag, bad, [0, 1], [0], False)
project('reject-work-limit', source, [0, 1], [0], False, maximumElementVisits=10)
project('reject-storage-limit', source, [0, 1], [0], False, maximumOutputBytes=1024)
report = dict(status='passed-software-ragged-interoperability', checks=checks,
              rootForms={k: v.form.to_dict()['class'] for k, v in layouts.items()},
              sourceSHA256=sha(source), binarySHA256=sha(a.binary), checkerSHA256=sha(Path(__file__)),
              anndata=version('anndata'), awkward=version('awkward'))
(a.out / 'report.json').write_text(json.dumps(report, indent=2)); print(json.dumps(report))
