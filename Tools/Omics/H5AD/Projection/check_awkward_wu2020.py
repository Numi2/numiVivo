#!/usr/bin/env python3
"""Check the complete public Wu 2020 3k AIRR modality through native projection."""
import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import anndata as ad
import awkward as ak
import h5py
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--out', type=Path, required=True)
a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=False)
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
assert hashlib.md5(a.source.read_bytes()).hexdigest() == '12c57c790f8a403751304c9de5a18cbf'
source = a.out / 'airr.h5ad'
# Mechanical extraction of the complete embedded AnnData group, not a count
# conversion or claim that native H5MU import supports AIRR records.
with h5py.File(a.source) as f, h5py.File(source, 'x') as g:
    original = f['mod/airr']
    for key, value in original.attrs.items(): g.attrs[key] = value
    for key in original: original.copy(key, g, name=key)
obj = ad.read_h5ad(source); assert obj.shape == (3000, 0)
records = ak.to_list(obj.obsm['airr'])
cases = []
for tag, rows in [('identity', list(range(3000))), ('reverse-and-repeat', list(range(2999, -1, -1)) + [0, 2999])]:
    plan = dict(schemaVersion=1, source={'bytes': list(bytes.fromhex(sha(source)))},
                provenance='Complete public Scirpy Wu 2020 3k AIRR modality; interoperability only', observationIndices=rows)
    pp = a.out / (tag + '.json'); pp.write_text(json.dumps(plan)); dest = a.out / tag
    result = subprocess.run([str(a.binary), 'singlecell-h5ad-project', str(source), '--plan', str(pp), '--output', str(dest)], capture_output=True, text=True)
    (a.out / (tag + '.log')).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    got = ad.read_h5ad(dest / 'projected.h5ad')
    assert ak.to_list(got.obsm['airr']) == [records[i] for i in rows]
    assert str(ak.type(got.obsm['airr'])).split(' * ', 1)[1] == str(ak.type(obj.obsm['airr'])).split(' * ', 1)[1]
    pd.testing.assert_frame_equal(got.obs, obj.obs.iloc[rows])
    pd.testing.assert_frame_equal(got.var, obj.var)
    np.testing.assert_array_equal(got.X, obj.X[rows])
    with h5py.File(source) as f, h5py.File(dest / 'projected.h5ad') as g:
        for key in f['obsm/airr']:
            x, y = f['obsm/airr'][key][()], g['obsm/airr'][key][()]
            assert x.dtype == y.dtype and x.tobytes() == y.tobytes()
    replay = subprocess.run([str(a.binary), 'singlecell-h5ad-project-verify', str(dest)], capture_output=True, text=True)
    (a.out / (tag + '-replay.log')).write_text(replay.stdout + replay.stderr)
    assert replay.returncode == 0, replay.stderr
    cases.append(dict(case=tag, outputCells=len(rows), outputSHA256=sha(dest / 'projected.h5ad'), nativeReplay=True, allAIRRRecordsExact=True))
report = dict(status='passed-real-receptor-interoperability-not-function-prediction', cells=3000,
              chains=sum(map(len, records)), fields=ak.fields(obj.obsm['airr']), cases=cases,
              sourceURL='https://exampledata.scverse.org/scirpy/wu2020_3k.h5mu',
              upstreamSHA256=sha(a.source), extractedH5ADSHA256=sha(source),
              binarySHA256=sha(a.binary), checkerSHA256=sha(Path(__file__)),
              anndata=version('anndata'), awkward=version('awkward'))
(a.out / 'report.json').write_text(json.dumps(report, indent=2)); print(json.dumps(report))
