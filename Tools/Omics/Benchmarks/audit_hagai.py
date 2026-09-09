#!/usr/bin/env python3
"""Record why the pinned pertpy Hagai release is ineligible for raw-count DE."""
import argparse
import hashlib
import json
from pathlib import Path
import h5py
import numpy as np
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
sha='41b98069fe5d88193cec2b1dcb0eb43843e4de83c8aa22ea29cef1bb65ba064e'
with a.source.open('rb') as stream:
    assert hashlib.file_digest(stream,'sha256').hexdigest()==sha
with h5py.File(a.source) as f:
    assert 'raw' not in f and len(f['layers'])==0
    values=f['X/data'][:10000]
    fractional=values[values!=np.floor(values)]
    assert len(fractional)>0
    report=dict(schemaVersion=1,status='ineligible-for-raw-count-DE',sourceSHA256=sha,
        sourceURL='https://ndownloader.figshare.com/files/46978846',publication='https://doi.org/10.1038/s41586-018-0657-2',
        reason='X contains noninteger expression and no raw or counts layer is supplied',
        inspectedStoredValues=len(values),fractionalValuesInInspectedPrefix=len(fractional),
        exampleFractionalValues=fractional[:5].astype(float).tolist(),
        policy='Do not round expression, invert normalization or synthesize counts. Acquire original raw data before count-model benchmarking.')
with a.out.open('x') as stream:json.dump(report,stream,indent=2);stream.write('\n')
print(json.dumps(report,indent=2))
