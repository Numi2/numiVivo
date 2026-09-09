#!/usr/bin/env python3
"""Compare two complete runs, including frozen references and per-cell probabilities."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('first',type=Path);p.add_argument('second',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
checks=[]
for donor in ['human1','human2','human3','human4']:
    for name in ['frozen-reference.npz','knn15-predictions.npz','balancedLogistic-predictions.npz','balancedLogistic-model.npz']:
        with np.load(a.first/donor/name,allow_pickle=False) as x, np.load(a.second/donor/name,allow_pickle=False) as y:
            assert x.files==y.files
            for key in x.files:
                np.testing.assert_array_equal(x[key],y[key],err_msg=donor+'/'+name+'/'+key)
        checks.append(donor+'/'+name)
    for name in ['metrics.json','native/report.json','train/projected.h5ad','query/projected.h5ad']:
        def sha(root):
            with (root/donor/name).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
        assert sha(a.first)==sha(a.second),donor+'/'+name
        checks.append(donor+'/'+name)
assert (a.first/'results.json').read_bytes()==(a.second/'results.json').read_bytes()
a.out.write_text(json.dumps({'passed':True,'checks':checks,'comparison':'Exact arrays, metrics, native reports and projected H5AD bytes'},indent=2)+'\n')
