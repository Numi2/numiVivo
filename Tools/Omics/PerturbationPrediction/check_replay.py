#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('first',type=Path);p.add_argument('second',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=[]
for source in sorted(a.first.rglob('*')):
 if not source.is_file():continue
 relative=source.relative_to(a.first);other=a.second/relative
 if source.suffix=='.npz':
  with np.load(source,allow_pickle=False) as x,np.load(other,allow_pickle=False) as y:
   assert x.files==y.files
   for key in x.files:np.testing.assert_array_equal(x[key],y[key],err_msg=str(relative)+'/'+key)
 else:assert source.read_bytes()==other.read_bytes(),str(relative)
 checks.append(str(relative))
assert sorted(str(x.relative_to(a.first)) for x in a.first.rglob('*') if x.is_file())==sorted(str(x.relative_to(a.second)) for x in a.second.rglob('*') if x.is_file())
a.out.write_text(json.dumps({'passed':True,'exactArraysAndMetadata':checks},indent=2)+'\n')
