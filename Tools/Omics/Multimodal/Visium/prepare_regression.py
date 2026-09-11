#!/usr/bin/env python3
"""Create a declared two-spot structural input; not a biological benchmark."""
import argparse,json
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--plan',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
(a.out/'outs/spatial').mkdir(parents=True,exist_ok=False)
plan=json.loads(a.plan.read_text());plan['counts']['id']='visium-structural';plan['counts']['evidence']='synthetic';plan['counts']['sourceDescription']='Two-spot structural control; no biological evidence'
(a.out/'plan.json').write_text(json.dumps(plan)+'\n')
(a.out/'outs/spatial/tissue_positions_list.csv').write_text('b1,1,0,2,11,21\nb0,1,0,0,10,20\noutside,0,2,0,50,60\n')
with h5py.File(a.out/'outs/filtered_feature_bc_matrix.h5','w') as f:
 m=f.create_group('matrix');m['shape']=np.array([2,2],dtype='int64');m['indptr']=np.array([0,1,3],dtype='int64');m['indices']=np.array([0,0,1],dtype='int64');m['data']=np.array([2,3,5],dtype='int64');m['barcodes']=np.array([b'b0',b'b1'])
 g=m.create_group('features');g['id']=np.array([b'gene0',b'gene1']);g['name']=np.array([b'G0',b'G1']);g['feature_type']=np.array([b'Gene Expression',b'Gene Expression']);g['genome']=np.array([b'GRCh38',b'GRCh38'])
