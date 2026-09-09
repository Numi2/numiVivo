#!/usr/bin/env python3
"""Reconstruct retained Scanpy reference outputs through direct Scanorama assembly."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import scanorama
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert not a.out.exists()
d=json.loads((a.root/'checks.json').read_text());assert d['status']=='full-cohort-measured'
path=a.inputs/('inputs.npz' if d['cohort']=='ding' else 'evaluation-inputs.npz');assert hashlib.sha256(path.read_bytes()).hexdigest()==d['inputSHA256'];inp=np.load(path,allow_pickle=False);x=inp['scores'];batch=inp['methods' if d['cohort']=='ding' else 'donors'];scale=d.get('scaleValue',1.0);levels=np.unique(batch);rows=[np.flatnonzero(batch==b) for b in levels];datasets=[x[ix].copy()/scale for ix in rows]
# Recompute matching, rather than trusting the retained anchor file alone.
alignments,matches=scanorama.find_alignments(datasets,knn=20,approx=False,alpha=.1,verbose=0)
old=json.loads((a.root/'alignment-order.json').read_text());assert alignments==[(v['source'],v['reference']) for v in old]
anchors=np.load(a.root/'anchors.npz');assert set(anchors.files)=={f'{i}-{j}' for i,j in matches}
for (i,j),pairs in matches.items():np.testing.assert_array_equal(anchors[f'{i}-{j}'],np.array([(rows[i][s],rows[j][t]) for s,t in sorted(pairs)],dtype=np.int64).reshape(-1,2))
corrected=scanorama.assemble(datasets,knn=20,sigma=15,approx=False,alpha=.1,ds_names=list(levels),batch_size=256,verbose=0,alignments=alignments,matches=matches)
y=np.empty_like(x)
for ix,v in zip(rows,corrected):y[ix]=v*scale
old=np.load(a.root/'scores.npz');np.testing.assert_array_equal(old['sourceRow'],np.arange(len(x)));np.testing.assert_array_equal(old['groupedSourceRow'],np.concatenate(rows));np.testing.assert_array_equal(old['scores'],y)
a.out.write_text(json.dumps(dict(status='passed',cohort=d['cohort'],scaleValue=scale,allAnchorsAndAlignmentOrderExactlyReconstructed=True,allSourceRowsExactlyRestored=True,allOutputCoordinatesExactlyReconstructed=True,maximumCoordinateError=float(np.max(np.abs(old['scores']-y))),qualification='Independent invocation of the same pinned reference library through direct assembly agrees bit-for-bit with the public Scanpy wrapper. This is replay/adapter verification, not an independent algorithm implementation or biological acceptance.'),indent=2)+'\n')
