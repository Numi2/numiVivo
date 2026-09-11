#!/usr/bin/env python3
"""Independently verify all partitioned counts, RNA aggregates, QC and identities."""
import argparse,gzip,hashlib,json
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();r=a.root;prepared=r/'prepared-complete'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
ref=np.load(prepared/'reference.npz');inv=json.loads((prepared/'inventory.json').read_text());semantic=hashlib.sha256();visited=0
barcodes=gzip.open(r/'GSM4367979_exp1-5.barcodes.tsv.gz','rt').read().splitlines();features=[x.split('\t') for x in gzip.open(r/'GSM4367979_exp1-5.features.tsv.gz','rt').read().splitlines()]
with h5py.File(prepared/'original.h5ad') as src,h5py.File(r/'rna-projection/projected.h5ad') as rn,h5py.File(r/'guide-projection/projected.h5ad') as gd:
 assert src['obs/_index'].asstr()[:].tolist()==barcodes
 assert src['var/_index'].asstr()[:].tolist()==[x[0] for x in features]
 assert src['var/feature_type'].asstr()[:].tolist()==[x[2] for x in features]
 masks=[ref['rnaMask'],~ref['rnaMask']];outputs=[rn,gd]
 for out,mask in zip(outputs,masks):
  assert out['var/_index'].asstr()[:].tolist()==[x[0] for x,m in zip(features,mask) if m]
  for column in src['obs']:np.testing.assert_array_equal(src['obs/'+column].asstr()[:],out['obs/'+column].asstr()[:])
  assert tuple(out['X'].attrs['shape'])==(len(barcodes),int(mask.sum()))
 sp=src['X/indptr'][:];ops=[o['X/indptr'][:] for o in outputs]
 for start in range(0,len(barcodes),256):
  stop=min(start+256,len(barcodes));lo,hi=map(int,sp[[start,stop]])
  indices=src['X/indices'][lo:hi];values=src['X/data'][lo:hi];rows=np.repeat(np.arange(start,stop)+1,np.diff(sp[start:stop+1]))
  triples=np.column_stack((indices+1,rows,values)).astype('<i8');semantic.update(triples.tobytes());visited+=len(values)
  for out,mask,op in zip(outputs,masks,ops):
   keep=mask[indices];expected_indices=(np.cumsum(mask)-1)[indices[keep]];expected_values=values[keep]
   ol,oh=map(int,op[[start,stop]]);actual_indices=out['X/indices'][ol:oh];actual_values=out['X/data'][ol:oh]
   expected_rows=rows[keep];actual_rows=np.repeat(np.arange(start,stop)+1,np.diff(op[start:stop+1]))
   # Canonical coordinates permit valid sparse storage order differences.
   e=np.lexsort((expected_indices,expected_rows));v=np.lexsort((actual_indices,actual_rows))
   np.testing.assert_array_equal(expected_rows[e],actual_rows[v]);np.testing.assert_array_equal(expected_indices[e],actual_indices[v]);np.testing.assert_array_equal(expected_values[e],actual_values[v])
 assert visited==inv['entries'] and semantic.hexdigest()==inv['semanticTriplesSHA256']
report=json.loads((r/'aggregate/report.json').read_text());bulk=report['pseudobulk'];m=bulk['matrix'];groups=bulk['groups'];look={name:i for i,name in enumerate(inv['groups'])}
assert report['canonicalNonzeros']==inv['rnaEntries'] and len(report['metadata']['cells'])==inv['cells']
assert bulk['featureIDs']==[x[0] for x,keep in zip(features,ref['rnaMask']) if keep]
assert [x['barcode'] for x in report['metadata']['cells']]==barcodes
for key,name in [('totalCounts','rnaTotals'),('detectedFeatures','rnaDetected'),('mitochondrialCounts','mitochondrial')]:np.testing.assert_array_equal([x[key] for x in report['quality']],ref[name])
for i,g in enumerate(groups):
 key=g['cellGroup']+'|'+g['condition'];expected_row=look[key];lo,hi=m['rowOffsets'][i:i+2]
 observed=np.zeros(int(ref['rnaMask'].sum()),dtype=np.int64);observed[np.asarray(m['featureIndices'][lo:hi],dtype=np.int64)]=m['counts'][lo:hi]
 np.testing.assert_array_equal(observed,ref['counts'][expected_row,ref['rnaMask']]);np.testing.assert_array_equal(g['sourceCellIndices'],np.flatnonzero(ref['group']==expected_row))
 assert g['biologicalReplicateID']=='K562-pooled-replication-unresolved' and 'donorID' not in g
assert len(groups)==len(look)
result=dict(status='passed',cells=inv['cells'],rnaFeatures=inv['rnaFeatures'],guideFeatures=inv['guideFeatures'],sourceEntries=visited,rnaEntries=inv['rnaEntries'],groups=len(groups),allPartitionCoordinatesAndCountsExact=True,sourceSemanticSHA256=semantic.hexdigest(),aggregateCountsAndMembershipExact=True,cellQCExact=True,reportSHA256=sha(r/'aggregate/report.json'),predictionFitted=False)
(r/'independent-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
