#!/usr/bin/env python3
"""Execute and independently verify the complete confident five-gemgroup cohort."""
import argparse,gzip,json,os,subprocess,time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();r=a.root;binary=r/'runtime/numivivo-omics';rows=[]
for label,args in [('selected-aggregate',['singlecell-h5ad-pseudobulk',r/'rna-projection/projected.h5ad','--plan',r/'selected-plan.json','--output',r/'selected-aggregate']),('selected-verify',['singlecell-h5ad-pseudobulk-verify',r/'selected-aggregate'])]:
 assert not (r/(label+'.log')).exists()
 with (r/(label+'.log')).open('wb') as err,(r/(label+'.stdout')).open('wb') as out:
  start=time.monotonic();q=subprocess.run(['/usr/bin/time','-l',str(binary),*map(str,args)],stdout=out,stderr=err,env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY='/Users/n/numivivo-multiassay-hdf5-20260909/libhdf5.dylib'))
 assert q.returncode==0,(label,(r/(label+'.log')).read_text());rows.append(dict(label=label,returnCode=q.returncode,seconds=time.monotonic()-start))
ref=np.load(r/'selected-reference.npz');selected=ref['selected'];inv=json.loads((r/'prepared-complete/inventory.json').read_text());look={x:i for i,x in enumerate(inv['groups'])};report=json.loads((r/'selected-aggregate/report.json').read_text());bulk=report['pseudobulk'];m=bulk['matrix']
np.testing.assert_array_equal(report['sourceObservationIndices'],selected);assert report['sourceCellCount']==40997 and len(report['metadata']['cells'])==32829
barcodes=gzip.open(r/'GSM4367979_exp1-5.barcodes.tsv.gz','rt').read().splitlines();assert [x['barcode'] for x in report['metadata']['cells']]==[barcodes[i] for i in selected]
for field,key in [('totalCounts','totals'),('detectedFeatures','detected'),('mitochondrialCounts','mitochondrial')]:np.testing.assert_array_equal([x[field] for x in report['quality']],ref[key][selected])
assert report['canonicalNonzeros']==int(ref['detected'][selected].sum())
for i,g in enumerate(bulk['groups']):
 row=look[g['cellGroup']+'|'+g['condition']];lo,hi=m['rowOffsets'][i:i+2];actual=np.zeros(33694,dtype=np.int64);actual[np.asarray(m['featureIndices'][lo:hi],dtype=np.int64)]=m['counts'][lo:hi]
 np.testing.assert_array_equal(actual,ref['counts'][row]);np.testing.assert_array_equal(g['sourceCellIndices'],np.flatnonzero(ref['group'][selected]==row))
 assert g['biologicalReplicateID']=='K562-pooled-replication-unresolved' and 'donorID' not in g
assert len(bulk['groups'])==160
(r/'selected-checks.json').write_text(json.dumps(dict(status='passed',runs=rows,sourceCells=40997,selectedCells=32829,features=33694,groups=160,canonicalNonzeros=report['canonicalNonzeros'],allCountsQualityAndMembershipExact=True,nativeReplayPassed=True,predictionFitted=False),indent=2)+'\n');print('Complete confident cohort: native counts, QC, membership and replay exact')
