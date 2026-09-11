#!/usr/bin/env python3
"""Stream original GEO integer MEX into H5AD; retain every source cell and feature."""
from decimal import Decimal
import argparse,csv,gzip,hashlib,itertools,json,time
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('root',type=Path);a=p.parse_args();r=a.root;out=r/'prepared-complete';out.mkdir(exist_ok=False)
def dump(name,value):(out/name).write_text(json.dumps(value,sort_keys=True,separators=(',',':'))+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def strings(group,name,values):
 d=group.create_dataset(name,data=np.asarray(values,dtype=object),dtype=h5py.string_dtype());d.attrs['encoding-type']='string-array';d.attrs['encoding-version']='0.2.0'
barcodes=gzip.open(r/'GSM4367979_exp1-5.barcodes.tsv.gz','rt').read().splitlines();features=[x.split('\t') for x in gzip.open(r/'GSM4367979_exp1-5.features.tsv.gz','rt').read().splitlines()]
assert len(barcodes)==len(set(barcodes))==40997 and len(features)==33758
ids=[x[0] for x in features];names=[x[1] for x in features];assert len(set(ids))==len(ids) and set(x[2] for x in features)=={'Gene Expression','CRISPR Guide Capture'}
assignments=list(csv.DictReader(gzip.open(r/'GSM4367979_exp1-5.cell_identities.csv.gz','rt')))
assert len({x['cell_barcode'] for x in assignments})==len(assignments)
by_barcode={x['cell_barcode']:x for x in assignments}
barcode_set=set(barcodes)
extra_assignments=[x for x in assignments if x['cell_barcode'] not in barcode_set]
dump('assignments-outside-matrix.json',extra_assignments)
platforms=[];guides=[];selected=[];reasons=[]
for i,b in enumerate(barcodes):
 gem=b.rsplit('-',1)[1];assert gem in ['1','2','3','4','5'];platforms.append('gemgroup-'+gem)
 x=by_barcode.get(b);guide=x['guide_identity'] if x else 'numivivo-unassigned';guides.append(guide)
 if x:
  multiplicity=Decimal(x['number_of_cells']);assert multiplicity.is_finite() and multiplicity>=0 and multiplicity==int(multiplicity)
  assert x['gemgroup']==gem and x['good_coverage'] in ['True','False']
 reason='missing-assignment' if x is None else ('low-guide-coverage' if x['good_coverage']!='True' else ('multiple-guide-calls' if Decimal(x['number_of_cells'])!=1 else 'included'))
 reasons.append(reason)
 if reason=='included':selected.append(i)
labels=[g+'|'+guide for g,guide in zip(platforms,guides)];unique=sorted(set(labels));index={v:i for i,v in enumerate(unique)}
group=np.asarray([index[x] for x in labels],dtype=np.int64)
reference=np.zeros((len(unique),len(features)),dtype=np.int64);totals=np.zeros(len(barcodes),dtype=np.int64);detected=np.zeros(len(barcodes),dtype=np.int64)
mito=np.asarray([name.startswith('MT-') for name in names]);mitototals=np.zeros(len(barcodes),dtype=np.int64)
rna=np.asarray([x[2]=='Gene Expression' for x in features]);assert np.count_nonzero(rna)==33694
rna_totals=np.zeros(len(barcodes),dtype=np.int64);rna_detected=np.zeros(len(barcodes),dtype=np.int64)
source=r/'GSM4367979_exp1-5.matrix.mtx.gz';source_record=json.loads((r/'matrix-retrieval.json').read_text());assert sha(source)==source_record['SHA256']
start=time.monotonic();semantic=hashlib.sha256();entry_count=0;last_cell=-1;last_features=np.array([],dtype=np.int64);maximum=0
with h5py.File(out/'original.h5ad','x') as h,gzip.open(source,'rt') as stream:
 assert next(stream).strip()=='%%MatrixMarket matrix coordinate integer general'
 line=next(stream)
 while line.startswith('%'):line=next(stream)
 ng,nc,nnz=map(int,line.split());assert (ng,nc,nnz)==(33758,40997,129839577)
 h.attrs['encoding-type']='anndata';h.attrs['encoding-version']='0.1.0'
 for axis,axis_ids,columns in [('obs',barcodes,{'sample':labels,'gemgroup':platforms,'guide_identity':guides,'selection_reason':reasons}),('var',ids,{'gene_symbol':names,'feature_type':[x[2] for x in features]})]:
  frame=h.create_group(axis);frame.attrs['encoding-type']='dataframe';frame.attrs['encoding-version']='0.2.0';frame.attrs['_index']='_index';frame.attrs['column-order']=np.asarray(list(columns),dtype=h5py.string_dtype());strings(frame,'_index',axis_ids)
  for name,values in columns.items():strings(frame,name,values)
 x=h.create_group('X');x.attrs['encoding-type']='csr_matrix';x.attrs['encoding-version']='0.1.0';x.attrs['shape']=np.array([nc,ng],dtype=np.int64)
 data=x.create_dataset('data',shape=(nnz,),dtype='<u4',chunks=(65536,),compression='gzip',compression_opts=1)
 indices=x.create_dataset('indices',shape=(nnz,),dtype='<i4',chunks=(65536,),compression='gzip',compression_opts=1)
 while True:
  lines=list(itertools.islice(stream,262144))
  if not lines:break
  block=np.fromstring(''.join(lines),sep=' ',dtype=np.int64).reshape(-1,3);assert len(block)==len(lines)
  feature=block[:,0]-1;cell=block[:,1]-1;values=block[:,2]
  assert np.all((feature>=0)&(feature<ng)) and np.all((cell>=0)&(cell<nc)) and np.all((values>0)&(values<=2**32-1))
  assert cell[0]>=last_cell and np.all(np.diff(cell)>=0)
  assert len(np.unique(cell*ng+feature))==len(block), 'duplicate source coordinates within block'
  if cell[0]==last_cell:assert len(np.intersect1d(last_features,feature[cell==last_cell]))==0, 'duplicate source coordinates across blocks'
  tail=feature[cell==cell[-1]]
  last_features=np.concatenate((last_features,tail)) if cell[-1]==last_cell else tail.copy()
  last_cell=int(cell[-1]);maximum=max(maximum,int(values.max()))
  end=entry_count+len(block);assert end<=nnz
  data[entry_count:end]=values;indices[entry_count:end]=feature
  np.add.at(reference,(group[cell],feature),values);np.add.at(totals,cell,values);np.add.at(detected,cell,1);np.add.at(mitototals,cell,values*mito[feature])
  np.add.at(rna_totals,cell,values*rna[feature]);np.add.at(rna_detected,cell,rna[feature].astype(np.int64))
  semantic.update(np.asarray(block,dtype='<i8').tobytes());entry_count=end
 assert entry_count==nnz
 pointers=np.concatenate(([0],np.cumsum(detected)));x.create_dataset('indptr',data=pointers)
 for name in ['uns','obsm','varm','obsp','varp','layers']:
  q=h.create_group(name);q.attrs['encoding-type']='dict';q.attrs['encoding-version']='0.1.0'
np.savez_compressed(out/'reference.npz',counts=reference,totals=totals,detected=detected,mitochondrial=mitototals,group=group,selected=np.asarray(selected,dtype=np.int64),rnaMask=rna,rnaTotals=rna_totals,rnaDetected=rna_detected)
samples=[dict(id=label,condition=label.split('|',1)[1],batchID=label.split('|',1)[0],biologicalReplicateID='K562-pooled-replication-unresolved',organism='NCBITaxon:9606') for label in unique]
mapping=dict(schemaVersion=1,id='replogle2020-upr-five-gemgroups',evidence='measured',sourceDescription='Complete GEO GSM4367979 exp1-5 integer matrix; technical gemgroups are not biological replicates; original guide labels retained',countUnit='umiCount',matrixPath='X',samples=samples,sampleColumn='sample',groupColumn='gemgroup',featureNameColumn='gene_symbol',mitochondrialFeatureIDs=[i for i,m in zip(ids,mito) if m])
dump('plan.json',dict(schemaVersion=1,mapping=mapping,contrasts=[]))
for assay,mask in [('rna',rna),('guide',~rna)]:
 dump(assay+'-projection-plan.json',dict(schemaVersion=1,source=dict(bytes=list(bytes.fromhex(sha(out/'original.h5ad')))),provenance='All original cells; exact original feature-type partition: '+assay,featureIndices=np.flatnonzero(mask).tolist(),maximumElementVisits=500000000,**({'maximumOutputBytes':2147483648} if assay=='rna' else {})))
dump('inventory.json',dict(cells=nc,features=ng,entries=nnz,rnaFeatures=int(rna.sum()),guideFeatures=int((~rna).sum()),rnaEntries=int(rna_detected.sum()),rnaUMI=int(rna_totals.sum()),maximumCount=maximum,totalUMI=int(totals.sum()),selectedCells=len(selected),groups=unique,reasonCounts={x:reasons.count(x) for x in sorted(set(reasons))},guideLabels=sorted(set(guides)),semanticTriplesSHA256=semantic.hexdigest(),sourceSHA256=source_record['SHA256'],h5adSHA256=sha(out/'original.h5ad'),seconds=time.monotonic()-start,allOriginalCellsAndFeaturesRetained=True,countsScoredForPrediction=False))
print(json.dumps(json.loads((out/'inventory.json').read_text()) | {'groups':len(unique)}))
