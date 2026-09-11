from pathlib import Path
import sys,json,gzip,hashlib,time,collections
r=Path(__file__).parent;sys.path.insert(0,str(r/'deps'))
import indexed_gzip,h5py,numpy as np,pandas as pd,anndata as ad
from scipy import sparse
out=r/'handoff';out.mkdir(exist_ok=False)
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
checks=json.loads((r/'source-count-validation.json').read_text());assert checks['status']=='passed-count-domain-and-CSR'
o=json.loads(gzip.decompress((r/'obs-metadata.json.gz').read_bytes()));v=json.loads(gzip.decompress((r/'var-metadata.json.gz').read_bytes()))
selected=np.array(o['ct3'])=='B';features=np.array(v['genome'])=='GRCh38';si=np.flatnonzero(selected);fi=np.flatnonzero(features);assert len(si)==15272 and len(fi)==20303
blocks=[];start=time.time()
with indexed_gzip.IndexedGzipFile(str(r/'GSE181897_concat.4.raw.h5ad.gz'),index_file=str(r/'source.gzidx'),readbuf_size=1024*1024) as gz,h5py.File(gz,'r',rdcc_nbytes=64*1024*1024) as f:
 ptr=f['X/indptr'][:].astype(np.int64)
 for begin in range(0,len(selected),1024):
  end=min(len(selected),begin+1024)
  if not selected[begin:end].any():continue
  lo,hi=map(int,(ptr[begin],ptr[end]));values=f['X/data'][lo:hi];cols=f['X/indices'][lo:hi]
  assert np.all(np.isfinite(values)&(values>=0)&(values==np.floor(values))) and values.max(initial=0)<=2**24
  x=sparse.csr_matrix((values.astype(np.uint32),cols,ptr[begin:end+1]-lo),shape=(end-begin,len(features)))
  blocks.append(x[selected[begin:end]][:,features])
counts=sparse.vstack(blocks,format='csr');counts.sort_indices();assert counts.shape==(len(si),len(fi))
with np.load(r/'source-cell-library-counts.npz') as z:np.testing.assert_array_equal(np.asarray(counts.sum(axis=1,dtype=np.uint64)).ravel(),z['RNA'][si])
obs=pd.DataFrame({k:np.array(val)[si] for k,val in o.items() if k!='_index'},index=pd.Index(np.array(o['_index'])[si],dtype=object));var=pd.DataFrame({k:np.array(val)[fi] for k,val in v.items() if k!='_index'},index=pd.Index(np.array(v['_index'])[fi],dtype=object))
obs['native_sample']=['|'.join((d,c,b)) for d,c,b in zip(obs['exp_id'],obs['cond'],obs['batch'])]
samples=[]
for sample in sorted(set(obs['native_sample'])):
 d,c,b=sample.split('|');donor='GSE181897:exp_id:'+d;samples.append(dict(id=sample,donorID=donor,biologicalReplicateID=donor,condition='source-code:'+c,batchID='GSE181897:pool:'+b,organism='NCBITaxon:9606'))
obj=ad.AnnData(counts,obs=obs,var=var);obj.write_h5ad(out/'B-source-codes-RNA.h5ad',compression='gzip')
mapping=dict(schemaVersion=1,id='GSE181897-B-all-source-condition-codes',evidence='measured',sourceDescription='All released ct3=B cells and genome=GRCh38 features from original GEO raw H5AD. External indexed-gzip extraction, preserving original condition codes without biological translation. Includes B_Naive, B_Mem and PB author ct2 labels. Feature_types uniformly Gene Expression is contradicted by genome BD99AbSeq for antibody features; those are excluded from RNA denominator.',countUnit='umiCount',matrixPath='X',sampleColumn='native_sample',groupColumn='ct3',featureNameColumn='_index',samples=samples)
# Feature index is already the primary source symbol; no separate name column needed.
mapping.pop('featureNameColumn')
write(out/'aggregate-plan.json',dict(schemaVersion=1,mapping=mapping,contrasts=[]))
keys=list(zip(obs['exp_id'],obs['cond']));unique=sorted(set(keys));lookup={k:i for i,k in enumerate(unique)};assignment=np.array([lookup[k] for k in keys]);incidence=sparse.csr_matrix((np.ones(len(keys),dtype=np.uint64),(assignment,np.arange(len(keys)))),shape=(len(unique),len(keys)));bulk=(incidence@counts).tocsr();sparse.save_npz(out/'reference-pseudobulk.npz',bulk)
np.savez_compressed(out/'source-selection.npz',cellIndices=si,featureIndices=fi,groupAssignment=assignment)
write(out/'reference-groups.json',[dict(expID=k[0],conditionCode=k[1],cells=int(np.count_nonzero(assignment==i)),libraryCounts=int(bulk[i].sum()),sourceCellIndices=np.flatnonzero(assignment==i).tolist()) for i,k in enumerate(unique)])
prior=json.loads((Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs/panel.json')).read_text());symbols=set(var.index);panel=[x for x in prior if x.removeprefix('symbol|') in symbols];write(out/'candidate-panel.json',panel);write(out/'unsupported-prior-panel.json',[x for x in prior if x not in panel])
# Counts identify available strata, not which code denotes IFN-beta or control.
crosstab=[]
for donor in sorted(set(obs['exp_id']),key=int):crosstab.append(dict(expID=donor,conditionCellCounts=dict(collections.Counter(obs.loc[obs['exp_id']==donor,'cond']))))
write(out/'source-donor-condition-counts.json',crosstab)
write(out/'preparation.json',dict(status='complete-source-code-handoff',sourceSHA256=sha(r/'GSE181897_concat.4.raw.h5ad.gz'),sourceCountValidationSHA256=sha(r/'source-count-validation.json'),protocolSHA256=sha(r/'PROTOCOL.md'),scriptSHA256=sha(Path(__file__)),cells=counts.shape[0],features=counts.shape[1],countRecords=counts.nnz,groups=len(unique),donors=len(set(obs['exp_id'])),sampleDefinitions=len(samples),RNAUMIs=int(counts.sum(dtype=np.uint64)),priorPanel=len(prior),candidatePanel=len(panel),unsupportedPanel=len(prior)-len(panel),conditionMeaningQualified=False,fitStarted=False,elapsedSeconds=time.time()-start))
write(out/'input-freeze.json',dict(schemaVersion=1,files={str(p.relative_to(out)):sha(p) for p in sorted(out.iterdir()) if p.is_file()},fitStarted=False,conditionMeaningQualified=False))
print((out/'preparation.json').read_text())
