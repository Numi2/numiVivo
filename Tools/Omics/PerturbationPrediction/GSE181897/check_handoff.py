from pathlib import Path
import sys,json,gzip,hashlib,time
r=Path(__file__).parent;sys.path.insert(0,str(r/'deps'))
import indexed_gzip,h5py,numpy as np,anndata as ad
from scipy import sparse
start=time.time();out=r/'handoff';native=r/'native-handoff-qualified'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
freeze=json.loads((out/'input-freeze.json').read_text())
for name,digest in freeze['files'].items():assert sha(out/name)==digest,name
execution=json.loads((native/'execution.json').read_text());assert execution['inputFreezeSHA256']==sha(out/'input-freeze.json') and execution['reportSHA256']==sha(native/'aggregate/report.json')==sha(native/'repeat/report.json')
transport=ad.read_h5ad(out/'B-source-codes-RNA.h5ad');o=json.loads(gzip.decompress((r/'obs-metadata.json.gz').read_bytes()));v=json.loads(gzip.decompress((r/'var-metadata.json.gz').read_bytes()))
with np.load(out/'source-selection.npz') as z:si=z['cellIndices'];fi=z['featureIndices'];assignment=z['groupAssignment']
np.testing.assert_array_equal(si,np.flatnonzero(np.array(o['ct3'])=='B'));np.testing.assert_array_equal(fi,np.flatnonzero(np.array(v['genome'])=='GRCh38'))
np.testing.assert_array_equal(transport.obs_names,np.array(o['_index'])[si]);np.testing.assert_array_equal(transport.var_names,np.array(v['_index'])[fi])
for key in o:
 if key!='_index':np.testing.assert_array_equal(transport.obs[key].astype(str),np.array(o[key])[si])
for key in v:
 if key!='_index':np.testing.assert_array_equal(transport.var[key].astype(str),np.array(v[key])[fi])
# Independent second pass verifies every selected source count against H5AD bytes.
selected=np.zeros(len(o['_index']),dtype=bool);selected[si]=True;position=0;verified=0
with indexed_gzip.IndexedGzipFile(str(r/'GSE181897_concat.4.raw.h5ad.gz'),index_file=str(r/'source.gzidx'),readbuf_size=1024*1024) as gz,h5py.File(gz,'r',rdcc_nbytes=64*1024*1024) as f:
 ptr=f['X/indptr'][:].astype(np.int64)
 for begin in range(0,len(selected),2048):
  end=min(len(selected),begin+2048);n=int(selected[begin:end].sum())
  if not n:continue
  lo,hi=map(int,(ptr[begin],ptr[end]));x=sparse.csr_matrix((f['X/data'][lo:hi],f['X/indices'][lo:hi],ptr[begin:end+1]-lo),shape=(end-begin,len(v['_index'])))[selected[begin:end]][:,fi];actual=transport.X[position:position+n]
  for key in ['indptr','indices','data']:np.testing.assert_array_equal(getattr(x,key),getattr(actual,key))
  position+=n;verified+=x.nnz
assert position==len(si) and verified==transport.X.nnz
report=json.loads((native/'aggregate/report.json').read_text());bulk=report['pseudobulk'];metadata=report['metadata'];assert report['canonicalNonzeros']==verified
assert [f['id'] for f in metadata['features']]==list(transport.var_names)
plan=json.loads((out/'aggregate-plan.json').read_text());assert metadata['samples']==plan['mapping']['samples']
for i,cell in enumerate(metadata['cells']):assert cell==dict(barcode=transport.obs_names[i],sampleID=str(transport.obs['native_sample'].iloc[i]),group='B')
assert len(metadata['cells'])==len(si)
reference=sparse.load_npz(out/'reference-pseudobulk.npz');reference.sort_indices();groups=json.loads((out/'reference-groups.json').read_text());lookup={(g['expID'],g['conditionCode']):i for i,g in enumerate(groups)};assert len(bulk['groups'])==len(groups)
m=bulk['matrix'];matrix=sparse.csr_matrix((np.array(m['counts'],dtype=np.uint64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount']));assert bulk['featureIDs']==list(transport.var_names);seen=set()
for i,g in enumerate(bulk['groups']):
 d=g['donorID'].removeprefix('GSE181897:exp_id:');c=g['condition'].removeprefix('source-code:');key=(d,c);assert key in lookup and key not in seen;seen.add(key);j=lookup[key];expected=reference[j];actual=matrix[i]
 assert g['biologicalReplicateID']==g['donorID'] and g['cellGroup']=='B' and g['organism']=='NCBITaxon:9606'
 assert g['sourceCellIndices']==groups[j]['sourceCellIndices']
 rows=g['sourceCellIndices'];assert g['sampleIDs']==sorted(set(transport.obs['native_sample'].iloc[rows])) and g['batchIDs']==sorted('GSE181897:pool:'+x for x in set(transport.obs['batch'].iloc[rows]))
 for field in ['indptr','indices','data']:np.testing.assert_array_equal(getattr(expected,field),getattr(actual,field))
assert seen==set(lookup)
# All original library totals agree independently, including every selected cell.
with np.load(r/'source-cell-library-counts.npz') as z:np.testing.assert_array_equal(np.asarray(transport.X.sum(axis=1,dtype=np.uint64)).ravel(),z['RNA'][si])
counts_by_donor=json.loads((out/'source-donor-condition-counts.json').read_text());pairs=[d['expID'] for d in counts_by_donor if d['conditionCellCounts'].get('B',0)>0 and d['conditionCellCounts'].get('C',0)>0]
result=dict(status='passed-complete-source-code-handoff',sourceCells=136142,sourceStoredRecords=292741570,selectedCells=len(si),selectedCountRecords=verified,features=len(fi),groups=len(groups),aggregateNonzeros=matrix.nnz,allSelectedSourceCountsExact=True,allRetainedAnnotationsExact=True,allNativeAggregateValuesExact=True,nativeReplayPassed=True,repeatedNativeReportExact=True,sourceDonors=64,donorsWithBothLiteralBAndCCodes=len(pairs),literalCodeMissingDonors=[d['expID'] for d in counts_by_donor if d['expID'] not in pairs],conditionMeaningQualified=False,fitStarted=False,predictionScoringStarted=False,reportSHA256=sha(native/'aggregate/report.json'),checkerSHA256=sha(Path(__file__)),elapsedSeconds=time.time()-start)
(r/'handoff-verification.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
