#!/usr/bin/env python3
"""Run native whole-cohort aggregation and verify source handoff independently."""
import argparse,json,os,re,subprocess,time
from pathlib import Path
import anndata as ad,h5py,numpy as np
from scipy import sparse
from download import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();root=a.study;src=root/'prepared';out=root/'native-aggregation';out.mkdir(exist_ok=False)
 freeze=json.loads((src/'input-freeze.json').read_text())
 for name,digest in freeze['files'].items():assert sha(src/name)==digest,name
 execution=dict(status='running',binarySHA256=sha(a.binary),hdf5SHA256=sha(Path(os.environ['NUMIVIVO_HDF5_LIBRARY'])),inputFreezeSHA256=sha(src/'input-freeze.json'),runnerSHA256=sha(Path(__file__)),commands=[]);write(out/'execution.json',execution)
 for cohort,file in [('query','query-full-QC.h5ad'),('kang','kang-full.h5ad')]:
  for op,args in [('aggregate',['singlecell-h5ad-pseudobulk',src/file,'--plan',src/(cohort+'-aggregate.json'),'--output',out/cohort]),('verify',['singlecell-h5ad-pseudobulk-verify',out/cohort])]:
   started=time.time();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(out/(cohort+'-'+op+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr);execution['commands'].append(dict(cohort=cohort,operation=op,exitCode=r.returncode,args=list(map(str,args)),seconds=time.time()-started,maximumResidentBytes=int(peak[1]) if peak else None));write(out/'execution.json',execution);assert r.returncode==0,(cohort,op,r.stderr)
  report=json.loads((out/cohort/'report.json').read_text());bulk=report['pseudobulk'];m=bulk['matrix'];actual=sparse.csr_matrix((np.array(m['counts'],dtype=np.uint64),m['featureIndices'],m['rowOffsets']),shape=(m['cellCount'],m['featureCount'])).toarray()
  with np.load(src/(cohort+'-library-counts.npz'),allow_pickle=False) as z:
   assert z['featureIDs'].tolist()==bulk['featureIDs'];ids=z['sampleIDs'].tolist();counts=z['counts']
   for row,g in enumerate(bulk['groups']):
    indices=[ids.index(s) for s in g['sampleIDs']];np.testing.assert_array_equal(actual[row],counts[indices].sum(axis=0,dtype=np.uint64));assert g['cellGroup']=='QC-PBMC'
  source=ad.read_h5ad(src/file,backed='r');assert len(report['metadata']['cells'])==len(source)
  for i,c in enumerate(report['metadata']['cells']):assert c['barcode']==source.obs_names[i] and c['sampleID']==str(source.obs['native_sample'].iloc[i])
  source.file.close();print(json.dumps(dict(cohort=cohort,cells=len(report['metadata']['cells']),groups=len(bulk['groups']),allAggregatesExact=True)),flush=True)
 # Independent original-source versus consolidated H5AD comparison, every admitted record.
 audit=json.loads((src/'query-source-audit.json').read_text());transport=ad.read_h5ad(src/'query-full-QC.h5ad',backed='r');checked=0
 for rec in audit['files']:
  with np.load(src/'qc'/(rec['accession']+'.npz'),allow_pickle=False) as q:mask=q['selected'];source_barcodes=q['barcodes']
  target_row=rec['transportStart']
  with h5py.File(root/'sources'/rec['path'],'r') as f:
   g=f['matrix'];np.testing.assert_array_equal(g['barcodes'].asstr()[:],source_barcodes);np.testing.assert_array_equal(g['features/id'].asstr()[:],transport.var_names);ptr=g['indptr'][:]
   for start in range(0,len(mask),4096):
    end=min(len(mask),start+4096);keep=mask[start:end]
    if not keep.any():continue
    lo,hi=map(int,(ptr[start],ptr[end]));x=sparse.csr_matrix((g['data'][lo:hi],g['indices'][lo:hi],ptr[start:end+1]-lo),shape=(end-start,transport.n_vars))[keep];y=transport.X[target_row:target_row+x.shape[0],:]
    np.testing.assert_array_equal(x.indptr,y.indptr);np.testing.assert_array_equal(x.indices,y.indices);np.testing.assert_array_equal(x.data,y.data);np.testing.assert_array_equal(transport.obs_names[target_row:target_row+x.shape[0]],[rec['accession']+'|'+v for v in source_barcodes[start:end][keep]])
    checked+=x.nnz;target_row+=x.shape[0]
  assert target_row==rec['transportEnd']
 assert checked==audit['admittedRecords'];transport.file.close()
 execution['status']='passed-native-aggregation-source-handoff-and-reconstruction';execution['allSelectedSourceRecordsExact']=checked;execution['allNativeAggregatesAndCellIdentitiesExact']=True;write(out/'execution.json',execution)
 write(out/'aggregation-freeze.json',dict(status='completed',executionSHA256=sha(out/'execution.json'),files={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()},fitStarted=False));print(json.dumps(dict(status=execution['status'],sourceRecordsCompared=checked)),flush=True)
if __name__=='__main__':main()
