"""Read complete original sparse cells once; freeze exact depth-binned counts.

The development gene panel is inherited unchanged, including unavailable genes.
All genes still contribute to each cell's original full RNA denominator.
"""
from pathlib import Path
import collections,gzip,hashlib,json,os,sys,time
import h5py,numpy as np
from scipy import sparse
root=Path(sys.argv[1]);origin=sys.argv[2];parent=Path('/Users/home/numivivo-count-calibration-20260911')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8<<20):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
protocol=load(root/'protocol.json');meta=load(parent/'inputs'/(origin+'.json'));query_input=load(parent/'query-input.json.gz');genes=sorted({q['featureID'] for q in query_input['queries']});assert len(genes)==16
source=Path(meta['sourcePath']);state={'status':'running','pid':os.getpid(),'startedUnix':time.time(),'origin':origin,'sourceCells':0,'nonzeros':0,'groupsChecked':0}
def save(): (root/(origin+'-prepare-state.json')).write_text(json.dumps(state,indent=2)+'\n')
save();assert sha(source)==meta['sourceSHA256'];signature=(source.stat().st_ino,source.stat().st_size,source.stat().st_mtime_ns)
report=load(parent/(origin+'-native.json.gz'));assert report['metadataSHA256']==sha(parent/'inputs'/(origin+'.json'))
models={m['conditionID']:{f['featureID']:f for f in m['features']} for m in report['models']}
indices=[meta['featureIDs'].index(g) for g in genes];lookup=np.full(len(meta['featureIDs']),-1,np.int64);lookup[indices]=np.arange(len(genes))
all_depths=[];all_groups=[];all_rows=[];panel_indices=[];panel_values=[];panel_ptr=[0];prepared=[]
try:
 with h5py.File(source,'r') as h:
  assert h['X'].attrs['shape'].tolist()==meta['sourceShape'];ptr=h['X/indptr'][:]
  for gi,group in enumerate(meta['groups']):
   full_counts=np.zeros(len(meta['featureIDs']),np.uint64);bins={};rows=sorted(group['sourceRows']);start=0
   while start<len(rows):
    end=start+1
    while end<len(rows) and end-start<512 and rows[end]==rows[end-1]+1:end+=1
    row0=rows[start];row1=rows[end-1]+1;first=int(ptr[row0]);last=int(ptr[row1]);offsets=ptr[row0:row1+1]-first
    ix=h['X/indices'][first:last].astype(np.int64);yy=h['X/data'][first:last].astype(np.uint64)
    for k,row in enumerate(range(row0,row1)):
     lo,hi=map(int,offsets[k:k+2]);ii=ix[lo:hi];y=yy[lo:hi]
     assert len(ii)>0 and np.all(ii[1:]>ii[:-1]) and np.all(y>0)
     depth=int(y.sum());assert 0<depth<=1_000_000_000
     full_counts[ii]+=y;mask=lookup[ii]>=0;selected=lookup[ii[mask]];counts=y[mask]
     if depth not in bins:bins[depth]=[0,np.zeros(len(genes),np.uint64)]
     bins[depth][0]+=1;bins[depth][1][selected]+=counts
     all_depths.append(depth);all_groups.append(gi);all_rows.append(row);panel_indices.extend(selected.tolist());panel_values.extend(counts.tolist());panel_ptr.append(len(panel_values))
     state['sourceCells']+=1;state['nonzeros']+=len(ii)
    start=end
   assert full_counts.tolist()==report['moments'][gi]['counts']
   depths=sorted(bins);prepared.append({'donorID':group['donorID'],'conditionID':group['conditionID'],'libraryCounts':depths,'cellsPerLibrary':[bins[d][0] for d in depths],
    'geneCounts':[bins[d][1].tolist() for d in depths]})
   assert sum(bins[d][0] for d in depths)==report['moments'][gi]['cells']
   assert sum(d*bins[d][0] for d in depths)==report['moments'][gi]['libraryCounts']
   state['groupsChecked']+=1;save();print(json.dumps(state),flush=True)
 assert signature==(source.stat().st_ino,source.stat().st_size,source.stat().st_mtime_ns)
 assert state['sourceCells']==report['cells'] and state['nonzeros']==report['nonzeros']
 records=[];donors=sorted({g['donorID'] for g in prepared})
 for j,gene in enumerate(genes):
  pair=[]
  for donor in donors:
   data={'donorID':donor}
   for condition,key in [('control','control'),('IFNB','treated')]:
    group=next(g for g in prepared if g['donorID']==donor and g['conditionID']==condition)
    data[key]={'libraryCounts':group['libraryCounts'],'cellsPerLibrary':group['cellsPerLibrary'],'geneCountsPerLibrary':[a[j] for a in group['geneCounts']]}
   pair.append(data)
  c=models['control'][gene];t=models['IFNB'][gene]
  records.append({'featureID':gene,'controlCellDispersion':c.get('cellDispersion'),'treatedCellDispersion':t.get('cellDispersion'),'controlCalibrationStatus':c['status'],'treatedCalibrationStatus':t['status'],'pairs':pair})
 queries=[]
 for q in query_input['queries']:
  bins={}
  for y,depth in zip(q['counts'],q['libraryCounts']):
   if depth not in bins:bins[depth]=[0,0]
   bins[depth][0]+=1;bins[depth][1]+=y
  depths=sorted(bins)
  queries.append({'id':q['id'],'featureID':q['featureID'],'donorID':q['donorID'],'plannedLibraryCounts':q['plannedLibraryCounts'],'control':{'libraryCounts':depths,'cellsPerLibrary':[bins[d][0] for d in depths],'geneCountsPerLibrary':[bins[d][1] for d in depths]}})
 matrix=sparse.csr_matrix((np.array(panel_values,np.uint32),np.array(panel_indices,np.int32),np.array(panel_ptr,np.int64)),shape=(len(all_depths),len(genes)));matrix.sort_indices()
 np.savez_compressed(root/(origin+'-individual-cells.npz'),data=matrix.data,indices=matrix.indices,indptr=matrix.indptr,shape=np.array(matrix.shape),libraryCounts=np.array(all_depths,np.uint32),groupIndices=np.array(all_groups,np.int32),sourceRows=np.array(all_rows,np.int32),featureIDs=np.array(genes))
 refs={'groups':[{'donorID':g['donorID'],'conditionID':g['conditionID']} for g in meta['groups']],'queries':query_input['queries']}
 with (root/(origin+'-reference-identities.json.gz')).open('xb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(json.dumps(refs,sort_keys=True,separators=(',',':')).encode())
 bindings={'sourceSHA256':meta['sourceSHA256'],'parentMetadataSHA256':sha(parent/'inputs'/(origin+'.json')),'parentCalibrationSHA256':sha(parent/(origin+'-native.json.gz')),'queryInputSHA256':sha(parent/'query-input.json.gz'),'protocolSHA256':sha(root/'protocol.json'),'individualCellsSHA256':sha(root/(origin+'-individual-cells.npz'))}
 obj={'origin':origin,'genes':records,'queries':queries,'gridSizes':protocol['gridRefinements'],'sourceBindings':bindings}
 raw=json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
 with (root/(origin+'-input.json.gz')).open('xb') as f,gzip.GzipFile(filename='',fileobj=f,mode='wb',mtime=0) as g:g.write(raw)
 state.update(status='completed',inputSHA256=hashlib.sha256(raw).hexdigest(),inputCompressedSHA256=sha(root/(origin+'-input.json.gz')),sourceBindings=bindings,genes=genes,depthBins=sum(len(g['libraryCounts']) for g in prepared),individualPanelNonzeros=int(matrix.nnz))
except BaseException as e:state.update(status='failed',error=repr(e));raise
finally:state.update(finishedUnix=time.time(),seconds=time.time()-state['startedUnix']);save()
