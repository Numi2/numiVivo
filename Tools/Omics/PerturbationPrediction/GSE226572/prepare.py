#!/usr/bin/env python3
"""Validate every released sparse record, apply frozen QC, prepare full native inputs."""
import argparse,json,sys,time
from pathlib import Path
import anndata as ad,h5py,numpy as np,pandas as pd
from scipy import sparse
from download import sha,write

def strings(dataset):return np.asarray(dataset.asstr()[:],dtype=str)
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--kang',type=Path,required=True);a=p.parse_args();root=a.study;out=root/'prepared';out.mkdir(exist_ok=False);qc=out/'qc';qc.mkdir()
 freeze=json.loads((root/'source-roster-freeze.json').read_text());assert sha(Path(__file__).with_name('PROTOCOL.md'))==freeze['protocolSHA256']
 downloads=json.loads((root/'downloads.json').read_text());assert downloads['status']=='completed' and len(downloads['files'])==24
 blocks=[];observations=[];summaries=[];sample_defs=[];library_counts=[];all_features=None;row_cursor=0
 for rec in downloads['files']:
  path=root/'sources'/rec['path'];assert sha(path)==rec['SHA256'];started=time.time()
  with h5py.File(path,'r') as f:
   g=f['matrix'];feature_ids=strings(g['features/id']);names=strings(g['features/name']);types=strings(g['features/feature_type']);genomes=strings(g['features/genome']);barcodes=strings(g['barcodes']);shape=tuple(map(int,g['shape'][:]));ptr=g['indptr'][:].astype(np.int64)
   assert shape==(len(feature_ids),len(barcodes)) and len(set(feature_ids))==len(feature_ids) and len(set(barcodes))==len(barcodes)
   assert len(ptr)==len(barcodes)+1 and ptr[0]==0 and ptr[-1]==len(g['data'])==len(g['indices']) and np.all(np.diff(ptr)>=0)
   features=dict(ids=feature_ids.tolist(),names=names.tolist(),types=types.tolist(),genomes=genomes.tolist())
   if all_features is None:all_features=features
   else:assert features==all_features,'Query source feature axes differ'
   rna=types=='Gene Expression';assert rna.all(),'Undeclared additional modality requires admission review'
   mito=np.char.startswith(names,'MT-');assert mito.any();n=len(barcodes);totals=np.zeros(n,dtype=np.uint64);detected=np.zeros(n,dtype=np.uint32);mito_totals=np.zeros(n,dtype=np.uint64);selected=np.zeros(n,dtype=bool);pieces=[];source_records=0
   for start in range(0,n,4096):
    end=min(n,start+4096);lo,hi=map(int,(ptr[start],ptr[end]));values=g['data'][lo:hi];indices=g['indices'][lo:hi]
    assert np.all(np.isfinite(values)&(values>=0)&(values==np.floor(values))) and values.max(initial=0)<2**32
    assert np.all((indices>=0)&(indices<len(feature_ids)))
    x=sparse.csr_matrix((values.astype(np.uint32),indices,ptr[start:end+1]-lo),shape=(end-start,len(feature_ids)))
    assert x.has_canonical_format and x.has_sorted_indices
    totals[start:end]=np.asarray(x.sum(axis=1,dtype=np.uint64)).ravel();detected[start:end]=np.asarray((x>0).sum(axis=1)).ravel();mito_totals[start:end]=np.asarray(x[:,mito].sum(axis=1,dtype=np.uint64)).ravel()
    # Integer inequality includes exactly20 percent; no floating-point boundary ambiguity.
    assert totals[start:end].max(initial=0)<2**60
    keep=(detected[start:end]>=750)&(mito_totals[start:end]*5<=totals[start:end]);selected[start:end]=keep
    if keep.any():pieces.append(x[keep])
    source_records+=len(values)
   assert selected.any(),'Frozen QC leaves an unsupported library'
   counts=sparse.vstack(pieces,format='csr');assert counts.shape[0]==int(selected.sum())
   np.testing.assert_array_equal(np.asarray(counts.sum(axis=1,dtype=np.uint64)).ravel(),totals[selected])
   np.savez_compressed(qc/(rec['accession']+'.npz'),barcodes=barcodes,totalCounts=totals,detectedGenes=detected,mitochondrialCounts=mito_totals,selected=selected)
   bulk=np.asarray(counts.sum(axis=0,dtype=np.uint64)).ravel();library_counts.append(bulk);blocks.append(counts)
   observations.append(pd.DataFrame(dict(native_sample=[rec['accession']]*len(counts.indptr[:-1]),population=['QC-PBMC']*counts.shape[0],source_barcode=barcodes[selected]),index=pd.Index([rec['accession']+'|'+b for b in barcodes[selected]],dtype=object)))
   donor='GSE226572:'+rec['donor'];condition='control' if rec['hours']==0 else 'IFNB:'+str(rec['hours'])+'h'
   sample_defs.append(dict(id=rec['accession'],donorID=donor,biologicalReplicateID=donor,condition=condition,batchID=rec['accession'],organism='NCBITaxon:9606'))
   summaries.append(dict(rec,sourceBarcodes=n,sourceRecords=source_records,sourceUMIs=int(totals.sum()),selectedCells=int(selected.sum()),selectedRecords=counts.nnz,selectedUMIs=int(bulk.sum()),transportStart=row_cursor,transportEnd=row_cursor+counts.shape[0],below750Genes=int((detected<750).sum()),above20PercentMito=int((mito_totals*5>totals).sum()),seconds=time.time()-started));row_cursor+=counts.shape[0]
   print(json.dumps({k:summaries[-1][k] for k in ['accession','sourceBarcodes','selectedCells','sourceRecords']}),flush=True)
 assert all_features is not None
 counts=sparse.vstack(blocks,format='csr');obs=pd.concat(observations);assert counts.shape[0]==len(obs)==row_cursor and obs.index.is_unique
 var=pd.DataFrame(dict(source_symbol=all_features['names'],feature_type=all_features['types'],genome=all_features['genomes']),index=pd.Index(all_features['ids'],dtype=object))
 ad.AnnData(counts,obs=obs,var=var).write_h5ad(out/'query-full-QC.h5ad',compression='gzip');np.savez_compressed(out/'query-library-counts.npz',counts=np.vstack(library_counts),featureIDs=np.array(all_features['ids']),featureNames=np.array(all_features['names']),sampleIDs=np.array([r['accession'] for r in downloads['files']]))
 mapping=dict(schemaVersion=1,id='GSE226572-complete-initial-QC-PBMC',evidence='measured',sourceDescription='All24rawGEOlibraries and allbarcodes passing declared initial count QC only. Not author curated clusters, validated singlets or authoritative cell identities. External10x-to-H5AD conversion; source symbols and Ensembl IDs retained.',countUnit='umiCount',matrixPath='X',sampleColumn='native_sample',groupColumn='population',featureNameColumn='source_symbol',samples=sample_defs)
 write(out/'query-aggregate.json',dict(schemaVersion=1,mapping=mapping,contrasts=[]));write(out/'query-source-audit.json',dict(status='passed-all-source-sparse-records',files=summaries,sourceBarcodes=sum(s['sourceBarcodes'] for s in summaries),sourceRecords=sum(s['sourceRecords'] for s in summaries),admittedCells=row_cursor,admittedRecords=counts.nnz,featureCount=len(all_features['ids']),mitochondrialIDs=np.array(all_features['ids'])[np.char.startswith(np.array(all_features['names']),'MT-')].tolist(),preparerSHA256=sha(Path(__file__))))
 # Use every previously admitted Kang cell; cell-type labels do not select rows.
 assert sha(a.kang)=='1d48c1ff10bcfad7c1bc75863ce702a491c6c71b08e921a5e9754a986a9d19f3'
 kang=ad.read_h5ad(a.kang);assert len(kang)==24673;kang.obs['population']='QC-PBMC';kang.write_h5ad(out/'kang-full.h5ad',compression='gzip')
 # Reuse the previously source-qualified sample definitions, including literal donor/condition fields.
 previous=Path('/Users/home/numivivo-reference-cross-study-20260911/inputs/kang-query.json');km=json.loads(previous.read_text())['mapping'];km['id']='Kang-all-admitted-PBMC';km['groupColumn']='population';km['sourceDescription']='Every24673source-admittedKangcell, without cell-type selection. Full source features and RNA denominators. Previously inspected training cohort.'
 write(out/'kang-aggregate.json',dict(schemaVersion=1,mapping=km,contrasts=[]))
 unique=sorted(set(kang.obs['native_sample'].astype(str)));kb=np.vstack([np.asarray(kang.X[np.asarray(kang.obs['native_sample'].astype(str)==s)].sum(axis=0,dtype=np.uint64)).ravel() for s in unique]);np.savez_compressed(out/'kang-library-counts.npz',counts=kb,featureIDs=np.asarray(kang.var_names,dtype=str),featureNames=np.asarray(kang.var_names,dtype=str),sampleIDs=np.array(unique))
 write(out/'input-freeze.json',dict(status='completed-preparation-before-native',protocolSHA256=freeze['protocolSHA256'],downloadsSHA256=sha(root/'downloads.json'),originalKangSHA256=sha(a.kang),files={str(f.relative_to(out)):sha(f) for f in sorted(out.rglob('*')) if f.is_file()},fitStarted=False,scoringStarted=False));print('Complete source preparation frozen',flush=True)
if __name__=='__main__':main()
