from remote_io import *
import numpy as np
r=json.loads((ROOT/'sources/roster.json').read_text());codes=np.load(ROOT/'sources/roster-codes.npz');selected=codes['selected'];cats=r['categories']
expected={key:[] for key in ['donor','cytokine','treatment']}
for sample_name in cats['sample']:
 donor,condition=sample_name.split('_',1)
 expected['donor'].append(cats['donor'].index(donor));expected['cytokine'].append(cats['cytokine'].index(condition));expected['treatment'].append(cats['treatment'].index('PBS' if condition=='PBS' else 'cytokine'))
assert len(np.unique(codes['sample']))==1092
for key in expected:assert np.array_equal(np.array(expected[key])[codes['sample']],codes[key])
assert len(set(r['sourceFeatures']))==40352
barcodes=[];nonzeros=[];totals=[];samples=[];axes=[];group_ids=[];genes=[]
for i in range(3456):
 p=ROOT/f'axes/run-{i:04d}.npz';a=np.load(p);lo,hi=int(a['lo']),int(a['hi']);assert lo==len(barcodes) and np.array_equal(a['sourceRows'],selected[lo:hi]);n=np.diff(a['indptr'])
 assert (n>=0).all() and (a['tscp_count']>=0).all()
 genes+=a['gene_count'].tolist()
 barcodes+=a['barcodes'].tolist();nonzeros+=n.tolist();totals+=a['tscp_count'].tolist()
 sid=np.unique(codes['sample'][selected[lo:hi]]);assert len(sid)==1 # original contiguous runs have one source sample
 group_ids.append(int(sid[0]));axes.append(dict(run=i,lo=lo,hi=hi,sourceFirst=int(selected[lo]),sourceLast=int(selected[hi-1]),entryFirst=int(a['indptr'][0]),entryLast=int(a['indptr'][-1]),sample=cats['sample'][int(sid[0])],axisSHA256=hashlib.sha256(p.read_bytes()).hexdigest()))
assert len(barcodes)==len(set(barcodes))==725031
used=sorted(set(codes['sample'][selected].tolist()));assert len(used)==24
for s in used:
 sid=cats['sample'][s];donor,condition=sid.split('_',1);samples.append(dict(id=sid,biologicalReplicateID=donor,donorID=donor,condition=condition,batchID='Parse-original-release',organism='Homo sapiens'))
metadata=dict(id='Parse10M-IFNB-PBS-all-released-cells',evidence='measured',sourceDescription='Original Parse 10M PBMC release, complete literal IFN-beta/PBS cohort. CC BY-NC 4.0, Parse Biosciences. All original RNA features; no new cell filtering. Source annotations are not new author QC or replicate qualification.',countUnit='umiCount',samples=samples,features=[dict(id=x,name=x,mitochondrial=False) for x in r['sourceFeatures']],cells=[dict(barcode=b,sampleID=cats['sample'][int(s)]) for b,s in zip(barcodes,codes['sample'][selected])])
plan=dict(schemaVersion=1,metadata=metadata,sourceDeclaration=json.dumps(dict(url=URL,ETag=ETAG,bytes=SIZE,selection='Every source row whose literal cytokine is PBS or IFN-beta; original row indices and range hashes retained in companion evidence; no full-file SHA256 claimed',license='CC BY-NC 4.0',dose='Unresolved; no prediction scoring',rowTotals='Original gene_count differs from retained X cardinality; historical tscp_count not declared as an exact X total. Actual streamed totals will be independently compared',mitochondrialAnnotation='Not supplied; mitochondrial fractions unavailable'),sort_keys=True),rowNonzeros=nonzeros)
out=ROOT/'prepared';out.mkdir(exist_ok=True);data=json.dumps(plan,separators=(',',':')).encode();(out/'plan.json').write_bytes(data)
(out/'runs.json').write_text(json.dumps(axes,indent=2)+'\n');(out/'freeze.json').write_text(json.dumps(dict(status='count-admission-only',cells=len(barcodes),features=len(r['sourceFeatures']),groups=len(samples),nonzeros=sum(nonzeros),sourceHistoricalTscpCounts=sum(totals),sourceGeneCountMismatches=int(np.count_nonzero(np.array(nonzeros)!=np.array(genes))),planSHA256=hashlib.sha256(data).hexdigest(),sourceETag=ETAG,sourceBytes=SIZE,predictionFitOrScoring=False,filtering='None beyond literal condition selection; all 12 donors and both conditions'),indent=2)+'\n');print((out/'freeze.json').read_text())
