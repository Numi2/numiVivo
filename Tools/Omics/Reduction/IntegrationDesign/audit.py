from pathlib import Path
import json,gzip,hashlib,collections
import numpy as np,h5py,anndata as ad
s=Path(__file__).parent;repo=Path('/Users/home/numivivo-expression-memory-20260911');base=repo/'Tools/Omics/Benchmarks/HIRISA/evidence/2026-09-10-full-clustering';retained=Path('/Users/home/numivivo-cluster-export-20260912');manifest=json.loads((retained/'source-manifest.json').read_text());records={v['storedPath']:v for v in manifest['records']}
def frozen(name):
 p=retained/name;b=p.read_bytes();rec=records[name];assert hashlib.sha256(b).hexdigest()==rec['storedSHA256'];decoded=gzip.decompress(b);assert hashlib.sha256(decoded).hexdigest()==rec['sourceSHA256'];return json.loads(decoded)
plan=frozen('run/pca-plan.json.gz');receipt=frozen('bundle/input/input/receipt.json.gz');assert bytes(receipt['plan']['bytes']).hex()==records['run/pca-plan.json.gz']['sourceSHA256'];source=Path('/Users/home/numivivo-hiris-20260910/hirisa.h5ad');h=hashlib.sha256()
with source.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
assert h.hexdigest()==bytes(receipt['source']['bytes']).hex();samples=plan['mapping']['samples'];lookup={v['id']:dict(v) for v in samples};assert len(lookup)==len(samples)==131
with h5py.File(source,'r') as f:
 accession=ad.io.read_elem(f['obs/'+plan['mapping']['sampleColumn']]);enrichment=ad.io.read_elem(f['obs/geo_enrichment']);assert len(accession)==len(enrichment)==1612594
 counts=collections.Counter();labels=collections.defaultdict(set)
 for a,e in zip(accession,enrichment):counts[str(a)]+=1;labels[str(a)].add(str(e))
assert set(counts)==set(lookup)
for sid,v in lookup.items():assert len(labels[sid])==1;v['cells']=counts[sid];v['sourceEnrichment']=next(iter(labels[sid]))
def matrix(rows,fields):
 columns=['intercept']+[f+'='+level for f in fields for level in sorted({r[f] for r in rows})[1:]]
 x=np.array([[1]+[int(r[f]==level) for f in fields for level in sorted({z[f] for z in rows})[1:]] for r in rows],dtype=float);return x,columns
def rank_qr(x):
 # Independent column-pivoted modified Gram-Schmidt, with reorthogonalization.
 a=x.copy();rank=0
 for k in range(min(a.shape)):
  pivot=k+int(np.argmax(np.sum(a[:,k:]**2,axis=0)));a[:,[k,pivot]]=a[:,[pivot,k]];norm=np.linalg.norm(a[:,k])
  if norm<1e-10:break
  q=a[:,k]/norm
  for _ in range(2):a[:,k+1:]-=np.outer(q,q@a[:,k+1:])
  rank+=1
 return rank
scopes={'all':list(lookup.values())}
for label in sorted({v['sourceEnrichment'] for v in lookup.values()}):scopes['enrichment:'+label]=[v for v in lookup.values() if v['sourceEnrichment']==label]
result={'scope':'Source-metadata design audit, not integration execution or biological preservation qualification. Unpenalized categorical design rank does not establish separate causal effects. Ridge can make a redundant solve unique without making its effect decomposition data-identifiable.','sourceSHA256':h.hexdigest(),'planSHA256':records['run/pca-plan.json.gz']['sourceSHA256'],'samples':list(lookup.values()),'scopes':{}}
for name,rows in scopes.items():
 designs={}
 for fields in [('donorID',),('batchID',),('donorID','batchID'),('donorID','batchID','condition'),('condition','sourceEnrichment'),('donorID','batchID','sourceEnrichment'),('donorID','batchID','condition','sourceEnrichment')]:
  x,columns=matrix(rows,fields);singular=np.linalg.svd(x,compute_uv=False);rank=int(np.sum(singular>1e-10));assert rank==rank_qr(x);weighted=x*np.sqrt(np.array([v['cells'] for v in rows]))[:,None];assert np.linalg.matrix_rank(weighted)==rank
  designs['+'.join(fields)]={'rows':len(rows),'columns':columns,'rank':rank,'columnCount':len(columns),'rankDeficiency':len(columns)-rank,'singularValues':singular.tolist(),'cellWeightedRankMatches':True}
 nuisance=designs['donorID+batchID'];full=designs['donorID+batchID+condition'];conditions=sorted({v['condition'] for v in rows});result['scopes'][name]={'libraries':len(rows),'cells':sum(v['cells'] for v in rows),'donors':len({v['donorID'] for v in rows}),'batches':len({v['batchID'] for v in rows}),'conditions':conditions,'designs':designs,'conditionRankIncrement':full['rank']-nuisance['rank'],'allConditionMainEffectDimensionsEstimable':full['rank']-nuisance['rank']==len(conditions)-1,'batchesWithMultipleDonors':[b for b in sorted({v['batchID'] for v in rows}) if len({v['donorID'] for v in rows if v['batchID']==b})>1]}
 protected=designs['condition+sourceEnrichment'];combined=designs['donorID+batchID+condition+sourceEnrichment']
 result['scopes'][name]['protectedConditionEnrichmentDimensions']=protected['rank']-1
 result['scopes'][name]['protectedDimensionsAfterNuisance']=combined['rank']-nuisance['rank']
 result['scopes'][name]['protectedDimensionsOverlappingNuisance']=protected['rank']-1-(combined['rank']-nuisance['rank'])
 aliases={}
 for label in sorted({v['sourceEnrichment'] for v in rows}):
  batches=sorted({v['batchID'] for v in rows if v['sourceEnrichment']==label})
  if len({v['sourceEnrichment'] for v in rows})>1 and all((v['sourceEnrichment']==label)==(v['batchID'] in batches) for v in rows):aliases[label]=batches
 result['scopes'][name]['exactEnrichmentIndicatorBatchAliases']=aliases
 print(name,{k:v for k,v in result['scopes'][name].items() if k!='designs'},'nuisance rank',nuisance['rank'],'/',nuisance['columnCount'],flush=True)
result['status']='pass-source-bindings-all-library-memberships-and-dual-rank-checks';(s/'design-audit.json').write_text(json.dumps(result,indent=2)+'\n')
