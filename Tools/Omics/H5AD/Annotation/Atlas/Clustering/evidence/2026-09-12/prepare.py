import gzip,json,hashlib,pathlib,h5py,anndata as ad,numpy as np
s=pathlib.Path(__file__).parent
source=pathlib.Path('/Users/home/numivivo-hiris-20260910/hirisa.h5ad')
def read(p): return gzip.decompress((s/p).read_bytes())
def sha(b): return hashlib.sha256(b).hexdigest()
manifest=json.loads((s/'source-manifest.json').read_text())
for rec in manifest['records']:
 p=s/rec['storedPath']
 if p.exists():
  assert sha(p.read_bytes())==rec['storedSHA256']
  assert sha(gzip.decompress(p.read_bytes()))==rec['sourceSHA256']
r=json.loads(read('bundle/result.json.gz'));receipt=json.loads(read('bundle/receipt.json.gz'));graph=json.loads(read('bundle/input/receipt.json.gz'));pca=json.loads(read('bundle/input/input/receipt.json.gz'))
def fp(x):return bytes(x['bytes']).hex()
assert fp(receipt['result'])==sha(read('bundle/result.json.gz'))
assert fp(receipt['input'])==sha(read('bundle/input/receipt.json.gz'))
assert fp(graph['input'])==sha(read('bundle/input/input/receipt.json.gz'))
h=hashlib.sha256()
with source.open('rb') as f:
 for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
assert h.hexdigest()==fp(pca['source'])
with h5py.File(source,'r') as f:
 assert fp(pca['plan'])==sha(read('run/pca-plan.json.gz'))
 mapping=json.loads(read('run/pca-plan.json.gz'))['mapping']
 ids=f['obs/_index'].asstr()[:]; samples=ad.io.read_elem(f['obs/'+mapping['sampleColumn']])
 assert len(r['cells'])==len(ids)==len(r['labels'])==1612594
 assert all(c['barcode']==ids[i] and c['sampleID']==samples[i] for i,c in enumerate(r['cells']))
assert np.bincount(r['labels']).tolist()==r['clusterSizes']
plan={'schemaVersion':1,'source':pca['source'],'provenance':'Native graph communities, not authoritative cell types or biological outcome predictions. All 1612594 source row barcode/sample identities independently matched in order. Clustering result SHA256 '+fp(receipt['result'])+'; clustering receipt SHA256 '+sha(read('bundle/receipt.json.gz'))+'.','edits':[{'path':'obs/numivivo_graph_community','mode':'add','value':{'int64':{'shape':[len(ids)],'values':r['labels']}}}]}
(s/'plan.json').write_text(json.dumps(plan,separators=(',',':'))+'\n')
(s/'identity-verification.json').write_text(json.dumps({'status':'PASS','cells':len(ids),'communities':len(r['clusterSizes']),'sourceSHA256':h.hexdigest(),'resultSHA256':fp(receipt['result']),'allBarcodeSamplePairsMatchInOrder':True,'receiptChainVerified':True},indent=2)+'\n')
print((s/'identity-verification.json').read_text())
