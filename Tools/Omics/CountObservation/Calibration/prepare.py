"""Bind every original training cell and full RNA denominator before calibration."""
from pathlib import Path
import gzip,hashlib,json,sys,time
import numpy as np
import h5py
root=Path(sys.argv[1]);out=root/'inputs';out.mkdir(exist_ok=False)
prior=Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2,allow_nan=False)+'\n')
protocol={'createdUnix':time.time(),'scope':'Full original Kang and HIRISA B-cell donor-condition training cells; no query treated outcomes','estimator':'Within-donor unbiased normalized-rate and distinct-cell moments; pooled gene-specific NB2 dispersion per condition, Poisson boundary at zero; equal-donor Gamma new-donor-rate prior with estimated measurement variance removed from donor variance and training-mean uncertainty retained','sourceScope':'All 15706 Kang and 18082 HIRISA RNA features; normalization uses all source genes','sourceCells':{'Kang':2651,'HIRISA':119513},'geneFilter':None,'calibrationOutOfDomain':'Unavailable with exact reason; never clamp to supported kernel bounds or replace with a default prior','queryApplication':'All 62 GSE181897 control donors and 16 previously hash-selected panel genes; known development data; unavailable calibration retained','independentBiologicalValidation':False,'modelAssumptions':['Independent cells within each donor-condition','Common gene dispersion across donors within condition','Cell rate independent of RNA depth','Exchangeable donors for rate prior','Estimated hyperparameters treated as fixed'],'referenceTolerance':{'relativeMoments':2e-8,'relativeCalibration':2e-7},'predictorChanged':False}
write(root/'protocol.json',protocol)
sourcepaths={'Kang':Path('/Users/home/numivivo-kang-gamma-product-20260909/original.h5ad'),'HIRISA':Path('/Users/home/numivivo-hiris-20260910/hirisa.h5ad')}
mapping=json.loads((prior/'feature-mapping.json').read_bytes());bindings=json.loads((prior/'source-bindings.json').read_bytes())
paths={'Kang':Path('/Users/home/numivivo-kang-gamma-product-20260909/report.json'),'HIRISA':Path('/Users/home/numivivo-hiris-20260910/inference-inputs/cohorts/02.json.gz')}
def strings(d):return d.asstr()[:].tolist()
def col(f,path):
 d=f[path]
 if isinstance(d,h5py.Dataset):return strings(d)
 if d.attrs.get('encoding-type')=='nullable-string-array':
  assert not np.any(d['mask'][:]);return strings(d['values'])
 labels=np.array(strings(d['categories']),object);codes=d['codes'][:];assert np.all(codes>=0);return labels[codes].tolist()
for origin,source in sourcepaths.items():
 assert sha(paths[origin])==bindings[str(paths[origin])]
 raw=paths[origin].read_bytes();cohort=json.loads(gzip.decompress(raw) if paths[origin].suffix=='.gz' else raw)
 qualified=json.loads((prior/(origin+'-cohort.json')).read_bytes())
 sourcehash=sha(source)
 with h5py.File(source,'r') as h:
  assert h['X'].attrs['encoding-type']=='csr_matrix'
  ids=col(h,'obs/'+h['obs'].attrs['_index']);assert len(ids)==len(set(ids))
  lookup={x:i for i,x in enumerate(ids)}
  cols=col(h,'var/'+h['var'].attrs['_index'])
  fm=sorted([x for x in mapping if x['study']==origin],key=lambda x:x['sourceIndex'])
  assert cols==[x['sourceID'] for x in fm]
  samples=col(h,'obs/'+('native_sample' if origin=='Kang' else 'geo_accession'))
  donor_col=col(h,'obs/geo_donor') if origin=='HIRISA' else None
  treatment_col=col(h,'obs/geo_treatment') if origin=='HIRISA' else None
  groups=[]
  for g in qualified['sourceGroups']:
   local=g['sourceCellIndices'];barcodes=[cohort['metadata']['cells'][i]['barcode'] for i in local]
   rows=[lookup[b] for b in barcodes]
   assert all(samples[i]==g['sampleIDs'][0] for i in rows)
   if origin=='HIRISA':
    assert all(donor_col[i]==g['donorID'] and treatment_col[i]==g['condition'] for i in rows)
   groups.append(dict(donorID=origin+':'+g['donorID'],conditionID='control' if g['condition'] in ['ctrl','none'] else 'IFNB',sourceRows=rows))
  assert sum(len(g['sourceRows']) for g in groups)==protocol['sourceCells'][origin]
  meta={'origin':origin,'sourcePath':str(source),'sourceSHA256':sourcehash,'sourceBytes':source.stat().st_size,'cohortSHA256':sha(paths[origin]),'featureMappingSHA256':sha(prior/'feature-mapping.json'),'featureIDs':[x['transportID'] for x in fm],'groups':groups,'sourceShape':h['X'].attrs['shape'].tolist()}
 write(out/(origin+'.json'),meta)
 print(json.dumps({'origin':origin,'sourceCells':sum(len(g['sourceRows']) for g in groups),'sourceSHA256':sourcehash,'sourceBytes':source.stat().st_size}),flush=True)
write(out/'freeze.json',{'protocolSHA256':sha(root/'protocol.json'),'prepareSHA256':sha(__file__),'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()}})
