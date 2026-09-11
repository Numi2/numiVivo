"""Verify retained bytes, repeatability, and all donor/gene diagnostic arithmetic."""
import argparse,gzip,hashlib,json
from pathlib import Path,PurePosixPath
import numpy as np
from scipy.stats import t
p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);a=p.parse_args();r=a.evidence
manifest=json.loads((r/'manifest.json').read_bytes())
actual={str(p.relative_to(r)) for p in r.rglob('*') if p.is_file()}
assert actual=={'manifest.json',*manifest['files']}
for name,item in manifest['files'].items():
 path=PurePosixPath(name);assert not path.is_absolute() and '..' not in path.parts
 b=(r/name).read_bytes();assert len(b)==item['bytes'] and hashlib.sha256(b).hexdigest()==item['SHA256']
arrays=np.load(r/'results/moments.npz',allow_pickle=False)
assert arrays['control'].shape==arrays['treated'].shape==(62,11800)
models=json.loads(gzip.decompress((r/'model-moments.json.gz').read_bytes()))
summary=json.loads((r/'results/summary.json').read_bytes());records=json.loads((r/'results/donors.json').read_bytes())
assert len(records)==124
observed=arrays['treated']-arrays['control'];noise=arrays['controlVariance']+arrays['treatedVariance']
assert np.all(noise>=0)
for origin,model in models.items():
 assert model['featureIDs']==arrays['panel'].tolist()
 mean=np.array(model['meanResponse']);v=np.array([np.nan if x is None else x for x in model['donorResponseVariances']]);available=np.isfinite(v)
 critical=float(t.ppf(.975,model['donors']-1))
 np.testing.assert_allclose(critical,model['studentCriticalValue'],rtol=1e-11,atol=1e-11)
 base=v*(1+1/model['donors']);half=critical*np.sqrt(base)
 residual=observed-mean
 s=next(x for x in summary if x['origin']==origin)
 np.testing.assert_allclose(s['rawResponseMSE'],np.mean(residual**2),rtol=0,atol=1e-14)
 np.testing.assert_allclose(s['squaredDonorAverageResidual'],np.mean(np.mean(residual,axis=0)**2),rtol=0,atol=1e-14)
 np.testing.assert_allclose(s['rawResponseMSE'],s['squaredDonorAverageResidual']+s['betweenDonorResidualVariance'],rtol=0,atol=1e-12)
 for di,donor in enumerate(arrays['donors'].tolist()):
  record=next(x for x in records if x['origin']==origin and x['donor']==donor)
  control,truth=arrays['control'][di],arrays['treated'][di]
  wide=critical*np.sqrt(base+noise[di]);assert np.all(wide[available]>=half[available])
  original=np.mean((np.abs(residual[di])<=half)[available])
  raw=np.mean((np.abs(residual[di])<=wide)[available])
  treated=np.mean(((truth>=np.maximum(0,control+mean-wide))&(truth<=np.maximum(0,control+mean+wide)))[available])
  np.testing.assert_allclose([record['originalRawCoverage'],record['observedCellWideningRawCoverage'],record['observedCellWideningTreatedCoverage']],[original,raw,treated],rtol=0,atol=1e-14)
repeat=json.loads((r/'repeat-check.json').read_bytes())
for item in repeat['records']:
 b=(r/'results'/item['path']).read_bytes();assert hashlib.sha256(b).hexdigest()==item['SHA256'] and item['repeatByteIdentical']
for name in ['pipeline.json']:
 one=json.loads((r/'results'/name).read_bytes());two=json.loads((r/'repeat-pipeline.json').read_bytes())
 assert one['scriptSHA256']==two['scriptSHA256']==hashlib.sha256((r/'diagnose.py').read_bytes()).hexdigest()
 assert one['independentBootstrapCountChecks']==two['independentBootstrapCountChecks']==248
print(json.dumps(dict(status='passed',files=len(manifest['files']),donors=62,features=11800,originDonorComparisons=124,
 completeCellBootstrapRepeats=2,newPredictor=False,independentBiologicalValidation=False)))
