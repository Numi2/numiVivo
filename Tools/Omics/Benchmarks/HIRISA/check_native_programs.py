#!/usr/bin/env python3
"""Compare every native program coordinate, detection and total with frozen RNA."""
import argparse,json,resource,platform,time
from pathlib import Path
import numpy as np
from check_integration_response import sha,blocks


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ['bundle','reference','metadata','out']:p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();assert not a.out.exists();began=time.monotonic()
 receipt_sha=sha(a.bundle/'receipt.json');receipt=json.loads((a.bundle/'receipt.json').read_text());receipt={k:bytes(v['bytes']).hex() if isinstance(v,dict) and set(v)=={'bytes'} else v for k,v in receipt.items()};model=json.loads((a.bundle/'model.json').read_text());reference=json.loads((a.reference/'reference.json').read_text());plan=json.loads((a.bundle/'plan.json').read_text())
 assert receipt['schemaVersion']==1 and receipt['matrixFormat']=='complete-row-major-u32-row-u32-program-f64-le/v1'
 assert receipt['integerFormat']=='complete-row-major-u32-row-u32-column-u64-le/v1'
 assert receipt['source']==reference['bindings']['source'] and receipt['metadata']==sha(a.metadata)
 assert plan['matchFeatureNames'] and plan['normalizationTarget']==10000
 for n,k in [('metadata.json','metadata'),('model.json','model'),('scores.bin','scores'),('detected-members.bin','detectedMembers'),('total-counts.bin','totalCounts'),('plan.json','plan')]:assert sha(a.bundle/n)==receipt[k],n
 assert sha(a.bundle/'original.h5ad')==receipt['source']
 for n,identity in reference['payloads'].items():assert (a.reference/n).stat().st_size==identity['bytes'] and sha(a.reference/n)==identity['SHA256']
 assert model['cells']==reference['cells']==1612594 and model['canonicalNonzeros']==reference['entries']
 assert model['sourcePasses']==2 and len(model['programs'])==len(reference['programs'])==2
 for native,measured in zip(model['programs'],reference['programs']):
  assert native['definition']==measured['definition']
  assert native['featureIndices']==measured['featureIndices'] and native['missingFeatureIDs']==measured['missingSymbols']
  np.testing.assert_allclose(native['effectiveWeights'],measured['effectiveWeights'],rtol=1e-14,atol=1e-16)
  assert abs(native['weightCoverage']-measured['weightCoverage'])<1e-14
 scores=np.load(a.reference/'scores.npy',mmap_mode='r');detected=np.load(a.reference/'detected.npy',mmap_mode='r');totals=np.load(a.reference/'total-counts.npy',mmap_mode='r')
 assert model['updates']==int(detected.sum())
 maximum=0.
 for first,last,values in blocks(a.bundle/'scores.bin',len(scores),2,8192):
  expected=np.asarray(scores[first:last]);valid=np.isfinite(expected)
  maximum=max(maximum,float(np.max(np.abs(values[valid]-expected[valid]))))
  np.testing.assert_allclose(values[valid],expected[valid],rtol=1e-12,atol=1e-12)
  assert (values[~valid]==0).all()
 dtype=np.dtype([('row','<u4'),('column','<u4'),('value','<u8')])
 for name,expected,width in [('detected-members.bin',detected,2),('total-counts.bin',totals,1)]:
  path=a.bundle/name;assert path.stat().st_size==len(scores)*width*16
  records=np.memmap(path,dtype=dtype,mode='r')
  for first in range(0,len(scores),8192):
   last=min(first+8192,len(scores));v=records[first*width:last*width]
   np.testing.assert_array_equal(v['row'],np.repeat(np.arange(first,last),width))
   np.testing.assert_array_equal(v['column'],np.tile(np.arange(width),last-first))
   np.testing.assert_array_equal(v['value'].reshape(last-first,width),np.asarray(expected[first:last]).reshape(last-first,width))
 assert sha(a.bundle/'original.h5ad')==receipt['source']
 for n,identity in reference['payloads'].items():assert sha(a.reference/n)==identity['SHA256']
 assert sha(a.bundle/'receipt.json')==receipt_sha
 # Recheck all payload identities after numerical reads.
 for n,k in [('metadata.json','metadata'),('model.json','model'),('scores.bin','scores'),('detected-members.bin','detectedMembers'),('total-counts.bin','totalCounts'),('plan.json','plan')]:assert sha(a.bundle/n)==receipt[k]
 result=dict(status='passed',cells=len(scores),programs=2,maximumScoreError=maximum,allCellMetadataExact=True,allDetectionsExact=True,allLibraryTotalsExact=True,definitionsAndResolvedIndicesExact=True,receiptSHA256=sha(a.bundle/'receipt.json'),referenceSHA256=sha(a.reference/'reference.json'),checkerSHA256=sha(Path(__file__)),seconds=time.monotonic()-began,maximumResidentBytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if platform.system()=='Darwin' else 1024),scope='Every original cell/program against separately frozen SciPy/RNA scores, exact detections and totals; original IDs and definitions retained. Native replay separate. Descriptive expression arithmetic, not biological calibration.')
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
