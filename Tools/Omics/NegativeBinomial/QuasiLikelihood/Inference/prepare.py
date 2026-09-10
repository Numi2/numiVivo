#!/usr/bin/env python3
"""Bind real counts and all expected native upstream stages for QL inference."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--abundance-root',type=Path,required=True)
p.add_argument('--moderation-root',type=Path,required=True)
p.add_argument('--ql-root',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();sha=lambda b:hashlib.sha256(b).hexdigest()
runs=json.loads((a.abundance_root/'family/complete.json').read_text())
prior={(r['case'],r['method']):r for r in json.loads((a.moderation_root/'family-final/complete.json').read_text())}
assert len(runs)==len(prior)==58
a.out.mkdir(parents=True,exist_ok=True);receipts=[]
for run in runs:
    assert run['status']=='passed';key=(run['case'],run['method']);r=prior[key];assert r['status']=='passed'
    relative=Path(*key)
    raw=(a.abundance_root/'family'/relative/'native.json.gz').read_bytes();assert sha(raw)==run['outputSHA256']
    native=json.loads(gzip.decompress(raw))['fit'];global_fit=native['globalFit'];assert native['completed']
    raw=(a.ql_root/run['case']/'input.json.gz').read_bytes();assert sha(raw)==run['inputSHA256']
    original=json.loads(gzip.decompress(raw))
    raw=(a.moderation_root/'family-final'/relative/'native.json.gz').read_bytes();assert sha(raw)==r['nativeSHA256']
    moderation=json.loads(gzip.decompress(raw));assert moderation['featureIndices']==original['featureIndices']
    payload=dict(featureIndices=original['featureIndices'],counts=original['counts'],design=original['design'],
        offsets=original['offsets'],contrast=original['contrast'],
        trendDispersions=[f['dispersion'] for f in global_fit['initialFits']],
        expectedRefits=global_fit['refittedFits'],expectedModeration=moderation['fit'],
        expectedResidualDF=[f['degreesOfFreedom'] for f in global_fit['adjustedResiduals']],
        expectedAbundance=[f['log2CountsPerMillion'] for f in native['abundanceFits']])
    out=a.out/relative;assert not out.exists();out.mkdir(parents=True)
    raw=gzip.compress(json.dumps(payload,separators=(',',':'),allow_nan=False).encode(),mtime=0)
    (out/'input.json.gz').write_bytes(raw)
    receipt=dict(case=key[0],method=key[1],genes=len(original['featureIndices']),inputSHA256=sha(raw),
        upstreamCountSHA256=run['inputSHA256'],upstreamAbundanceSHA256=run['outputSHA256'],
        upstreamModerationSHA256=r['nativeSHA256'])
    (out/'input-receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');receipts.append(receipt)
    print('prepared',*key,receipt['genes'],flush=True)
(a.out/'complete.json').write_text(json.dumps(receipts,sort_keys=True,indent=2)+'\n')
