#!/usr/bin/env python3
"""Extract exact native QL residual and abundance inputs without refitting counts."""
import argparse,gzip,hashlib,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--abundance-root',type=Path,required=True)
parser.add_argument('--ql-root',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args();sha=lambda raw:hashlib.sha256(raw).hexdigest()
runs=json.loads((args.abundance_root/'family/complete.json').read_text());assert len(runs)==58
args.out.mkdir(parents=True,exist_ok=True);receipts=[]
for run in runs:
    assert run['status']=='passed'
    path=args.abundance_root/'family'/run['case']/run['method']/'native.json.gz'
    raw=path.read_bytes();assert sha(raw)==run['outputSHA256']
    n=json.loads(gzip.decompress(raw))['fit'];assert n['completed'] and not n['failures']
    original=args.ql_root/run['case']/'input.json.gz';raw=original.read_bytes();assert sha(raw)==run['inputSHA256']
    feature_ids=json.loads(gzip.decompress(raw))['featureIndices']
    residuals=n['globalFit']['adjustedResiduals'];abundance=n['abundanceFits']
    assert len(residuals)==len(abundance)==len(feature_ids)
    assert all(r.get('quasiDispersion') is not None for r in residuals)
    payload=dict(featureIndices=feature_ids,variances=[r['quasiDispersion'] for r in residuals],
                 degreesOfFreedom=[r['degreesOfFreedom'] for r in residuals],
                 abundance=[r['log2CountsPerMillion'] for r in abundance])
    data=gzip.compress(json.dumps(payload,separators=(',',':'),allow_nan=False).encode(),mtime=0)
    dest=args.out/run['case']/run['method'];dest.mkdir(parents=True,exist_ok=True)
    if (dest/'input.json.gz').exists():assert (dest/'input.json.gz').read_bytes()==data
    else:(dest/'input.json.gz').write_bytes(data)
    receipt=dict(case=run['case'],method=run['method'],genes=len(feature_ids),inputSHA256=sha(data),
                 upstreamNativePath=str(path),upstreamNativeSHA256=run['outputSHA256'],
                 upstreamCountPath=str(original),upstreamCountSHA256=run['inputSHA256'],
                 minimumDF=min(payload['degreesOfFreedom']),maximumDF=max(payload['degreesOfFreedom']),
                 zeroVariances=sum(x==0 for x in payload['variances']))
    (dest/'input-receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');receipts.append(receipt)
    print('extracted',run['case'],run['method'],flush=True)
(args.out/'complete.json').write_text(json.dumps(receipts,sort_keys=True,indent=2)+'\n')
