#!/usr/bin/env python3
"""Freeze exact full-support source families and extract donor pseudobulk inputs."""
import argparse,gzip,hashlib,json,math,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--lrt-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
sha=lambda b:hashlib.sha256(b).hexdigest()
source=json.loads((a.lrt_root/'protocol.json').read_text())
cases=[c for c in source['cases'] if (c['family']=='null' and c['id'].endswith('-default')) or (c['family']=='treatment' and c['id']!='kang-active')];assert len(cases)==29
assert not (a.out/'protocol.json').exists();a.out.mkdir(parents=True,exist_ok=True)
protocol=dict(host=source['host'],cases=cases,protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md').read_bytes()),sourceProtocolSHA256=sha((a.lrt_root/'protocol.json').read_bytes()))
(a.out/'protocol.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n')
for c in cases:
 raw=subprocess.check_output(['ssh',source['host'],'cat',c['source']]);assert sha(raw)==c['sourceSHA256'];raw=gzip.decompress(raw);assert sha(raw)==c['logicalSHA256'];s=json.loads(raw)
 old=s['contrasts'][0];design=old['design'];assert not old['request'].get('negativeBinomialOptions',{}).get('zeroTotalDonorPolicy')
 features=[f for f in old['features'] if f['status']=='tested'];idx=[f['featureIndex'] for f in features];positions={g:j for j,g in enumerate(idx)}
 csr=s['pseudobulk']['matrix'];counts=[]
 for i in design['sourcePseudobulkIndices']:
  row=[0]*len(idx)
  for k in range(csr['rowOffsets'][i],csr['rowOffsets'][i+1]):
   g=csr['featureIndices'][k]
   if g in positions:row[positions[g]]=csr['counts'][k]
  counts.append(row)
 libconstant=sum(math.log(n) for n in design['libraryCounts'])/len(design['libraryCounts'])
 offsets=[math.log(v)+libconstant for v in design['sizeFactorValues']]
 diag=old['negativeBinomial']['features'];lrt=json.loads((a.lrt_root/c['id']/'check.json').read_text())
 data=dict(case=c['id'],featureIndices=idx,featureIDs=[f['featureID'] for f in features],counts=list(map(list,zip(*counts))),design=design['rows'],contrast=design['contrast'],offsets=offsets,offsetConstant=libconstant,nativeTrend=[diag[i]['trendDispersion'] for i in idx],nativeWaldCalls=lrt['originalBHCalls'],nativeLRTCalls=lrt['likelihoodRatioBHCalls'],sourceSHA256=c['sourceSHA256'])
 d=a.out/c['id'];d.mkdir();packed=gzip.compress(json.dumps(data,allow_nan=False,separators=(',',':')).encode(),mtime=0);(d/'input.json.gz').write_bytes(packed)
 (d/'input-receipt.json').write_text(json.dumps(dict(inputSHA256=sha(packed),logicalSHA256=sha(gzip.decompress(packed)),sourceSHA256=c['sourceSHA256'],genes=len(idx),observations=len(counts),columns=len(design['contrast'])),sort_keys=True,indent=2)+'\n')
 print(c['id'],len(idx),flush=True)
