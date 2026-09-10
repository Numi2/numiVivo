#!/usr/bin/env python3
"""Re-read hashed source CSR independently and verify each frozen reference input."""
import argparse,gzip,hashlib,json,math,subprocess
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
sha=lambda b:hashlib.sha256(b).hexdigest();protocol=json.loads((root/'protocol.json').read_text());results=[]
for rec in protocol['cases']:
 raw=subprocess.check_output(['ssh',protocol['host'],'cat',rec['source']]);assert sha(raw)==rec['sourceSHA256'];expanded=gzip.decompress(raw);assert sha(expanded)==rec['logicalSHA256'];s=json.loads(expanded);old=s['contrasts'][0];design=old['design']
 path=root/rec['id'];receipt=json.loads((path/'input-receipt.json').read_text());raw=(path/'input.json.gz').read_bytes();assert sha(raw)==receipt['inputSHA256'];i=json.loads(gzip.decompress(raw));indices=[g['featureIndex'] for g in old['features'] if g['status']=='tested'];assert indices==i['featureIndices'];assert [s['metadata']['features'][g]['id'] for g in indices]==i['featureIDs']
 m=s['pseudobulk']['matrix'];matrix=csr_matrix((np.asarray(m['counts'],dtype=np.uint64),m['featureIndices'],m['rowOffsets']),shape=(len(s['pseudobulk']['groups']),len(s['metadata']['features'])))
 expected=matrix[np.asarray(design['sourcePseudobulkIndices'])][:,indices].T.toarray();assert np.array_equal(expected,np.asarray(i['counts'],dtype=np.uint64));assert i['design']==design['rows'] and i['contrast']==design['contrast']
 offset=np.asarray(i['offsets']);original=np.log(design['sizeFactorValues']);constant=float(np.mean(np.log(design['libraryCounts'])));error=float(np.max(np.abs(offset-original-constant)));assert error<1e-12
 assert i['nativeTrend']==[old['negativeBinomial']['features'][g]['trendDispersion'] for g in indices]
 result=dict(case=rec['id'],status='passed',sourceSHA256=rec['sourceSHA256'],inputSHA256=receipt['inputSHA256'],genes=len(indices),countEntries=int(expected.size),offsetShiftMaximumError=error);results.append(result);print(rec['id'],'passed',flush=True)
(root/'input-checks.json').write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
