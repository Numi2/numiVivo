#!/usr/bin/env python3
"""Independently reconstruct GO descriptors and Jaccard geometry from captured JSON."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re
import numpy as np

SOURCE='e70fcd49808cab8d724de8d5a332940911206e1c8ef44cc7b568d048ed795c85'

def digest(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()

def verify(root,identities,expected_receipt):
    assert digest(root/'receipt.json')==expected_receipt
    receipt=json.loads((root/'receipt.json').read_text())
    assert receipt['sourceSHA256']==SOURCE and digest(identities)==receipt['identitiesSHA256']
    assert receipt['controlsVerified'] is False and receipt['experimentalGuideAssignmentsVerified'] is False
    assert receipt['expressionOutcomesRead'] is False
    listed=[]
    for entry in receipt['files']:
        relative=Path(entry['path'])
        assert not relative.is_absolute() and '..' not in relative.parts
        path=root/relative
        assert not path.is_symlink() and path.is_file()
        assert digest(path)==entry['sha256'] and path.stat().st_size==entry['bytes']
        listed.append(entry['path'])
    assert len(listed)==len(set(listed))
    assert set(listed)=={str(f.relative_to(root)) for f in root.rglob('*') if f.is_file() and f.name!='receipt.json'}
    builds=[json.loads((root/name).read_text())['build_version'] for name in ('metadata-before.json','metadata-after.json')]
    assert builds==[receipt['buildVersion']]*2
    identity=json.loads(identities.read_text()); by_symbol=defaultdict(list)
    for i,name in enumerate(identity['names']): by_symbol[name].append(identity['features'][i])
    prefixes=defaultdict(list);unknown=[]
    for guide in identity['groups']:
        parts=guide.rsplit('_',1)
        if len(parts)==2 and re.fullmatch('[A-Za-z0-9-]+',parts[0]) and re.fullmatch('p(DS|BA)[0-9]+',parts[1]): prefixes[parts[0]].append(guide)
        else:unknown.append(guide)
    candidates=json.loads((root/'candidates.json').read_text())
    assert candidates['unresolvedSourceLabels']==unknown
    assert candidates['experimentalAssignmentsVerified'] is False
    requests=json.loads((root/'requests.json').read_text())
    assert [r['target'] for r in requests]==sorted(prefixes)
    rebuilt={}; coverage=[]; overlap=[]
    for request in requests:
        target=request['target']; assert len(by_symbol[target])==1
        gene=by_symbol[target][0]
        assert request['sourceFeatureID']==gene and request['sourceGuides']==sorted(prefixes[target])
        assert request['assignment']=='unverified-guide-prefix'
        item={k:request[k] for k in ['target','sourceFeatureID','sourceGuides','assignment','status']}
        if request['status']=='annotation-not-found':
            assert request['httpStatus']==404
        else:
            assert request['status']=='retrieved'
            raw=json.loads((root/request['response']).read_text())
            assert request['responseSHA256']==digest(root/request['response'])
            fields=raw.get('ensembl',[]); fields=[fields] if isinstance(fields,dict) else fields
            valid=raw.get('taxid')==9606 and gene in [f.get('gene') for f in fields]
            retained=[]; excluded=0
            if valid:
                for group in ['BP','MF','CC']:
                    records=raw.get('go',{}).get(group,[])
                    records=[records] if isinstance(records,dict) else records
                    for record in records:
                        q=record.get('qualifier',''); q=q if isinstance(q,list) else q.split('|')
                        if 'NOT' in q or record['id'] in ['GO:0008150','GO:0003674','GO:0005575']:
                            excluded+=1;continue
                        assert re.fullmatch('GO:[0-9]{7}',record['id'])
                        retained.append({**record,'namespace':group})
            distinct=sorted({r['id'] for r in retained})
            status=('supported' if distinct else 'no-usable-GO-annotations') if valid else 'identity-mismatch'
            item.update(status=status,currentSymbol=raw.get('symbol'),terms=len(distinct),excludedAnnotations=excluded)
            if distinct:rebuilt[target]=dict(terms=distinct,annotations=retained)
            cited=[]
            for record in retained:
                publications=record.get('pubmed',[])
                publications=publications if isinstance(publications,list) else [publications]
                if any(str(p)=='27984733' for p in publications):cited.append(record)
            if cited:overlap.append(dict(target=target,sourceFeatureID=gene,annotations=cited))
        coverage.append(item)
    assert coverage==json.loads((root/'coverage.json').read_text())
    assert rebuilt==json.loads((root/'annotations.json').read_text())
    assert overlap==json.loads((root/'study-citation-overlap.json').read_text())
    supported=sorted(rebuilt); vocabulary=sorted({t for d in rebuilt.values() for t in d['terms']})
    sets=[set(rebuilt[t]['terms']) for t in supported]
    pairwise=np.array([[len(a & b)/len(a | b) for b in sets] for a in sets])
    binary=np.array([[term in labels for term in vocabulary] for labels in sets],dtype=np.int64)
    intersection=binary@binary.T; sizes=binary.sum(axis=1)
    independent=intersection/(sizes[:,None]+sizes[None,:]-intersection)
    np.testing.assert_array_equal(pairwise,independent)
    minimum=float(np.linalg.eigvalsh(pairwise).min());assert minimum>=-1e-10
    off_diagonal=pairwise.copy();np.fill_diagonal(off_diagonal,0)
    no_shared=[supported[i] for i in np.flatnonzero(off_diagonal.max(axis=1)==0)]
    assert len(supported)==receipt['supportedCandidates'] and len(vocabulary)==receipt['terms']
    assert len(overlap)==receipt['candidatesWithStudyCitation']
    assert sum(len(x['annotations']) for x in overlap)==receipt['studyCitationAnnotations']
    return dict(status='passed-capture-reconstruction',receiptSHA256=expected_receipt,files=len(listed),
        candidates=len(requests),supportedCandidates=len(supported),terms=len(vocabulary),
        annotations=sum(len(x['annotations']) for x in rebuilt.values()),
        independentJaccardExact=True,minimumJaccardEigenvalue=minimum,
        candidatesWithoutSharedTerms=no_shared,studyCitationOverlap=overlap,
        unsupported=[dict(target=x['target'],status=x['status']) for x in coverage if x['status']!='supported'],
        controlsVerified=False,experimentalGuideAssignmentsVerified=False,predictionQualification=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--capture',type=Path,required=True);p.add_argument('--identities',type=Path,required=True)
    p.add_argument('--receipt-sha256',required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();assert not a.out.exists()
    result=verify(a.capture,a.identities,a.receipt_sha256)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n');print(a.out.read_text())
