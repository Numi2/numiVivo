#!/usr/bin/env python3
"""Capture fixed direct GO descriptors without reading expression outcomes.

Guide prefixes are candidates until experimental guide/control identities are
verified. This command neither fits a model nor declares those identities proven.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import re
import urllib.error
import urllib.parse
import urllib.request

from prepare_ingestion import SOURCE, sha, write

ROOTS={'GO:0008150','GO:0003674','GO:0005575'}
PMID=27984733

def get(url):
    with urllib.request.urlopen(url,timeout=30) as response: raw=response.read(4_000_001)
    if len(raw)>4_000_000: raise ValueError('response byte limit')
    json.loads(raw)
    return raw

def candidates(identity):
    features=defaultdict(list)
    for symbol,feature in zip(identity['names'],identity['features']): features[symbol].append(feature)
    guides=defaultdict(list); unresolved=[]
    for label in identity['groups']:
        match=re.fullmatch(r'([A-Za-z0-9-]+)_p(?:DS|BA)[0-9]+',label)
        if match: guides[match.group(1)].append(label)
        else: unresolved.append(label)
    records=[]
    for target in sorted(guides):
        ids=features[target]
        record=dict(target=target,sourceGuides=sorted(guides[target]),assignment='unverified-guide-prefix')
        if len(ids)==1: record.update(sourceFeatureID=ids[0],status='resolved-source-feature')
        else: record.update(status='unresolved-source-symbol')
        records.append(record)
    return records,unresolved

def terms(value,gene):
    ensembl=value.get('ensembl',[])
    if isinstance(ensembl,dict): ensembl=[ensembl]
    ids={e['gene'] for e in ensembl if isinstance(e,dict) and isinstance(e.get('gene'),str)}
    if value.get('taxid')!=9606 or gene not in ids: return 'identity-mismatch',[],[],[]
    included=[]; excluded=[]
    for namespace in ('BP','MF','CC'):
        records=value.get('go',{}).get(namespace,[])
        if isinstance(records,dict): records=[records]
        for entry in records:
            qualifier=entry.get('qualifier','')
            qualifiers=qualifier if isinstance(qualifier,list) else qualifier.split('|')
            record=dict(entry,namespace=namespace)
            if 'NOT' in qualifiers or entry.get('id') in ROOTS:
                excluded.append(record); continue
            if not re.fullmatch(r'GO:[0-9]{7}',entry.get('id','')): raise ValueError('invalid GO identifier')
            included.append(record)
    distinct=sorted({e['id'] for e in included})
    return ('supported' if distinct else 'no-usable-GO-annotations'),distinct,included,excluded

def mentions_study(entry):
    values=entry.get('pubmed',[])
    if not isinstance(values,list): values=[values]
    return str(PMID) in [str(v) for v in values]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--identities',type=Path,required=True)
    p.add_argument('--ingestion-manifest',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    manifest=json.loads(a.ingestion_manifest.read_text())
    expected=next(e['logicalSHA256'] for e in manifest['files'] if e['logicalPath']=='prepared/identities.json')
    assert sha(a.identities)==expected
    assert next(e['sha256'] for e in manifest['externalFiles'] if e['path']=='source/original.h5ad')==SOURCE
    identity=json.loads(a.identities.read_text()); records,unresolved=candidates(identity)
    assert len(records)==90 and all(r['status']=='resolved-source-feature' for r in records)
    a.out.mkdir(parents=True,exist_ok=False); (a.out/'responses').mkdir()
    write(a.out/'candidates.json',dict(candidates=records,unresolvedSourceLabels=unresolved,
        originalSHA256=SOURCE,identitiesSHA256=expected,experimentalAssignmentsVerified=False))
    before=get('https://mygene.info/v3/metadata');(a.out/'metadata-before.json').write_bytes(before)
    def fetch(record):
        gene=record['sourceFeatureID'];url='https://mygene.info/v3/gene/'+urllib.parse.quote(gene)+'?fields=symbol,taxid,ensembl.gene,entrezgene,go'
        request=dict(record,requestURL=url,retrievedUTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
        try: raw=get(url)
        except urllib.error.HTTPError as e:
            if e.code==404:return dict(request,status='annotation-not-found',httpStatus=404)
            raise
        path=a.out/'responses'/(gene+'.json');path.write_bytes(raw)
        return dict(request,status='retrieved',response=str(path.relative_to(a.out)),responseSHA256=sha(path),bytes=len(raw))
    with ThreadPoolExecutor(max_workers=4) as pool: requests=list(pool.map(fetch,records))
    write(a.out/'requests.json',requests)
    after=get('https://mygene.info/v3/metadata');(a.out/'metadata-after.json').write_bytes(after)
    build=json.loads(before)['build_version'];assert build==json.loads(after)['build_version']
    coverage=[];annotations={};overlap=[]
    for request in requests:
        record={k:request[k] for k in ['target','sourceFeatureID','sourceGuides','assignment','status']}
        if request['status']=='retrieved':
            path=a.out/request['response'];assert sha(path)==request['responseSHA256']
            value=json.loads(path.read_text())
            status,distinct,included,excluded=terms(value,request['sourceFeatureID'])
            record.update(status=status,currentSymbol=value.get('symbol'),terms=len(distinct),excludedAnnotations=len(excluded))
            if distinct: annotations[request['target']]=dict(terms=distinct,annotations=included)
            cited=[e for e in included if mentions_study(e)]
            if cited: overlap.append(dict(target=request['target'],sourceFeatureID=request['sourceFeatureID'],annotations=cited))
        coverage.append(record)
    write(a.out/'coverage.json',coverage);write(a.out/'annotations.json',annotations);write(a.out/'study-citation-overlap.json',overlap)
    files=[dict(path=str(f.relative_to(a.out)),sha256=sha(f),bytes=f.stat().st_size) for f in sorted(a.out.rglob('*')) if f.is_file()]
    write(a.out/'receipt.json',dict(schemaVersion=1,status='captured-descriptors-experimental-assignments-unverified',
        sourceSHA256=SOURCE,identitiesSHA256=expected,ingestionManifestSHA256=sha(a.ingestion_manifest),
        implementationSHA256=sha(Path(__file__)),buildVersion=build,candidates=len(records),
        supportedCandidates=len(annotations),terms=len(set(t for x in annotations.values() for t in x['terms'])),
        studyPMID=PMID,candidatesWithStudyCitation=len(overlap),studyCitationAnnotations=sum(len(x['annotations']) for x in overlap),
        expressionOutcomesRead=False,controlsVerified=False,experimentalGuideAssignmentsVerified=False,files=files))
    print(json.dumps({k:v for k,v in json.loads((a.out/'receipt.json').read_text()).items() if k!='files'},indent=2))

if __name__=='__main__':main()
