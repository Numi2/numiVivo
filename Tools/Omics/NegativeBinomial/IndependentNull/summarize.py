#!/usr/bin/env python3
"""Retain complete null families, failures, call inventories and joint-family comparisons."""
import argparse, collections, csv, gzip, hashlib, json, math
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--available',action='store_true')
a=p.parse_args();root=a.root;rows=[];families={};sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
load=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
def number(v):
    try:r=float(v);return r if math.isfinite(r) else None
    except (TypeError,ValueError):return None
def bh(values):
    v=np.array(values);n=len(v)
    if not n:return []
    order=np.argsort(v,kind='stable');q=np.empty(n);q[order]=np.minimum(1,np.minimum.accumulate((v[order]*n/np.arange(1,n+1))[::-1])[::-1]);return q.tolist()
def inventory(table):
    p=[r['p'] for r in table.values() if r['p'] is not None];q=[r['q'] for r in table.values() if r['q'] is not None]
    return dict(tested=len(p),adjusted=len(q),rawPBelow001=sum(x<.01 for x in p),rawPBelow005=sum(x<.05 for x in p),
        rawPFractionBelow001=sum(x<.01 for x in p)/len(p) if p else None,rawPFractionBelow005=sum(x<.05 for x in p)/len(p) if p else None,
        bh001=sum(x<=.01 for x in q),bh005=sum(x<=.05 for x in q),bh010=sum(x<=.10 for x in q),anyBH005=any(x<=.05 for x in q))
for number_ in range(1,10):
    cohort=f'{number_:02d}';meta=load(root/'reference-inputs'/cohort/'input.json')
    for variant,methods in [('native',['wald','lrt','ql']),('native-fixed',['wald','lrt'])]:
        for method in methods:
            directory=root/variant/(cohort+'-'+method)
            if not (directory/'receipt.json').exists():assert a.available;continue
            receipt=load(directory/'receipt.json');assert sha(directory/'output.json.gz')==receipt['outputSHA256']
            out=load(directory/'output.json.gz');result=out.get('result');name=variant+'-'+method
            table={f['featureID']:dict(p=f.get('pValue'),q=f.get('adjustedPValue'),effect=f.get('log2FoldChange')) for f in result['features']} if result else {}
            families[(cohort,name)]=table
            diagnostics=result['negativeBinomial']['features'] if result else []
            row=dict(cohort=cohort,method=name,attempted=meta['attemptedFeatures'],eligible=meta['eligibleFeatures'],
                withheld=meta['attemptedFeatures']-sum(v['p'] is not None for v in table.values()),
                statuses=dict(collections.Counter(f['status'] for f in result['features'])) if result else {},
                errors=dict(collections.Counter(d['error'] for d in diagnostics if d.get('error'))),
                geneWiseLowerBoundaries=sum(d.get('geneWiseLowerBoundary',False) for d in diagnostics),
                geneWiseUpperBoundaries=sum(d.get('geneWiseUpperBoundary',False) for d in diagnostics),
                methodError=out.get('error'),returnCode=receipt['returnCode'],seconds=receipt['seconds'],
                familyFailures=result['negativeBinomial'].get('quasiLikelihood',{}).get('failures',[]) if result else [],
                outputSHA256=receipt['outputSHA256'],**inventory(table))
            rows.append(row)
    reference=root/'reference'/cohort;runs=load(reference/'runs.json');receipt=load(root/'reference'/(cohort+'-receipt.json'))
    for filename,digest in receipt['files'].items():assert sha(reference/filename)==digest
    for method,suffix in [('edgeR-QL','edgeR-QL'),('limma-voom','limma-voom'),('DESeq2','DESeq2'),('DESeq2-default','DESeq2-default-results')]:
        filename='native_size_factors-'+suffix+'.tsv';path=reference/filename;table={}
        if path.exists():
            with path.open() as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    table[row['featureID']]=dict(p=number(row.get('pValue',row.get('pvalue'))),
                        q=number(row.get('adjustedPValue',row.get('padj'))),effect=number(row.get('log2FoldChange')))
        run=runs['native_size_factors-'+('DESeq2' if method=='DESeq2-default' else method)]
        name='reference-'+method;families[(cohort,name)]=table
        rows.append(dict(cohort=cohort,method=name,attempted=meta['attemptedFeatures'],eligible=meta['eligibleFeatures'],
            withheld=meta['attemptedFeatures']-sum(v['p'] is not None for v in table.values()),returnCode=receipt['returnCode'],
            methodStatus=run['status'],methodError=run.get('error'),warnings=run['warnings'],messages=run['messages'],
            seconds=run['seconds'],**inventory(table)))
joint=[]
for cohort in [f'{i:02d}' for i in range(1,10)]:
    for variant in ['native','native-fixed']:
        names=[variant+'-wald',variant+'-lrt','native-ql','reference-edgeR-QL','reference-limma-voom','reference-DESeq2']
        if not all((cohort,name) in families for name in names):continue
        genes=sorted(set.intersection(*[{g for g,v in families[(cohort,name)].items() if v['p'] is not None} for name in names]))
        for name in names:
            subset={g:dict(families[(cohort,name)][g]) for g in genes}
            # A separate BH family, never silently replace the original full-family q.
            for g,q in zip(genes,bh([subset[g]['p'] for g in genes])):subset[g]['q']=q
            joint.append(dict(cohort=cohort,comparison=variant,method=name,jointFamilySize=len(genes),**inventory(subset)))
events={}
for method in sorted({r['method'] for r in rows}):
    selected=[r for r in rows if r['method']==method]
    events[method]=dict(cohorts=len(selected),cohortsWithBH005=sum(r['anyBH005'] for r in selected),
        descriptiveFraction=sum(r['anyBH005'] for r in selected)/len(selected),totalBH005Calls=sum(r['bh005'] for r in selected),
        tested=sum(r['tested'] for r in selected))
complete=len(rows)==81
summary=dict(status='complete-inventory' if complete else 'partial-inventory',rows=len(rows),events=events,
    scope='Nine disjoint donor cohorts in one selected study; sham labels, no biological intervention. Descriptive event fractions, not nine studies, universal FDR, power, interval coverage or production qualification.',
    sourceCheck=load(root/'source-check.json'),originalProtocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),
    followupSHA256=sha(Path(__file__).with_name('FOLLOWUP.md')))
prefix='available-' if a.available else ''
for name,value in [('summary.json',summary),('families.json',rows),('joint-families.json',joint)]:
    (root/(prefix+name)).write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
columns=['cohort','method','attempted','eligible','tested','withheld','rawPBelow001','rawPBelow005','bh001','bh005','bh010','anyBH005']
with (root/(prefix+'families.tsv')).open('w') as f:
    w=csv.DictWriter(f,fieldnames=columns,delimiter='\t',extrasaction='ignore',lineterminator='\n');w.writeheader();w.writerows(rows)
print(json.dumps(summary,indent=2))
