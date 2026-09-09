#!/usr/bin/env python3
"""Audit original author-filtered UMI matrices, without reconstructing counts.

This release is distinct from the ineligible normalized pertpy H5AD. No
experimental benchmark result is asserted by a count/metadata inventory.
"""
import argparse,gzip,hashlib,json
from pathlib import Path
import numpy as np
FILES={
'mouse1_unst_filtered_by_cell_cluster0.txt.gz':'6d534e16cc4b4c3a81ad64d93548c2f3ec9742bae45b568c1ad06463204ddf45',
'mouse2_unst_filtered_by_cell_cluster0.txt.gz':'90ef0c28b5112c55ee034565b3f54aa013c470d64b2e3a9d29f1305e0ebb95b5',
'mouse3_unst_filtered_by_cell_cluster0.txt.gz':'41f1c0c65d969f59b2acb17f8d655939e7a534d9dbc92293af0afa2e97cc28b5',
'mouse1_lps6_filtered_by_cell_cluster0.txt.gz':'f3272d159b54e3ef4bbffbb3655fbed762152ef043ebf16a02285f95517c3bd3',
'mouse2_lps6_filtered_by_cell_cluster0.txt.gz':'29ec96593aab0cae453e7e9ae6312ea3edd5d48c9af1c9f5444072bc38fbede1',
'mouse3_lps6_filtered_by_cell_cluster0.txt.gz':'f2305a30632a598843637f077dc1b586ed8998ebc11013e39d8e5b6a7250890c'}
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();results=[]
for name,sha in FILES.items():
    path=a.source_dir/name
    with path.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==sha
    genes=set();axis=hashlib.sha256();nnz=total=0;maximum=0
    with gzip.open(path,'rt') as f:
        barcodes=f.readline().split();assert len(set(barcodes))==len(barcodes)
        for line in f:
            gene,values=line.split(maxsplit=1);assert gene not in genes;genes.add(gene);axis.update((gene+'\n').encode())
            row=np.fromstring(values,sep=' ',dtype=np.float64)
            assert len(row)==len(barcodes) and np.isfinite(row).all() and (row>=0).all() and np.equal(row,np.floor(row)).all()
            assert row.max(initial=0)<=1_000_000,'audit exact-sum bound exceeded'
            nnz+=int(np.count_nonzero(row));total+=int(row.sum());maximum=max(maximum,int(row.max(initial=0)))
    record=dict(file=name,url='https://www.ebi.ac.uk/biostudies/files/E-MTAB-6754/'+name,sha256=sha,
        cells=len(barcodes),features=len(genes),featureAxisSHA256=axis.hexdigest(),nonzeros=nnz,totalUMIs=total,maximumUMI=maximum)
    results.append(record);print(json.dumps(record),flush=True)
assert len({r['featureAxisSHA256'] for r in results})==1
report=dict(status='original-UMI-matrices-audited-not-benchmarked',study='E-MTAB-6754',
    studyURL='https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6754',
    scope='Mouse source groups 1-3, unstimulated versus six-hour LPS; author QC and cluster0 selection retained',
    files=results,cells=sum(r['cells'] for r in results),nonzeros=sum(r['nonzeros'] for r in results),
    donorMapping='Source mouse1, mouse2 and mouse3 prefixes; verify sample metadata when preparing benchmark',
    qualification='Integer UMI/axis inventory only. Native import, count-size capacity, donor inference and expected-biology comparisons remain unqualified.')
with a.out.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
