#!/usr/bin/env python3
"""Freeze disjoint original-donor cohorts from metadata, without reading expression values."""
import argparse, hashlib, json
from pathlib import Path
import h5py
from anndata.io import read_elem

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root
sha=lambda s:hashlib.sha256(s.encode()).hexdigest()
protocol=Path(__file__).with_name('PROTOCOL.md')
source=json.loads((root/'source-receipt.json').read_text())
assert source['sha256']=='2f019ba5ae48fdf61e02f510aac66cf6035b44aa4c582d52efce9ebccad90ee3'
assert not (root/'cohorts.json').exists()
with h5py.File(root/'source.h5ad','r') as h:
    obs=read_elem(h['obs']);var=read_elem(h['raw/var'])
assert len(obs)==160632 and len(var)==32357 and obs.index.is_unique and var.index.is_unique
assert set(obs.AIFI_L1.astype(str))=={'B cell'}
donors=[]
for donor,g in obs.groupby('donor_id',observed=True):
    assert not g[['batch_id','sample.sampleKitGuid','sample.visitName']].isna().any().any()
    assert all(g[k].nunique()==1 for k in ['batch_id','sample.sampleKitGuid','sample.visitName','subject.ageGroup','sex','subject.cmv'])
    donors.append(dict(donorID=str(donor),batchID=str(g.batch_id.iloc[0]),
        sampleKitID=str(g['sample.sampleKitGuid'].iloc[0]),visit=str(g['sample.visitName'].iloc[0]),
        ageGroup=str(g['subject.ageGroup'].iloc[0]),sex=str(g.sex.iloc[0]),cmv=str(g['subject.cmv'].iloc[0]),
        cells=len(g),subtypeCounts={str(k):int(v) for k,v in g.AIFI_L2.value_counts().items()}))
assert len(donors)==108 and len({d['sampleKitID'] for d in donors})==108
donors.sort(key=lambda d:(d['batchID'],sha('numivivo-independent-null-v1|cohort|'+d['donorID']),d['donorID']))
cohorts=[];samples=[]
for i in range(9):
    cohortID=f'{i+1:02d}';members=donors[i*12:(i+1)*12];batches={}
    for d in members:batches.setdefault(d['batchID'],[]).append(d)
    odds=sorted((b for b,ds in batches.items() if len(ds)%2),key=lambda b:(sha('numivivo-independent-null-v1|odd|'+cohortID+'|'+b),b))
    assert len(odds)%2==0
    for batch,ds in batches.items():
        ds.sort(key=lambda d:(sha('numivivo-independent-null-v1|arm|'+d['donorID']),d['donorID']))
        start=odds.index(batch)%2 if batch in odds else 0
        for j,d in enumerate(ds):
            d['condition']=['shamA','shamB'][(start+j)%2];d['cohortID']=cohortID
            samples.append(dict(id=d['donorID'],donorID=d['donorID'],biologicalReplicateID=d['donorID'],
                batchID=d['batchID'],condition=d['condition'],organism='NCBITaxon:9606'))
    assert sum(d['condition']=='shamA' for d in members)==6
    cohorts.append(dict(id=cohortID,donors=members))
mapping=dict(schemaVersion=1,id='human-immune-health-atlas-independent-null',evidence='measured',
    sourceDescription='Complete CELLxGENE B/Plasma release b6986a7f-981e-4e04-93ed-57f444749b8b; original UMI counts; sham conditions assigned by frozen metadata-only donor randomization; no biological intervention.',
    countUnit='umiCount',matrixPath='raw/X',samples=sorted(samples,key=lambda d:d['id']),sampleColumn='donor_id',groupColumn='AIFI_L1',
    featureNameColumn='gene_symbols',mitochondrialFeatureIDs=[str(x) for x in var.index[var.mito.astype(bool)]])
plan=dict(schemaVersion=1,mapping=mapping,contrasts=[])
manifest=dict(source=source,protocolSHA256=hashlib.sha256(protocol.read_bytes()).hexdigest(),cohorts=cohorts,
    assignmentUsesExpression=False,sourceCells=len(obs),sourceFeatures=len(var))
(root/'cohorts.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
(root/'stream-plan.json').write_text(json.dumps(plan,sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(cohorts=9,donors=len(donors),cells=len(obs),features=len(var),
    designResidualDFBeforeCellFilter={c['id']:12-(len({d['batchID'] for d in c['donors']})+1) for c in cohorts}),indent=2))
