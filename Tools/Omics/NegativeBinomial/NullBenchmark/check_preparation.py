#!/usr/bin/env python3
"""Independent original-source split reconstruction and sparse count validation."""
import argparse,hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
from scipy import sparse
from prepare import SOURCES,PROTOCOL
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--source-root',type=Path,default=Path('/Users/home'));p.add_argument('--out',type=Path,required=True)
a=p.parse_args();assert not a.out.exists();checks=[]
for study in ['kang','hagai']:
 d=a.root/study;source=a.source_root/('numivivo-'+study+'-r-native-current-20260909')/'original.h5ad'
 with source.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==SOURCES[study]
 obj=ad.read_h5ad(source);old=json.loads(source.with_name('plan.json').read_text());mapping=old['mapping'];control=old['contrasts'][0]['controlCondition']
 original=obj.obs[mapping['sampleColumn']].astype(str).to_numpy();barcodes=obj.obs[mapping['barcodeColumn']].astype(str).to_numpy() if mapping.get('barcodeColumn') else obj.obs_names.to_numpy()
 specs={s['id']:s for s in mapping['samples']};selected=np.array([i for i,s in enumerate(original) if specs[s]['condition']==control])
 raw=obj.X.tocsc().astype(np.int64);full=np.asarray(raw[selected].sum(axis=0)).ravel()
 annotation=json.loads((d/'annotation-plan.json').read_text());assert bytes(annotation['source']['bytes']).hex()==SOURCES[study]
 signatures=[]
 for seed in range(1,11):
  s=d/str(seed);members=json.loads((s/'memberships.json').read_text());ref=np.load(s/'reference.npz');edit=annotation['edits'][seed-1]
  assert members['sourceRows']==selected.tolist() and np.array_equal(ref['sourceRows'],selected)
  assert members['originalBarcodes']==barcodes.tolist() and members['originalSamples']==original.tolist()
  expected={};all_rows=[]
  for sid in sorted(set(original[selected])):
   positions=np.flatnonzero(original==sid)
   ranked=sorted((hashlib.sha256('|'.join(['numivivo-null-v1',study,str(seed),sid,str(barcodes[i])]).encode('utf-8')).hexdigest(),int(i)) for i in positions)
   for arm,rank in [('shamA',0),('shamB',1)]:
    rows=sorted(i for _,i in ranked[rank::2]);target=sid+'__'+arm
    expected[target]=rows;assert members['groups'][target]['sourceRows']==rows
    assert members['groups'][target]['sample']==dict(specs[sid],id=target,condition=arm)
    all_rows.extend(rows)
  assert sorted(all_rows)==selected.tolist() and len(all_rows)==len(set(all_rows))
  values=edit['value']['categorical'];labels=[values['categories'][i] for i in values['codes']]
  truth=original.copy().astype(object)
  for sid,rows in expected.items():truth[rows]=sid
  assert labels==truth.tolist()
  observed=np.stack([np.asarray(raw[expected[sid]].sum(axis=0)).ravel() for sid in members['observationSampleIDs']])
  assert np.array_equal(observed,ref['counts']) and np.array_equal(observed.sum(axis=0),full)
  signatures.append(hashlib.sha256(json.dumps(labels).encode()).hexdigest())
  checks.append(dict(study=study,seed=seed,cells=len(selected),genes=obj.n_vars,umis=int(observed.sum()),exactMembershipAndCounts=True))
 assert len(set(signatures))==10
(a.out).write_text(json.dumps(dict(status='passed-independent-source-split-reconstruction',protocolSHA256=PROTOCOL,checks=checks),indent=2)+'\n')
print(json.dumps(dict(status='passed',splits=len(checks))))
