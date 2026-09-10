#!/usr/bin/env python3
"""Freeze all original annotation strata and donor folds without reading RNA values."""
import argparse,json,time
from pathlib import Path
import h5py
import numpy as np
from check_integration_response import sha,source_codes,one

def read(p):return json.loads(p.read_text())
def write(p,v):
 with p.open('x') as f:json.dump(v,f,sort_keys=True,indent=2,allow_nan=False);f.write('\n')
def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','repo','out'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();s=a.study;r=a.out
 assert not (r/'freeze.json').exists()
 source=s/'native-release/original.h5ad';metadata=s/'integration-full/native/metadata.json'
 assert sha(source)=='0873e698ebf8770a54e6dba09724ffbeda5e1a67dbf24d4e223552d4aca7969c'
 assert sha(metadata)=='f81771d3a3f2db3fc23a4ec21b6d52f160ed2b56f843e567dfdf443e8cfc2b7f'
 design=read(s/'design.json');samples=design['samples'];n=1612594
 lib=source_codes(source,metadata,samples,n,8192)
 strata=sorted({(one(x,'cell type'),one(x,'treatment')) for x in samples});donors=sorted({one(x,'subject id') for x in samples})
 sc=np.array([strata.index((one(x,'cell type'),one(x,'treatment'))) for x in samples],dtype=np.uint16)[lib]
 dc=np.array([donors.index(one(x,'subject id')) for x in samples],dtype=np.uint16)[lib]
 provisional={};labels=np.empty(n,dtype=np.uint16)
 with h5py.File(source) as h:
  for first in range(0,n,8192):
   values=h['obs/celltype.l2'].asstr()[first:first+8192]
   for j,v in enumerate(values):
    assert v and len(v)<=256
    if v not in provisional:provisional[v]=len(provisional)
    labels[first+j]=provisional[v]
 names=sorted(provisional);remap=np.array([names.index(k) for k in provisional],dtype=np.uint16);labels=remap[labels]
 counts=np.zeros((len(strata),len(donors),len(names)),dtype=np.int64)
 np.add.at(counts,(sc,dc,labels),1);assert counts.sum()==n
 np.savez_compressed(r/'rows.npz',library=lib,stratum=sc,donor=dc,label=labels,counts=counts)
 folds=[dict(id=f's{si:02d}-{donor}',stratum=si,donor=di,cells=int(counts[si,di].sum())) for si in range(len(strata)) for di,donor in enumerate(donors) if counts[si,di].sum()]
 assert sum(x['cells'] for x in folds)==n
 ledger=dict(cells=n,strata=[dict(preparation=x,treatment=y) for x,y in strata],donors=donors,labels=names,libraries=[x['accession'] for x in samples],folds=folds,counts=counts.tolist(),rare=((counts.sum(axis=1)/counts.sum(axis=(1,2))[:,None])<.01).tolist(),labelField='celltype.l2',labelsAuthoritative=False)
 write(r/'ledger.json',ledger)
 paths={'baseline':s/'integration-full/pca/scores.bin','native':s/'integration-full/native/scores.bin',**{'harmony-'+str(k):s/'integration-harmony'/('seed-'+str(k))/'scores.bin' for k in (7,19,41)}}
 expected={'baseline':read(s/'integration-harmony/seed-7/report.json')['bindings']['scores']['SHA256'],'native':bytes(read(s/'integration-full/native/receipt.json')['scores']['bytes']).hex(),**{'harmony-'+str(k):read(s/'integration-harmony'/('seed-'+str(k))/'report.json')['output']['SHA256'] for k in (7,19,41)}}
 for k,v in paths.items():assert sha(v)==expected[k] and v.stat().st_size==n*20*16
 write(r/'freeze.json',dict(schemaVersion=1,declaredAtUnix=time.time(),cells=n,components=20,folds=len(folds),sourceSHA256=sha(source),metadataSHA256=sha(metadata),designSHA256=sha(s/'design.json'),files={k:sha(r/k) for k in ('protocol.md','rows.npz','ledger.json')},matrices={k:dict(path=str(v.relative_to(s)),SHA256=expected[k]) for k,v in paths.items()},preparerSHA256=sha(Path(__file__)),metricsRead=False))
 print(dict(cells=n,strata=len(strata),labels=len(names),donors=len(donors),folds=len(folds),rareClassStrata=int(np.count_nonzero((counts.sum(axis=1)>0)&np.asarray(ledger['rare'])))))
if __name__=='__main__':main()
