#!/usr/bin/env python3
"""Exercise frozen PCA query lifecycle, sparse formats, feature alignment and rejection controls."""
import argparse,copy,hashlib,json,shutil,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
from pca_query_reference import check
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
def write(path,value):path.write_text(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False))
def run(args,error=None):
    r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True)
    (a.out/(str(len(commands))+'.log')).write_text(r.stdout+r.stderr)
    commands.append(dict(arguments=list(map(str,args)),exitCode=r.returncode,expectedRejection=error));write(a.out/'commands.json',commands)
    assert (r.returncode!=0 and error in r.stderr) if error else r.returncode==0,r.stderr
old=json.loads((a.fixtures/'plan.json').read_text());fit={k:old[k] for k in ['schemaVersion','mapping','reduction']};fit['featureNamespace']='fixture-gene-ID'
write(a.out/'fit.json',fit)
reference=a.out/'reference';run(['singlecell-h5ad-pca',a.fixtures/'csr.h5ad','--plan',a.out/'fit.json','--output',reference])
qp=dict(schemaVersion=1,mapping=copy.deepcopy(fit['mapping']),featureNamespace=fit['featureNamespace']);qp['mapping'].pop('groupColumn',None);write(a.out/'query.json',qp)
source=ad.read_h5ad(a.fixtures/'csr.h5ad');source.obs_names=['query-'+str(v) for v in source.obs_names]
checks=[];base=None
for kind in ['csr','csc','reordered']:
    data=source[:,::-1].copy() if kind=='reordered' else source.copy()
    data.X=data.X.tocsc() if kind=='csc' else data.X.tocsr()
    path=a.out/(kind+'.h5ad');data.write_h5ad(path)
    out=a.out/kind;run(['singlecell-h5ad-pca-query',path,'--plan',a.out/'query.json','--reference',reference,'--output',out])
    run(['singlecell-h5ad-pca-query-verify',out])
    result,scores=check(out,path);checks.append(result)
    if base is None:base=scores
    else:np.testing.assert_allclose(scores,base,rtol=1e-12,atol=1e-12)
run(['singlecell-h5ad-pca-query',a.out/'csr.h5ad','--plan',a.out/'query.json','--reference',reference,'--output',a.out/'repeat'])
assert (a.out/'repeat/receipt.json').read_bytes()==(a.out/'csr/receipt.json').read_bytes()
run(['singlecell-h5ad-pca-query',a.out/'csr.h5ad','--plan',a.out/'query.json','--reference',reference,'--output',a.out/'csr'],'destination exists')
for name,edit,error in [
    ('namespace',lambda v:v.update(featureNamespace='different'),'namespace or count unit mismatch'),
    ('units',lambda v:v['mapping'].update(countUnit='readCount'),'namespace or count unit mismatch'),
    ('organism',lambda v:v['mapping']['samples'][0].update(organism='different'),'gene universe and organism'),
    ('labels',lambda v:v['mapping'].update(groupColumn='label'),'label mapping'),
    ('budget',lambda v:v.update(maximumProjectionUpdates=1),'projection-update budget')]:
    bad=copy.deepcopy(qp);edit(bad);write(a.out/(name+'.json'),bad)
    run(['singlecell-h5ad-pca-query',a.out/'csr.h5ad','--plan',a.out/(name+'.json'),'--reference',reference,'--output',a.out/('reject-'+name)],error)
    assert not (a.out/('reject-'+name)).exists()
legacy=copy.deepcopy(fit);legacy.pop('featureNamespace');write(a.out/'legacy-fit.json',legacy)
run(['singlecell-h5ad-pca',a.fixtures/'csr.h5ad','--plan',a.out/'legacy-fit.json','--output',a.out/'legacy-reference'])
run(['singlecell-h5ad-pca-query',a.out/'csr.h5ad','--plan',a.out/'query.json','--reference',a.out/'legacy-reference','--output',a.out/'unbound-namespace'],'training namespace required')
missing=source[:,:-1].copy();missing.write_h5ad(a.out/'missing.h5ad')
run(['singlecell-h5ad-pca-query',a.out/'missing.h5ad','--plan',a.out/'query.json','--reference',reference,'--output',a.out/'missing'],'gene universe and organism')
run(['singlecell-h5ad-pca-query',a.fixtures/'csr.h5ad','--plan',a.out/'query.json','--reference',reference,'--output',a.out/'overlap'],'overlaps training cell identities')
for kind in ['coordinate','value','report','reference']:
    dst=a.out/('tamper-'+kind);shutil.copytree(a.out/'csr',dst)
    if kind=='reference':
        path=dst/'reference/loadings.bin';receiptPath=dst/'reference/receipt.json';key='loadings';error='PCA source reconstruction differs'
    else:
        path=dst/('report.json' if kind=='report' else 'scores.bin');receiptPath=dst/'receipt.json';key='report' if kind=='report' else 'scores';error='PCA query source reconstruction differs'
    if kind=='report':
        value=json.loads(path.read_text());value['projectionUpdates']+=1;write(path,value)
    else:
        value=bytearray(path.read_bytes());value[0 if kind=='coordinate' else 8]^=1;path.write_bytes(value)
    receipt=json.loads(receiptPath.read_text());receipt[key]={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(receiptPath,receipt)
    run(['singlecell-h5ad-pca-query-verify',dst],error)
assert not list(a.out.glob('.numivivo-pca*'))
write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(c['expectedRejection'] is not None for c in commands),formats=checks,replayExact=True,rehashedTamperingRejected=True,binarySHA256=hashlib.sha256(a.binary.read_bytes()).hexdigest()))
