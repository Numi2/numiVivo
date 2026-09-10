#!/usr/bin/env python3
"""Exercise the actual target-kernel CLI with numerical controls and rejected artifacts."""
import argparse
import copy
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
import numpy as np
from combinations import sha, write
from go_transfer import independent_prediction
from threadpoolctl import threadpool_limits


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    context=dict(id='synthetic-context',organism='NCBITaxon:9606',featureNamespace='synthetic-gene',perturbationNamespace='synthetic-target',countUnit='umiCount')
    fingerprint=dict(bytes=[0]*32)
    raw=np.array([[10,20,30,0],[30,10,20,1],[0,20,10,4],[10,1,50,0],[20,30,40,1]],dtype=np.uint64)
    names=['a','b','c','d'];terms=[['x','y'],['y','z'],['w','x'],[]]
    offsets=[0];columns=[];values=[]
    for row in raw:
        ix=np.flatnonzero(row);columns.extend(ix.tolist());values.extend(row[ix].tolist());offsets.append(len(values))
    training=dict(schemaVersion=1,selection=dict(schemaVersion=1,context=context,controlCondition='control',targets=[dict(id=t,condition=t) for t in names],provenance='Synthetic numerical control'),
                  source=fingerprint,sourceReport=fingerprint,evidence='synthetic',featureIDs=['g0','g1','g2','g3'],
                  matrix=dict(cellCount=5,featureCount=4,rowOffsets=offsets,featureIndices=columns,counts=values))
    plan=dict(schemaVersion=1,context=context,descriptorNamespace='synthetic-terms',source=fingerprint,provenance='Synthetic numerical control',
              targets=[dict(targetID=t,featureID='gene-'+t,status='available' if term else 'noData',terms=term) for t,term in zip(names,terms)],regularization=1,maximumWork=200000000)
    descriptors=[dict(targetID='held',featureID='gene-held',status='available',terms=['x','y']),
                 dict(targetID='missing',featureID='gene-missing',status='noData',terms=[]),
                 dict(targetID='unresolved',status='unresolvedIdentity',terms=[]),
                 dict(targetID='disjoint',featureID='gene-disjoint',status='available',terms=['unshared'])]
    query=dict(schemaVersion=1,context=context,descriptorNamespace=plan['descriptorNamespace'],source=fingerprint,queries=[dict(id=d['targetID'],descriptor=d) for d in descriptors])
    for name,value in [('training',training),('plan',plan),('query',query)]:write(a.out/(name+'.json'),value)
    commands=[]
    def run(name,args,success=True):
        with (a.out/(name+'.log')).open('w') as log:r=subprocess.run([str(a.binary.resolve()),*map(str,args)],stdout=log,stderr=subprocess.STDOUT)
        commands.append(dict(name=name,returnCode=r.returncode,expectedSuccess=success));write(a.out/'commands.json',commands)
        assert r.returncode==(0 if success else 65),(name,r.returncode)
    model=a.out/'model';prediction=a.out/'prediction'
    fit=['singlecell-target-kernel-fit',a.out/'training.json','--plan',a.out/'plan.json','--output',model]
    predict=['singlecell-target-kernel-predict',model,'--plan',a.out/'query.json','--output',prediction]
    run('fit',fit);run('verify',['singlecell-target-kernel-verify',model]);run('predict',predict);run('prediction-verify',['singlecell-target-kernel-prediction-verify',prediction])
    report=json.loads((prediction/'report.json').read_text());assert [r['status'] for r in report['queries']]==['predicted','noData','unresolvedIdentity','noSharedTerms']
    cpm=raw/raw.sum(axis=1,keepdims=True)*1e6;baseline=np.log1p(cpm[0]);y=np.log1p(cpm[1:4])-baseline
    sets=[set(t) for t in terms[:3]];kernel=np.array([[len(x&z)/len(x|z) for z in sets] for x in sets]);similarity=np.array([len(set(['x','y'])&x)/len(set(['x','y'])|x) for x in sets])
    errors=[]
    for field,response in [('expression',y),('shuffledExpression',np.roll(y,-1,axis=0))]:
        oracle=np.maximum(baseline+independent_prediction(kernel,similarity,response,1),0);native=np.array(report['queries'][0][field])
        np.testing.assert_allclose(native,oracle,atol=1e-12,rtol=1e-12);errors.append(float(np.max(abs(native-oracle))))
    for i,expected in enumerate([baseline,np.maximum(baseline+(np.log1p(cpm[1:])-baseline).mean(axis=0),0),np.maximum(baseline+y.mean(axis=0),0)]):
        np.testing.assert_allclose(report['baselines'][i]['expression'],expected,atol=1e-12,rtol=1e-12)
    for name in ['context','namespace','seen','known-alias','unsupported-training-alias','duplicate-query','empty-query','unknown-field','duplicate-terms','work']:
        changed=copy.deepcopy(query)
        if name=='context':changed['context']['id']='other'
        elif name=='namespace':changed['descriptorNamespace']='other'
        elif name=='seen':changed['queries'][0]['descriptor']['targetID']='a'
        elif name=='known-alias':changed['queries'][0]['descriptor']['featureID']='gene-b'
        elif name=='unsupported-training-alias':changed['queries'][0]['descriptor']['featureID']='gene-d'
        elif name=='duplicate-query':changed['queries'][1]['id']='held'
        elif name=='empty-query':changed['queries']=[]
        elif name=='unknown-field':changed['heldOutRNA']=[1,2]
        elif name=='duplicate-terms':changed['queries'][0]['descriptor']['terms']=['x','x']
        else:changed['maximumWork']=1
        path=a.out/(name+'.json');write(path,changed);destination=a.out/('reject-'+name)
        run('reject-'+name,['singlecell-target-kernel-predict',model,'--plan',path,'--output',destination],False);assert not destination.exists()
    for name in ['missing-training-descriptor','extra-training-descriptor','duplicate-feature','invalid-regularization','fit-work']:
        changed=copy.deepcopy(plan)
        if name=='missing-training-descriptor':changed['targets'].pop()
        elif name=='extra-training-descriptor':changed['targets'].append(dict(targetID='extra',featureID='gene-extra',status='available',terms=['x']))
        elif name=='duplicate-feature':changed['targets'][1]['featureID']='gene-a'
        elif name=='invalid-regularization':changed['regularization']=-1
        else:changed['maximumWork']=1
        path=a.out/(name+'.json');write(path,changed);destination=a.out/('reject-'+name)
        run('reject-'+name,['singlecell-target-kernel-fit',a.out/'training.json','--plan',path,'--output',destination],False);assert not destination.exists()
    run('reject-model-overwrite',fit,False);run('reject-prediction-overwrite',predict,False)
    for name in ['model','plan','query','report','implementation']:
        is_prediction=name in ['query','report'];copy_path=a.out/('tampered-'+name);shutil.copytree(prediction if is_prediction else model,copy_path)
        receipt_path=copy_path/'receipt.json';receipt=json.loads(receipt_path.read_text())
        if name=='implementation':receipt['implementation']['bytes']=[1]*32
        else:
            path=copy_path/(name+'.json');object=json.loads(path.read_text())
            if name=='model':object['inverse'][0]+=0.01
            elif name=='plan':object['targets'][0]['terms']=['w','z']
            elif name=='query':object['queries'][0]['id']='altered'
            else:object['queries'][0]['expression'][0]+=0.1
            # Swift receipts contain integer-only hashes; compact sorted encoding is canonical.
            path.write_text(json.dumps(object,sort_keys=True,separators=(',',':')))
            receipt['result' if name=='report' else name]=dict(bytes=list(hashlib.sha256(path.read_bytes()).digest()))
        receipt_path.write_text(json.dumps(receipt,sort_keys=True,separators=(',',':')))
        run('reject-rehashed-'+name,['singlecell-target-kernel-prediction-verify' if is_prediction else 'singlecell-target-kernel-verify',copy_path],False)
    assert not list(a.out.rglob('.numivivo-target-kernel-*'))
    write(a.out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(not r['expectedSuccess'] for r in commands),
          independentMaximumPredictionError=max(errors),binarySHA256=sha(a.binary),implementationSHA256=sha(Path(__file__)),scope='Synthetic numerical and artifact-integrity checks, not biological qualification.'))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()
