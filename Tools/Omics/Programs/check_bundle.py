#!/usr/bin/env python3
"""Prepare independent count fixtures or qualify native program bundles and replay."""
import argparse
import copy
import hashlib
import json
import subprocess
import time
from pathlib import Path
import numpy as np


def write(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(root):
    import anndata as ad
    import pandas as pd
    from scipy import sparse
    root.mkdir(exist_ok=False, parents=True)
    x=np.array([[2,1,7],[0,0,5],[0,0,0],[1,3,0]], dtype=np.uint32)
    mapping=dict(schemaVersion=1,id='program-bundle-control',evidence='synthetic',countUnit='umiCount',matrixPath='X',
                 sourceDescription='Numerical fixture only',sampleColumn='sample',featureNameColumn='name',
                 samples=[dict(id='s',biologicalReplicateID='d',condition='c',batchID='b',organism='fixture')])
    definition=dict(id='signed',organism='fixture',featureNamespace='fixture',sourceURI='urn:numivivo:numerical-fixture',
                    sourceVersion='1',sourceDescription='Numerical fixture only',members=[dict(featureID='A',weight=2),dict(featureID='B',weight=-1),dict(featureID='missing',weight=1)],minimumWeightCoverage=.75)
    plan=dict(schemaVersion=1,mapping=mapping,programs=dict(definitions=[definition]),normalizationTarget=1234)
    write(root/'plan.json',plan)
    for encoding in ['csr','csc','dense','names','ambiguous']:
        names=['A','B','C'] if encoding!='ambiguous' else ['A','A','C']
        ids=['e1','e2','e3'] if encoding in ['names','ambiguous'] else ['A','B','C']
        matrix=sparse.csc_matrix(x) if encoding=='csc' else x if encoding=='dense' else sparse.csr_matrix(x)
        ad.AnnData(matrix,obs=pd.DataFrame({'sample':['s']*4},index=['c0','c1','c2','c3']),var=pd.DataFrame({'name':names},index=ids)).write_h5ad(root/(encoding+'.h5ad'))
    total=x.sum(axis=1);normalized=np.log1p(x.astype(float)/np.where(total>0,total,1)[:,None]*1234)
    write(root/'expected.json',dict(scores=((2*normalized[:,0]-normalized[:,1])/3).tolist(),detected=(x[:,:2]>0).sum(axis=1).tolist(),totals=total.tolist()))
    write(root/'preparation.json',dict(status='prepared',numpy=np.__version__,anndata=ad.__version__,files={p.name:sha(p) for p in root.iterdir() if p.is_file()}))


def run(binary, inputs, out):
    out.mkdir(exist_ok=False, parents=True)
    commands=[]
    def command(*args, failure=False):
        began=time.monotonic();p=subprocess.run([str(binary.resolve()),*map(str,args)],capture_output=True,text=True)
        i=len(commands);(out/f'{i:02d}.stdout').write_text(p.stdout);(out/f'{i:02d}.stderr').write_text(p.stderr)
        commands.append(dict(arguments=list(map(str,args)),exitCode=p.returncode,expectedFailure=failure,seconds=time.monotonic()-began))
        write(out/'commands.json',commands)
        assert p.returncode==(65 if failure else 0),p.stderr[-4000:]
    expected=json.loads((inputs/'expected.json').read_text())
    dtype=np.dtype([('row','<u4'),('column','<u4'),('bits','<u8')])
    plan=json.loads((inputs/'plan.json').read_text());models=[]
    for encoding in ['csr','csc','dense','names']:
        current=copy.deepcopy(plan);current['matchFeatureNames']=encoding=='names';path=out/(encoding+'-plan.json');write(path,current)
        bundle=out/encoding
        command('singlecell-h5ad-programs',inputs/(encoding+'.h5ad'),'--plan',path,'--output',bundle)
        command('singlecell-h5ad-programs-verify',bundle)
        score=np.fromfile(bundle/'scores.bin',dtype=dtype);detected=np.fromfile(bundle/'detected-members.bin',dtype=dtype);totals=np.fromfile(bundle/'total-counts.bin',dtype=dtype)
        for records in [score,detected,totals]:
            np.testing.assert_array_equal(records['row'],np.arange(4));np.testing.assert_array_equal(records['column'],np.zeros(4,dtype=np.uint32))
        np.testing.assert_allclose(score['bits'].view('<f8'),expected['scores'],rtol=1e-14,atol=1e-14)
        np.testing.assert_array_equal(detected['bits'],expected['detected']);np.testing.assert_array_equal(totals['bits'],expected['totals'])
        model=json.loads((bundle/'model.json').read_text());models.append(model)
        assert model['sourcePasses']==2 and model['programs'][0]['missingFeatureIDs']==['missing']
        assert model['programs'][0]['weightCoverage']==.75
        assert [c['barcode'] for c in json.loads((bundle/'metadata.json').read_text())['cells']]==['c0','c1','c2','c3']
        if encoding=='names':assert json.loads((bundle/'metadata.json').read_text())['features'][0]['id']=='e1'
    assert all(m==models[0] for m in models)
    command('singlecell-h5ad-programs',inputs/'csr.h5ad','--plan',out/'csr-plan.json','--output',out/'repeat')
    assert (out/'repeat/receipt.json').read_bytes()==(out/'csr/receipt.json').read_bytes()
    command('singlecell-h5ad-programs',inputs/'csr.h5ad','--plan',out/'csr-plan.json','--output',out/'csr',failure=True)
    # The existing JSON workflow still shares exactly the same scoring arithmetic.
    legacy=dict(schemaVersion=1,mapping=plan['mapping'],programs=plan['programs'],reduction={'normalizationTarget':1234})
    # Avoid PCA on the tiny fixture: default normalization is 10000 instead.
    legacy.pop('reduction');write(out/'legacy-plan.json',legacy)
    command('singlecell-h5ad-pseudobulk',inputs/'csr.h5ad','--plan',out/'legacy-plan.json','--output',out/'legacy')
    ordinary=copy.deepcopy(plan);ordinary['normalizationTarget']=10000;write(out/'ordinary-plan.json',ordinary)
    command('singlecell-h5ad-programs',inputs/'csr.h5ad','--plan',out/'ordinary-plan.json','--output',out/'ordinary')
    legacy_result=json.loads((out/'legacy/report.json').read_text())['programs']
    actual=np.fromfile(out/'ordinary/scores.bin',dtype=dtype)['bits'].view('<f8')
    assert actual.tolist()==[0 if row[0] is None else row[0] for row in legacy_result['scores']]
    for name in ['scores.bin','detected-members.bin','total-counts.bin','metadata.json','model.json']:
        path=out/'csr'/name;original=path.read_bytes();path.write_bytes(original[:-1])
        command('singlecell-h5ad-programs-verify',out/'csr',failure=True);path.write_bytes(original)
    # A rehashed fabricated result must fail reconstruction, not just a hash check.
    path=out/'csr/scores.bin';original=path.read_bytes();fake=bytearray(original);fake[8]^=1;path.write_bytes(fake)
    receipt_path=out/'csr/receipt.json';receipt_bytes=receipt_path.read_bytes();receipt=json.loads(receipt_bytes)
    # VivoFingerprint retains the exact 32-byte array.
    receipt['scores']={'bytes':list(bytes.fromhex(sha(path)))};receipt_path.write_text(json.dumps(receipt,sort_keys=True,separators=(',',':')))
    command('singlecell-h5ad-programs-verify',out/'csr',failure=True)
    assert 'reconstruction differs' in (out/f'{len(commands)-1:02d}.stderr').read_text()
    path.write_bytes(original);receipt_path.write_bytes(receipt_bytes)
    for case in ['ambiguous','score-budget','update-budget','unknown','normalization','coverage']:
        bad=copy.deepcopy(plan);source=inputs/'csr.h5ad'
        if case=='ambiguous':bad['matchFeatureNames']=True;source=inputs/'ambiguous.h5ad'
        elif case=='score-budget':bad['programs']['maximumScoreValues']=1
        elif case=='update-budget':bad['programs']['maximumUpdates']=1
        elif case=='unknown':bad['cellSelection']=[0]
        elif case=='normalization':bad['normalizationTarget']=0
        elif case=='coverage':bad['programs']['definitions'][0]['minimumWeightCoverage']=1
        path=out/(case+'-plan.json');write(path,bad);target=out/('rejected-'+case)
        command('singlecell-h5ad-programs',source,'--plan',path,'--output',target,failure=True)
        assert not target.exists()
    assert not list(out.glob('.numivivo-programs-*'))
    command('singlecell-h5ad-programs-verify',out/'csr')
    write(out/'checks.json',dict(status='passed',commands=len(commands),expectedRejections=sum(c['expectedFailure'] for c in commands),binarySHA256=sha(binary),checkerSHA256=sha(Path(__file__)),scope='CSR/CSC/dense, exact names and original IDs, signed scores and empty-cell mask, original JSON arithmetic, repeat/replay, corruption and rehashed-result reconstruction rejection; controlled numerical fixtures only.'))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepare',type=Path);p.add_argument('--binary',type=Path);p.add_argument('--inputs',type=Path);p.add_argument('--out',type=Path);a=p.parse_args()
    if a.prepare:prepare(a.prepare)
    else:
        assert a.binary and a.inputs and a.out
        run(a.binary,a.inputs,a.out)


if __name__=='__main__':main()
