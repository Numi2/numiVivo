"""Run native models and freeze every prediction before opening query labels."""
import argparse, datetime, hashlib, json, os, subprocess, sys, time
from pathlib import Path
import numpy as np
from prepare import sha, write
ROOT = None
HDF5 = None
FINISH = False

def main():
    binary=ROOT/'numivivo-omics';inputs=ROOT/'inputs';execution=[]
    finish=FINISH
    assert not (ROOT/'prediction-freeze.json').exists(), 'predictions already frozen; do not rerun after reveal'
    freeze=json.loads((ROOT/'input-freeze.json').read_text())
    assert all(sha(inputs/name)==digest for name,digest in freeze['files'].items())
    assert sha(binary)==json.loads((ROOT/'protocol-freeze.json').read_text())['implementationSHA256']
    env=dict(os.environ,NUMIVIVO_HDF5_LIBRARY=str(HDF5))
    def run(name,args):
        if finish: return
        start=time.monotonic();r=subprocess.run(['/usr/bin/time','-lp',str(binary),*map(str,args)],capture_output=True,text=True,env=env)
        (ROOT/(name+'.stdout')).write_text(r.stdout);(ROOT/(name+'.stderr')).write_text(r.stderr)
        execution.append(dict(name=name,args=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-start))
        write(ROOT/'execution.json',execution)
        if r.returncode: raise RuntimeError(name+' failed; retained logs')
    for method in ['logistic','knn']:
        reference=ROOT/('reference-'+method)
        if not reference.exists():run('fit-'+method,['singlecell-reference-fit',inputs/'training.h5ad','--plan',inputs/('fit-'+method+'.json'),'--output',reference])
        run('verify-reference-'+method,['singlecell-reference-verify',reference])
        mapped=ROOT/('mapped-'+method)
        run('map-'+method,['singlecell-reference-map',inputs/'query.h5ad','--plan',inputs/'query.json','--reference',reference,'--output',mapped])
        run('verify-map-'+method,['singlecell-reference-map-verify',mapped])
    run('repeat-logistic',['singlecell-reference-map',inputs/'query.h5ad','--plan',inputs/'query.json','--reference',ROOT/'reference-logistic','--output',ROOT/'mapped-logistic-repeat'])
    model=json.loads((ROOT/'reference-logistic/model.json').read_text())
    report=json.loads((ROOT/'mapped-logistic/report.json').read_text())
    repeat=json.loads((ROOT/'mapped-logistic-repeat/report.json').read_text())
    assert report==repeat
    knn=json.loads((ROOT/'mapped-knn/report.json').read_text())
    train=np.array(model['reduction']['scores']);labels=np.array(model['referenceLabels']);classes=model['classes']
    centroids=np.array([train[labels==c].mean(axis=0) for c in classes]);majority=sorted(classes,key=lambda c:(-np.count_nonzero(labels==c),c))[0]
    # Field names are supplied by the qualified owner; no test labels are needed.
    rows=report['cells'];query=np.array([r['scores'] for r in rows])
    distances=((query[:,None,:]-centroids[None,:,:])**2).sum(axis=2)
    centroid=[classes[i] for i in distances.argmin(axis=1)]
    predictions=dict(cells=[r['cell'] for r in rows],nativeLogistic=[r.get('candidateLabel') or '__unmapped__' for r in rows],nativeKNN15=[r.get('candidateLabel') or '__unmapped__' for r in knn['cells']],nearestCentroid=centroid,trainingMajority=[majority]*len(rows),classes=classes)
    write(ROOT/'predictions.json',predictions)
    write(ROOT/'prediction-freeze.json',dict(createdUTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),inputFreezeSHA256=sha(ROOT/'input-freeze.json'),runnerSHA256=sha(__file__),hdf5SHA256=sha(str(HDF5)),queryLabelsRead=False,files={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'predictions.json',ROOT/'reference-logistic/model.json',ROOT/'reference-knn/model.json',ROOT/'mapped-logistic/report.json',ROOT/'mapped-knn/report.json',ROOT/'mapped-logistic-repeat/report.json']},nativeRepeatIdentical=True))
    print('All prediction bytes frozen; query labels not opened.')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True)
    p.add_argument('--hdf5',type=Path,required=True);p.add_argument('--finish-existing-native-run',action='store_true')
    a=p.parse_args();ROOT=a.study.resolve();HDF5=a.hdf5.resolve();FINISH=a.finish_existing_native_run;main()
