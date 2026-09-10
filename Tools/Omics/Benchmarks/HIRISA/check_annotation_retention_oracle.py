#!/usr/bin/env python3
"""Verify every annotation fit and query confusion using independent per-cell SVD."""
import argparse
from pathlib import Path
import numpy as np
from check_integration_response import blocks,sha
from prepare_annotation_retention import read,write


def direct_fit(x,y):
    classes,counts=np.unique(y,return_counts=True)
    if not len(classes):return None
    encoded=np.searchsorted(classes,y)
    weights=1/(len(classes)*counts[encoded])
    center=np.sum(x*weights[:,None],axis=0)
    centered=x-center
    scale=np.sqrt(np.sum(centered**2*weights[:,None],axis=0));scale[scale<1e-12]=1
    design=np.column_stack([np.ones(len(x)),centered/scale])
    target=np.zeros((len(x),len(classes)));target[np.arange(len(x)),encoded]=1
    penalty=np.zeros((x.shape[1],x.shape[1]+1));penalty[:,1:]=np.eye(x.shape[1])
    a=np.vstack([design*np.sqrt(weights)[:,None],penalty])
    b=np.vstack([target*np.sqrt(weights)[:,None],np.zeros((x.shape[1],len(classes)))])
    beta=np.linalg.lstsq(a,b,rcond=None)[0]
    coef=beta[1:]/scale[:,None];intercept=beta[0]-center@coef
    return classes,coef,intercept


def direct_confusion(x,y,classes,coef,intercept,total):
    confusion=np.zeros((total,total))
    for first in range(0,len(x),4096):
        values=x[first:first+4096]@coef+intercept
        chosen=values==np.max(values,axis=1)[:,None]
        rows,cols=np.nonzero(chosen)
        np.add.at(confusion,(y[first:first+4096][rows],classes[cols]),1/np.sum(chosen,axis=1)[rows])
    return confusion


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('study','root','result','out'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();r=a.root;freeze=read(r/'freeze.json');result=read(a.result)
    assert not result['erasure'] and result['freezeSHA256']==sha(r/'freeze.json')
    for name,h in freeze['files'].items():assert sha(r/name)==h
    path=a.study/freeze['matrices'][result['matrix']]['path'];assert sha(path)==result['matrixSHA256']==freeze['matrices'][result['matrix']]['SHA256']
    data=np.load(r/'rows.npz');sc,dc,lc=[data[k] for k in ('stratum','donor','label')]
    n,d=freeze['cells'],freeze['components'];all_x=np.empty((n,d))
    for first,last,x in blocks(path,n,d,8192):all_x[first:last]=x
    errors=[];queried=0
    for f in result['folds']:
        training=(sc==f['stratum'])&(dc!=f['donor']);query=(sc==f['stratum'])&(dc==f['donor'])
        C=data['counts'].shape[2]
        assert np.bincount(lc[training],minlength=C).tolist()==f['trainingCounts']
        assert np.bincount(lc[query],minlength=C).tolist()==f['queryCounts']
        model=direct_fit(all_x[training],lc[training]);queried+=int(query.sum())
        if model is None:
            assert f['model'] is None;continue
        classes,coef,intercept=model
        assert classes.tolist()==f['model']['classes']
        expected_coef=np.asarray(f['model']['weight']);expected_intercept=np.asarray(f['model']['intercept'])
        assert np.allclose(coef,expected_coef,rtol=1e-8,atol=1e-9),f['id']
        assert np.allclose(intercept,expected_intercept,rtol=1e-8,atol=1e-9),f['id']
        confusion=direct_confusion(all_x[query],lc[query],classes,coef,intercept,C)
        delta=float(np.max(np.abs(confusion-np.asarray(f['confusion']))))
        assert delta<1e-7,(f['id'],'confusion',delta)
        errors.append(dict(id=f['id'],coefficientMaxAbsoluteError=float(np.max(np.abs(coef-expected_coef))),interceptMaxAbsoluteError=float(np.max(np.abs(intercept-expected_intercept))),confusionMaxAbsoluteError=delta))
    assert queried==n and len(result['folds'])==freeze['folds']
    write(a.out,dict(status='passed',queryCells=queried,folds=len(result['folds']),checks=errors,resultSHA256=sha(a.result),freezeSHA256=sha(r/'freeze.json'),oracleSHA256=sha(Path(__file__)),scope='Every original query cell and every donor fit; direct class-weighted per-cell augmented SVD, independent confusion accumulation; no new biological qualification'))
    print(dict(status='passed',cells=queried,folds=len(result['folds']),maxCoefficientError=max(x['coefficientMaxAbsoluteError'] for x in errors),maxConfusionError=max(x['confusionMaxAbsoluteError'] for x in errors)))
if __name__=='__main__':main()
