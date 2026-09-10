#!/usr/bin/env python3
"""Measure complete-cohort donor-held-out retention of supplied annotation labels."""
import argparse
from pathlib import Path
import numpy as np
from check_integration_response import blocks,sha
from prepare_annotation_retention import read,write


def fit(count, sx, sxx):
    classes=np.flatnonzero(count>0)
    if not len(classes):return None
    means=sx[classes]/count[classes,None]
    center=means.mean(axis=0);second=(sxx[classes]/count[classes,None,None]).mean(axis=0)
    var=np.diag(second)-center**2
    assert var.min()>=-1e-10*max(1,float(np.abs(second).max()))
    scale=np.sqrt(np.maximum(var,0));scale[scale<1e-12]=1
    gram=(second-np.outer(center,center))/np.outer(scale,scale)
    rhs=(means-center).T/scale[:,None]/len(classes)
    weight=np.linalg.solve(gram+np.eye(len(center)),rhs)/scale[:,None]
    intercept=np.full(len(classes),1/len(classes))-center@weight
    assert np.isfinite(weight).all() and np.isfinite(intercept).all()
    return dict(classes=classes.tolist(),center=center.tolist(),scale=scale.tolist(),weight=weight.tolist(),intercept=intercept.tolist())


def training_stats(count,sx,sxx,stratum,held):
    keep=np.arange(count.shape[1])!=held
    return count[stratum,keep].sum(axis=0),sx[stratum,keep].sum(axis=0),sxx[stratum,keep].sum(axis=0)


def accumulate(confusion,x,labels,model):
    if model is None:return
    values=x@np.asarray(model['weight'])+model['intercept']
    assert np.isfinite(values).all()
    ties=values==values.max(axis=1,keepdims=True)
    fractions=ties/ties.sum(axis=1,keepdims=True)
    for col,label in enumerate(model['classes']):
        confusion[:,label]+=np.bincount(labels,weights=fractions[:,col],minlength=len(confusion))


def measures(confusion,truth):
    output=[];pred=confusion.sum(axis=0)
    for label,n in enumerate(truth):
        tp=confusion[label,label]
        recall=float(tp/n) if n else None
        precision=float(tp/pred[label]) if pred[label] else None
        f1=float(2*tp/(n+pred[label])) if n+pred[label] else None
        output.append(dict(label=label,queryCells=int(n),recall=recall,precision=precision,f1=f1))
    return output


def comparisons(baseline,result,ledger):
    assert len(baseline['folds'])==len(result['folds'])==len(ledger['folds'])
    groups={}
    for a,b in zip(baseline['folds'],result['folds']):
        assert (a['id'],a['stratum'],a['donor'],a['trainingCounts'],a['queryCounts'])==(b['id'],b['stratum'],b['donor'],b['trainingCounts'],b['queryCounts'])
        supported=sum(n>0 for n in b['trainingCounts'])
        for first,second in zip(a['metrics'],b['metrics']):
            if not first['queryCells']:continue
            label=first['label'];control=1/supported if b['trainingCounts'][label]>0 else 0.
            groups.setdefault((b['stratum'],label),[]).append(dict(donor=b['donor'],baselineRecall=first['recall'],candidateRecall=second['recall'],erasureRecall=control,queryCells=first['queryCells'],trainingCells=b['trainingCounts'][label],sufficient=b['trainingCounts'][label]>=20 and first['queryCells']>=20 and a['model'] is not None and b['model'] is not None))
    out=[]
    for (stratum,label),folds in sorted(groups.items()):
        sufficient=all(f['sufficient'] for f in folds)
        base=float(np.mean([f['baselineRecall'] for f in folds]));control=float(np.mean([f['erasureRecall'] for f in folds]))
        losses=[f['baselineRecall']-f['candidateRecall'] for f in folds]
        sensitive=base>=control+.05
        out.append(dict(stratum=stratum,label=label,labelName=ledger['labels'][label],rare=ledger['rare'][stratum][label],folds=folds,baselineMeanRecall=base,candidateMeanRecall=float(np.mean([f['candidateRecall'] for f in folds])),erasureMeanRecall=control,meanRecallLoss=float(np.mean(losses)),maximumFoldRecallLoss=max(losses),sufficientSupport=sufficient,controlSensitive=sensitive,qualified=bool(sufficient and sensitive and np.mean(losses)<=.05 and max(losses)<=.10),marginFailure=bool(sufficient and sensitive and (np.mean(losses)>.05 or max(losses)>.10))))
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('study','root','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--matrix',required=True);p.add_argument('--baseline',type=Path);p.add_argument('--erase',action='store_true')
    a=p.parse_args();r=a.root;freeze=read(r/'freeze.json');ledger=read(r/'ledger.json')
    for name,h in freeze['files'].items():assert sha(r/name)==h
    path=a.study/freeze['matrices'][a.matrix]['path'];assert sha(path)==freeze['matrices'][a.matrix]['SHA256']
    arrays=np.load(r/'rows.npz');sc,dc,lc=[arrays[k] for k in ('stratum','donor','label')];count=arrays['counts']
    S,D,C=count.shape;n=freeze['cells'];d=freeze['components']
    assert n==len(lc)==int(count.sum()) and count.tolist()==ledger['counts']
    sx=np.zeros((S,D,C,d));sxx=np.zeros((S,D,C,d,d))
    code=(sc.astype(np.int64)*D+dc)*C+lc
    for first,last,x in blocks(path,n,d,8192):
        if a.erase:x.fill(0)
        cc=code[first:last]
        for k in np.unique(cc):
            y=x[cc==k];sx.reshape(-1,d)[k]+=y.sum(axis=0);sxx.reshape(-1,d,d)[k]+=y.T@y
    folds=[];by_pair={}
    for f in ledger['folds']:
        s,donor=f['stratum'],f['donor'];train,train_x,train_xx=training_stats(count,sx,sxx,s,donor)
        model=fit(train,train_x,train_xx)
        record=dict(**f,trainingCounts=train.tolist(),queryCounts=count[s,donor].tolist(),model=model,status='measured' if model else 'unavailable-training')
        record['confusion']=np.zeros((C,C));by_pair[s,donor]=record;folds.append(record)
    for first,last,x in blocks(path,n,d,8192):
        if a.erase:x.fill(0)
        ss,dd,ll=sc[first:last],dc[first:last],lc[first:last]
        for code in np.unique(ss.astype(np.int64)*D+dd):
            s,donor=divmod(int(code),D);mask=(ss==s)&(dd==donor);f=by_pair[s,donor]
            accumulate(f['confusion'],x[mask],ll[mask],f['model'])
    for f in folds:
        confusion=f['confusion'];truth=np.asarray(f['queryCounts'])
        assert np.allclose(confusion.sum(axis=1),truth if f['model'] else np.zeros(C),rtol=0,atol=1e-7)
        f['metrics']=measures(confusion,truth);f['confusion']=confusion.tolist()
        f['accuracy']=float(np.trace(confusion)/truth.sum());f['balancedRecall']=float(np.mean([m['recall'] for m in f['metrics'] if m['recall'] is not None]))
        if a.erase and f['model']:
            expected=np.zeros((C,C));expected[:,f['model']['classes']]=truth[:,None]/len(f['model']['classes'])
            assert np.allclose(confusion,expected,rtol=1e-10,atol=1e-8)
    result=dict(schemaVersion=1,cells=n,folds=folds,foldCount=len(folds),matrix=a.matrix,erasure=a.erase,freezeSHA256=sha(r/'freeze.json'),matrixSHA256=sha(path),scorerSHA256=sha(Path(__file__)),labelsAuthoritative=False,scope='Full-cohort transductive preservation of supplied author l2 annotations; no prospective or independent biological qualification')
    result['comparisons']=comparisons(read(a.baseline) if a.baseline else result,result,ledger)
    result['completePreservationQualified']=all(x['qualified'] for x in result['comparisons'])
    write(a.out,result)
    print(dict(matrix=a.matrix,erasure=a.erase,folds=len(folds),comparisons=len(result['comparisons']),rare=sum(x['rare'] for x in result['comparisons']),sensitiveSufficient=sum(x['sufficientSupport'] and x['controlSensitive'] for x in result['comparisons']),marginFailures=sum(x['marginFailure'] for x in result['comparisons']),complete=result['completePreservationQualified']))

if __name__=='__main__':main()
