#!/usr/bin/env python3
"""Freeze HIRISA predictions, independently reconstruct fits, then score every fold."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import ijson
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def read(path): return json.loads(path.read_text())
def write(path,value): path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')
def fingerprint(path): return {'bytes':list(bytes.fromhex(sha(path)))}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root;bundle=root/'prediction-batch';frozen=read(root/'prediction-folds.json')
    transport=read(root/'prediction-transport-receipt.json');audit=read(root/'native-independent-verification.json')
    receipt=read(bundle/'receipt.json');plan=read(bundle/'plan.json')
    assert sha(root/'prediction-folds.json')==transport['frozenFoldsSHA256']
    assert sha(root/'prediction-batch-plan.json')==transport['planSHA256']
    assert plan==read(root/'prediction-batch-plan.json') and receipt['plan']==fingerprint(bundle/'plan.json')
    assert receipt['source']=={'bytes':list(bytes.fromhex(frozen['sourceH5ADSHA256']))}
    assert receipt['sourceReport']==plan['sourceReport']=={'bytes':list(bytes.fromhex(audit['reportSHA256']))}
    assert len(receipt['folds'])==len(plan['folds'])==len(frozen['folds'])==79
    assert audit['status']=='passed' and audit['nativeRawIngestionLinked']
    assert sha(root/'native-independent-verification.json')==transport['nativeVerificationSHA256']
    # Freeze every output BEFORE reading any count array, including targets.
    paths=[bundle/'plan.json',bundle/'receipt.json']
    for i,(status,fold) in enumerate(zip(receipt['folds'],frozen['folds'])):
        assert status['id']==fold['id']
        d=bundle/'folds'/f'{i:03d}';assert read(d/'receipt.json')==status
        paths.append(d/'receipt.json')
        if status['status']=='completed':
            for name in ('model','prediction'):
                path=d/(name+'.json');assert status[name]==fingerprint(path);paths.append(path)
        else:
            assert status['status']=='failed' and status.get('failure')
            assert not (d/'model.json').exists() and not (d/'prediction.json').exists()
    identities={str(path.relative_to(root)):sha(path) for path in paths}
    native_freeze=read(root/'prediction-native-output-freeze.json')
    assert native_freeze['targetsOpenedForScoring'] is False and native_freeze['folds']==79
    assert {str(path.relative_to(bundle)):identities[str(path.relative_to(root))] for path in paths}==native_freeze['files']
    freeze_path=root/'prediction-output-freeze.json'
    if freeze_path.exists(): assert read(freeze_path)['files']==identities
    else: write(freeze_path,dict(schemaVersion=1,createdAtUnix=time.time(),files=identities,foldCount=79,targetsOpenedForScoring=False))
    assert sha(root/'inference-inputs/counts.npz')==audit['inferenceCountsSHA256']
    with np.load(root/'inference-inputs/counts.npz') as archive:
        counts=archive['counts'];samples=archive['sampleIDs'].tolist();genes=archive['featureIDs'].tolist()
    assert len(genes)==18082 and len(samples)==131
    sample_rows={name:i for i,name in enumerate(samples)}
    # Query memberships are compared with the independently audited source report.
    report=root/'native-full/report.json';assert sha(report)==audit['reportSHA256']
    with report.open('rb') as f: groups=list(ijson.items(f,'pseudobulk.groups.item'))
    by_sample={g['sampleIDs'][0]:g for g in groups};results=[];methods=['noChange','meanResponse','medianResponse','contextRidge']
    for i,fold in enumerate(frozen['folds']):
        status=receipt['folds'][i];base={k:fold[k] for k in ('id','cohort','population','treatment','heldOutDonor')}
        if status['status']!='completed':
            results.append({**base,'status':'failed','failure':status['failure'],'scores':[]});continue
        d=bundle/'folds'/f'{i:03d}';model=read(d/'model.json');prediction=read(d/'prediction.json')
        assert model['method']==frozen['method'] and model['featureIDs']==prediction['featureIDs']==genes
        assert model['source']==status['trainingAggregate']==prediction['referenceSource']
        pairs=sorted(fold['trainingPairs'],key=lambda v:v['donor'])
        ids=[pair[key] for pair in pairs for key in ('control','treated')]
        native_fold=plan['folds'][i]
        assert [groups[j]['sampleIDs'][0] for j in native_fold['trainingGroupIndices']]==ids
        assert [groups[j]['sampleIDs'][0] for j in native_fold['queryGroupIndices']]==[fold['queryControlAccession']]
        assert fold['scoringTargetAccession'] not in ids+[fold['queryControlAccession']]
        assert model['trainingDonors']==[pair['donor'] for pair in pairs]
        assert len(prediction['predictions'])==1;query=prediction['predictions'][0]
        assert query['group']['sourceCellIndices']==by_sample[fold['queryControlAccession']]['sourceCellIndices']
        assert query['group']['donorID']==fold['heldOutDonor'] and query['group']['cellGroup']==fold['population']
        yy=counts[[sample_rows[s] for s in ids]].astype(float)
        logs=np.log1p(yy/yy.sum(axis=1,keepdims=True)*1e6);controls=logs[::2];response=logs[1::2]-controls
        paired=yy[::2]+yy[1::2];selected=np.flatnonzero((paired.sum(axis=0)>=10)&((paired>0).sum(axis=0)>=2))
        assert selected.tolist()==model['selectedFeatureIndices']
        centers=controls[:,selected].mean(axis=0);scales=controls[:,selected].std(axis=0)
        constant=np.all(controls[:,selected]==controls[0,selected],axis=0)
        centers[constant]=controls[0,selected[constant]];scales[constant|(scales==0)]=1
        contexts=(controls[:,selected]-centers)/scales/np.sqrt(len(selected))
        mean=response.mean(axis=0);median=np.median(response,axis=0)
        dual=np.linalg.solve(contexts@contexts.T+np.eye(len(pairs)),response-mean)
        errors={}
        def close(name,actual,expected):
            actual=np.asarray(actual);expected=np.asarray(expected)
            assert actual.shape==expected.shape and np.isfinite(actual).all(),name
            errors[name]=float(np.max(np.abs(actual-expected),initial=0))
            assert np.allclose(actual,expected,rtol=1e-9,atol=1e-9),(fold['id'],name,errors[name])
        for key,value in [('contextCenters',centers),('contextScales',scales),('contexts',contexts),
                          ('meanResponse',mean),('medianResponse',median),('dualCoefficients',dual)]:close(key,model[key],value)
        q_counts=counts[sample_rows[fold['queryControlAccession']]].astype(float)
        assert query['libraryCounts']==int(q_counts.sum())
        q=np.log1p(q_counts/q_counts.sum()*1e6);close('control',query['control'],q)
        context=(q[selected]-centers)/scales/np.sqrt(len(selected))
        expected={'noChange':np.zeros(len(genes)),'meanResponse':mean,'medianResponse':median,'contextRidge':mean+(context@contexts.T)@dual}
        estimates={e['baseline']:e for e in query['estimates']};assert set(estimates)==set(methods)
        for name in methods:
            e=estimates[name];predicted=np.maximum(0,q+expected[name])
            close(name+'-unclippedResponse',e['unclippedResponse'],expected[name]);close(name+'-treated',e['predictedTreated'],predicted)
            close(name+'-response',e['predictedResponse'],predicted-q);close(name+'-impliedCPMSum',e['impliedCPMSum'],np.expm1(predicted).sum())
        # Only now select this fold's held-out treated outcome for scoring.
        target_counts=counts[sample_rows[fold['scoringTargetAccession']]].astype(float)
        target=np.log1p(target_counts/target_counts.sum()*1e6);truth=target-q;scores=[]
        for family,indices in [('all-source-genes',np.arange(len(genes))),('training-selected-context-genes',selected)]:
            for name in methods:
                e=estimates[name];predicted=np.asarray(e['predictedTreated'])[indices];effect=np.asarray(e['predictedResponse'])[indices]
                err=predicted-target[indices];delta=effect-truth[indices]
                pearson=None
                if np.ptp(effect)>0 and np.ptp(truth[indices])>0:
                    pearson=float(np.corrcoef(effect,truth[indices])[0,1]);assert np.isfinite(pearson)
                scores.append(dict(family=family,method=name,features=len(indices),treatedRMSE=float(np.sqrt(np.mean(err**2))),
                    treatedMAE=float(np.mean(np.abs(err))),responseRMSE=float(np.sqrt(np.mean(delta**2))),
                    responseMAE=float(np.mean(np.abs(delta))),responsePearson=pearson,impliedCPMSum=e['impliedCPMSum']))
        results.append({**base,'status':'completed','scores':scores,'numericalMaximumAbsoluteErrors':errors,'trainingDonors':model['trainingDonors']})
        print(fold['id'],'passed',flush=True)
    summaries=[]
    for cohort in sorted({f['cohort'] for f in frozen['folds']}):
        rows=[r for r in results if r['cohort']==cohort];complete=all(r['status']=='completed' for r in rows)
        for family in ('all-source-genes','training-selected-context-genes'):
            for name in methods:
                values=[s for r in rows for s in r['scores'] if s['method']==name and s['family']==family]
                averages={}
                for metric in ('treatedRMSE','treatedMAE','responseRMSE','responseMAE','responsePearson','impliedCPMSum'):
                    metric_values=[s[metric] for s in values]
                    averages[metric]=float(np.mean(metric_values)) if complete and all(v is not None for v in metric_values) else None
                summaries.append(dict(cohort=cohort,population=rows[0]['population'],treatment=rows[0]['treatment'],family=family,method=name,
                    status='complete' if complete else 'incomplete',expectedFolds=len(rows),completedFolds=len(values),equalFoldMeans=averages))
    assert all(sha(root/path)==digest for path,digest in identities.items()),'Frozen outputs changed during scoring'
    result=dict(schemaVersion=1,status='passed' if all(r['status']=='completed' for r in results) else 'incomplete',
        folds=results,summaries=summaries,foldCount=len(results),completedFolds=sum(r['status']=='completed' for r in results),
        nativeRawIngestionLinked=True,outputFreezeSHA256=sha(freeze_path),nativeOutputFreezeSHA256=sha(root/'prediction-native-output-freeze.json'),countsSHA256=sha(root/'inference-inputs/counts.npz'),
        scorerSHA256=sha(Path(__file__)),scope='All frozen known-perturbation/new-donor folds; independent numerical reconstruction and empirical prediction errors. No unseen-perturbation identity, uncertainty, causal or general biological qualification.')
    write(root/'prediction-scores.json',result)
    print(json.dumps({'status':result['status'],'folds':len(results),'summaries':len(summaries)}))

if __name__=='__main__':main()
