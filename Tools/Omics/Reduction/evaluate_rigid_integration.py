#!/usr/bin/env python3
"""Evaluate a label-free within-batch isometry projection of frozen corrections.

This is a declared development experiment, not a qualified native method.
All candidate seeds and failures remain retained under the original margins.
"""
import argparse,hashlib,json,time
from importlib.metadata import version
from pathlib import Path
import numpy as np
from scipy.linalg import orthogonal_procrustes
from check_full_integration import matrix,metrics as donor_metrics
from evaluate_ding_integration import metrics as ding_metrics

def save(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def project(x,target,batches):
    """Minimize sum ||X_b Q_b + t_b - Y_b||^2 with Q_b.T Q_b = I.

    No labels, conditions, programs, acceptance gates or inferred types enter.
    Orthogonal transformations may include reflections; no scale is fitted.
    """
    assert x.shape==target.shape and x.ndim==2 and len(batches)==len(x)
    assert np.isfinite(x).all() and np.isfinite(target).all()
    result=np.empty_like(x);fits=[];d=x.shape[1]
    for batch in np.unique(batches):
        rows=np.flatnonzero(batches==batch);assert len(rows)>d
        mx=x[rows].mean(axis=0);my=target[rows].mean(axis=0)
        a=x[rows]-mx;b=target[rows]-my;cross=a.T@b
        u,s,vt=np.linalg.svd(cross,full_matrices=False);q=u@vt
        reference,scale=orthogonal_procrustes(a,b)
        np.testing.assert_allclose(q,reference,atol=1e-10,rtol=1e-10)
        np.testing.assert_allclose(s.sum(),scale,atol=1e-8,rtol=1e-12)
        residual=float(np.linalg.norm(q.T@q-np.eye(d),ord=2));assert residual<1e-10
        t=my-mx@q;result[rows]=x[rows]@q+t
        objective=float(np.sum((result[rows]-target[rows])**2))
        optimum=float(np.sum(a*a)+np.sum(b*b)-2*s.sum())
        assert abs(objective-optimum)<1e-8*(1+np.sum(a*a)+np.sum(b*b))
        translation_objective=float(np.sum((a-b)**2));assert objective<=translation_objective+1e-7
        # Orthogonality certifies every within-level pair distance, independent of
        # cell labels. Centered Gram reconstruction is checked without an n^2 array.
        np.testing.assert_allclose((result[rows]-my).T@(result[rows]-my),q.T@(a.T@a)@q,rtol=1e-10,atol=1e-8)
        fits.append(dict(batch=str(batch),cells=len(rows),sourceMean=mx.tolist(),targetMean=my.tolist(),orthogonalMatrix=q.tolist(),translation=t.tolist(),singularValues=s.tolist(),relativeSmallestSingularValue=float(s[-1]/s[0]),orthogonalityOperatorNorm=residual,objective=objective,independentOptimalObjective=optimum,translationOnlyObjective=translation_objective,maximumWithinBatchSquaredDistanceRelativeErrorBound=residual,unconstrainedDisplacementRMS=float(np.sqrt(np.mean(np.sum((target[rows]-x[rows])**2,axis=1)))),projectedDisplacementRMS=float(np.sqrt(np.mean(np.sum((result[rows]-x[rows])**2,axis=1))))))
    return result,fits

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',required=True,choices=['kang','hagai','ding'])
    for name in ['native','inputs','protocol','out']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    protocol=json.loads(a.protocol.read_text());assert protocol['status']=='declared-before-rigid-results'
    cohort=a.cohort
    if cohort=='ding':
        inp=np.load(a.inputs/'inputs.npz',allow_pickle=False);x=inp['scores'];batch=inp['methods'];baseline=json.loads((a.inputs/'baseline.json').read_text());input_path=a.inputs/'inputs.npz'
    else:
        inp=np.load(a.inputs/'evaluation-inputs.npz',allow_pickle=False);x=inp['scores'];batch=inp['donors'];baseline=json.loads((a.inputs/'checks.json').read_text())['baseline'];input_path=a.inputs/'evaluation-inputs.npz'
    np.testing.assert_array_equal(x,matrix(a.native/'pca/scores.bin',*x.shape))
    metadata=json.loads((a.native/'pca/metadata.json').read_text());samples={s['id']:s for s in metadata['samples']};field='batchID' if cohort=='ding' else 'donorID'
    assert [samples[c['sampleID']][field] for c in metadata['cells']]==batch.tolist()
    result=dict(status='in-progress',cohort=cohort,protocolSHA256=sha(a.protocol),inputSHA256=sha(input_path),metadataSHA256=sha(a.native/'pca/metadata.json'),cells=len(x),components=x.shape[1],versions={name:version(name) for name in ['numpy','scipy','scikit-learn']},baseline=baseline,runs=[])
    save(a.out/'checks.json',result);margin=protocol['gateMargins']
    for mode in protocol['modes']:
        for seed in protocol['seeds']:
            folder=(mode+'-'+str(seed)) if cohort=='ding' or mode=='fixed' else 'seed-'+str(seed)
            root=a.native/folder;label=mode+'-'+str(seed)
            target=matrix(root/'scores.bin',*x.shape);assert json.loads((root/'metadata.json').read_text())==metadata
            y,fits=project(x,target,batch);save(a.out/(label+'-transforms.json'),fits)
            np.savez_compressed(a.out/(label+'-scores.npz'),scores=y)
            print(cohort+' '+label+' projected',flush=True);start=time.perf_counter()
            if cohort=='ding':
                m=ding_metrics(y,inp,30,4,a.out,label)
                g=dict(mixingImproved=m['sameMethodExcess']<baseline['sameMethodExcess'],cellTypesPreserved=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'],everyTypeRecallPreserved=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items()),programsPreserved=all(m['programs'][p]['spearman']>=v['spearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()),withinStratumProgramsPreserved=all(m['programs'][p]['withinStratumSpearman']>=v['withinStratumSpearman']-margin['maximumProgramSpearmanLoss'] for p,v in baseline['programs'].items()))
            else:
                m=donor_metrics(y,inp['donors'],inp['conditions'],inp['cellTypes'],inp['program'],30,a.out/(label+'-neighbors.npz'))
                g=dict(mixingImproved=m['sameDonorExcess']<baseline['sameDonorExcess'],conditionPreserved=m['conditionBalancedAccuracy']>=baseline['conditionBalancedAccuracy']-margin['maximumConditionBalancedAccuracyLoss'],programPreserved=m['programSpearman']>=baseline['programSpearman']-margin['maximumProgramSpearmanLoss'],withinStratumProgramPreserved=m['withinStratumProgramSpearman']>=baseline['withinStratumProgramSpearman']-margin['maximumProgramSpearmanLoss'],completeClassifierStrata=not m['missingClassifierStrata'])
                if baseline['cellTypeBalancedAccuracy'] is not None:
                    g['cellTypesPreserved']=m['cellTypeBalancedAccuracy']>=baseline['cellTypeBalancedAccuracy']-margin['maximumCellTypeBalancedAccuracyLoss'];g['everyTypeRecallPreserved']=all(m['cellTypeRecall'][t]>=v-margin['maximumPerTypeRecallLoss'] for t,v in baseline['cellTypeRecall'].items())
            result['runs'].append(dict(mode=mode,seed=seed,sourceScoresSHA256=sha(root/'scores.bin'),projectedScoresSHA256=sha(a.out/(label+'-scores.npz')),transformSHA256=sha(a.out/(label+'-transforms.json')),metrics=m,gates=g,evaluationSeconds=time.perf_counter()-start));save(a.out/'checks.json',result)
            print(cohort+' '+label+' measured '+json.dumps(g),flush=True)
    result['status']='all-declared-runs-measured';result['allMeasuredGatesPassed']=all(all(r['gates'].values()) for r in result['runs']);result['qualification']=protocol['qualification'];save(a.out/'checks.json',result)
if __name__=='__main__':main()
