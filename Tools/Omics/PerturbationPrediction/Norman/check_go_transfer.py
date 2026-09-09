#!/usr/bin/env python3
"""Independent numerical, repeatability and excluded-outcome checks for GO transfer."""
import argparse
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from go_transfer import (LAMBDAS, METHODS, fit, independent_prediction, reference,
                         solve_weights, select)
from combinations import counts, load, sha, write

p=argparse.ArgumentParser(description=__doc__)
for key in ('reference','prepared','repeat-prepared','first','repeat','out'):
    p.add_argument('--'+key,type=Path,required=True)
a=p.parse_args()
with threadpool_limits(limits=1):
    ref=reference(a.reference)
    for path in a.prepared.iterdir():
        assert path.read_bytes()==(a.repeat_prepared/path.name).read_bytes(),path.name
    for name in ('predictions.npz','folds.json','receipt.json'):
        assert (a.first/name).read_bytes()==(a.repeat/name).read_bytes(),name
    prediction=load(a.first/'predictions.npz')
    desc=load(a.prepared/'kernel.npz')
    infos=json.loads((a.first/'folds.json').read_text())
    targets=sorted(t for t in ref['conditions'] if t!='control' and '_' not in t)
    assert len(targets)==105 and prediction['targetIDs'].tolist()==targets
    np.testing.assert_array_equal(ref['feature_ids'],prediction['featureIDs'])
    assert len(infos)==len(targets)
    for i,t in enumerate(targets):
        data=select(ref,t)
        assert infos[i]['trainingTargets']==data['trainingTargets']
        assert infos[i]['excludedConditions']==data['excluded']
        assert t not in infos[i].get('descriptorTrainingTargets',[])
        assert len(data['trainingTargets'])==104
        assert bool(prediction['supported'][i])==(t in desc['targetIDs'])
        for method in METHODS:
            row=prediction[method][i]
            if method in METHODS[2:] and not prediction['supported'][i]:
                assert np.isnan(row).all()
            else:
                assert np.isfinite(row).all() and (row>=0).all()
        for values in infos[i].get('weights',{}).values():
            assert abs(sum(values)-1)<1e-10
    errors=[]
    for target in ('AHR','ZNF318'):
        data=select(ref,target);names=desc['targetIDs'].tolist()
        use=[i for i,t in enumerate(data['trainingTargets']) if t in names]
        rows=[names.index(data['trainingTargets'][i]) for i in use]
        k=desc['kernel'][np.ix_(rows,rows)]
        raw=counts(data['trainingCounts']);ctrl=counts(data['control'])
        baseline=np.log1p(ctrl/ctrl.sum()*1e6)
        y=np.log1p(raw/raw.sum(axis=1,keepdims=True)*1e6)[use]-baseline
        for lam in LAMBDAS:
            _,loo=solve_weights(k,np.zeros(len(k)),lam)
            for i in (0,len(k)//2,len(k)-1):
                keep=np.arange(len(k))!=i
                expected=independent_prediction(k[np.ix_(keep,keep)],k[i,keep],y[keep],lam)
                actual=loo[i]@y
                np.testing.assert_allclose(actual,expected,rtol=1e-9,atol=1e-10)
                errors.append(float(np.max(abs(actual-expected))))
    mutation_targets=('AHR','IER5L','KIAA1804','ZNF318')
    for target in mutation_targets:
        i=targets.index(target);data=select(ref,target)
        changed=dict(ref,counts=ref['counts'].copy())
        names=ref['conditions'].tolist()
        changed['counts'][[names.index(t) for t in data['excluded']]]=1
        mutated=select(changed,target)
        values,baseline,top,info=fit(mutated,desc)
        for method,value in values.items():
            np.testing.assert_array_equal(prediction[method][i],value)
        np.testing.assert_array_equal(prediction['baseline'][i],baseline)
        np.testing.assert_array_equal(prediction['trainingTop1000'][i],top)
        recorded=dict(infos[i]);recorded.pop('excludedOutcomeMutationInputsExact')
        assert info==recorded
    out=dict(status='passed',allPreparationsAndPredictionsByteExact=True,
             folds=105,supportedTargets=int(prediction['supported'].sum()),
             allFoldMembershipsAndUnsupportedOutputsChecked=True,
             realExplicitLOORefits=len(errors),maximumLOOError=max(errors),
             excludedOutcomeMutationRefits=list(mutation_targets),
             predictionsSHA256=sha(a.first/'predictions.npz'),
             descriptorSHA256=sha(a.prepared/'kernel.npz'),
             implementationSHA256=sha(Path(__file__)),
             qualification='Numerical, source identity and outcome isolation evidence; independent biological validation remains separate.')
    if a.out.exists():raise ValueError('output already exists')
    write(a.out,out);print(json.dumps(out))
