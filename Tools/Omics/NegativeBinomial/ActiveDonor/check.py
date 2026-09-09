#!/usr/bin/env python3
"""Independent real-data checks of explicit gene-specific active-donor NB fits.

Counts must be genes x pseudobulks independently aggregated from the archived
H5AD (Bioconductor/prepare.py). No inferred source counts or synthetic evidence.
"""
import argparse, hashlib, json, warnings
from collections import Counter, defaultdict
from importlib.metadata import version
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.optimize import minimize_scalar
from scipy.special import erfc
from statsmodels.stats.multitest import multipletests

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--baseline',type=Path,required=True)
p.add_argument('--active',type=Path,required=True)
p.add_argument('--counts',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args()
assert not a.out.exists(),a.out
base=json.loads(a.baseline.read_text()); report=json.loads(a.active.read_text())
for key in base.keys()-{'contrasts'}: assert base[key]==report[key],key
b=base['contrasts'][0];c=report['contrasts'][0];d=c['design'];nb=c['negativeBinomial']
assert len(base['contrasts'])==len(report['contrasts'])==1
assert b['design']==d and b['negativeBinomial']['trend']==nb['trend']
request=dict(c['request']);opts=dict(request['negativeBinomialOptions'])
assert opts.pop('zeroTotalDonorPolicy')=='activeDonorProfile'
request['negativeBinomialOptions']=opts
assert request==b['request'] and request['design']=='pairedDonors'
assert not request['adjustForBatch'],'This reference qualifies the unbatched paired benchmarks only'
counts=pd.read_csv(a.counts,sep='\t',index_col=0,float_precision='round_trip')
obs=d['observations'];samples=[o['sampleIDs'][0] for o in obs]
assert all(len(o['sampleIDs'])==1 for o in obs)
assert list(counts.columns)==samples
assert list(counts.index)==report['pseudobulk']['featureIDs']
assert list(counts.index)==[f['featureID'] for f in c['features']]
yall=counts.to_numpy().T
# The independent Hagai TSV stores exact integers with a '.0' suffix.
# Require exact integrality and safe range before changing representation.
assert np.all(np.isfinite(yall)) and np.all(yall>=0) and np.all(yall<=2**53-1)
assert np.all(yall==np.floor(yall)), 'Fractional counts are not accepted'
yall=yall.astype(np.int64)
m=report['pseudobulk']['matrix']
# Full original count identity, including filtered genes and all observations.
for i,source in enumerate(d['sourcePseudobulkIndices']):
    row=np.zeros(yall.shape[1],dtype=np.int64)
    for k in range(m['rowOffsets'][source],m['rowOffsets'][source+1]):
        row[m['featureIndices'][k]]=m['counts'][k]
    assert np.array_equal(row,yall[i])
donors=np.array([o['donorID'] for o in obs]);conditions=np.array([o['condition'] for o in obs])
all_donors=sorted(set(donors));sf=np.array(d['sizeFactorValues']);offset=np.log(sf)
for donor in all_donors:
    mask=donors==donor
    assert mask.sum()==2 and set(conditions[mask])=={request['controlCondition'],request['treatmentCondition']}
summary=dict(status='in-progress',originalStatusCounts=dict(Counter(f['status'] for f in b['features'])),
    activeStatusCounts=dict(Counter(f['status'] for f in c['features'])),exactSourceCounts=True,
    exactGlobalDesign=True,exactFullSupportPrior=True,fullSupportFitsUnchanged=0,independentFitsChecked=0,
    supportOutcomes=Counter(),activeDonorCounts=Counter(),maxErrors=defaultdict(float),
    profileChecks=[],referenceWarnings=[],genes=[])
profile_groups=Counter();tested=[]
for i,(f,diag) in enumerate(zip(c['features'],nb['features'])):
    y=yall[:,i];old=b['features'][i];old_diag=b['negativeBinomial']['features'][i]
    eligible=int(y.sum())>=request['minimumFeatureCounts'] and np.count_nonzero(y)>=request['minimumExpressingPseudobulks']
    zero=[donor for donor in all_donors if np.all(y[donors==donor]==0)]
    resolution=diag.get('supportResolution')
    if not eligible or not zero:
        assert resolution is None
        assert diag==old_diag
        assert {k:v for k,v in f.items() if k!='adjustedPValue'}=={k:v for k,v in old.items() if k!='adjustedPValue'}
        if f['status']=='tested': summary['fullSupportFitsUnchanged']+=1
    else:
        assert resolution is not None and old['status']=='rankDeficientSupport'
        retained=np.flatnonzero(~np.isin(donors,zero));active=sorted(set(donors[retained]));n=len(retained)
        x=np.column_stack([np.ones(n),conditions[retained]==request['treatmentCondition']]+[donors[retained]==v for v in active[1:]])
        contrast=np.zeros(x.shape[1]);contrast[1]=1
        assert resolution['excludedDonorIDs']==zero and resolution['activeDonorIDs']==active
        assert resolution['retainedObservationIndices']==retained.tolist()
        assert resolution['excludedObservationIndices']==np.flatnonzero(np.isin(donors,zero)).tolist()
        mean=float(np.mean(y[retained]/sf[retained])) if n else 0.
        assert np.isclose(resolution['meanNormalizedCount'],mean,rtol=2e-15,atol=1e-14)
        if len(active)<request['minimumReplicatesPerCondition']:
            outcome='insufficientReplication';assert f['status']=='insufficientActiveDonors'
            assert 'finalFit' not in diag
        else:
            assert resolution['columnNames']==['intercept','treatment-minus-control']+['donor:'+v for v in active[1:]]
            assert np.array_equal(resolution['rows'],x) and np.array_equal(resolution['contrast'],contrast)
            assert resolution['residualDegreesOfFreedom']==n-x.shape[1]
            outcome='rankDeficientPositiveSupport' if np.linalg.matrix_rank(x[y[retained]>0])<x.shape[1] else 'ready'
            if outcome!='ready':assert f['status']=='rankDeficientSupport' and 'finalFit' not in diag
        assert resolution['outcome']==outcome
        summary['supportOutcomes'][outcome]+=1;summary['activeDonorCounts'][len(active)]+=1
        entry=dict(featureID=f['featureID'],outcome=outcome,activeDonors=len(active),status=f['status'])
        summary['genes'].append(entry)
        if 'finalFit' in diag:
            assert outcome=='ready'
            expected_target=nb['trend']['intercept']+nb['trend']['inverseMeanCoefficient']/mean
            assert np.isclose(diag['trendDispersion'],expected_target,rtol=2e-15,atol=1e-14)
            yy=y[retained];oo=offset[retained];alpha=diag['finalDispersion'];native=diag['finalFit']
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                fit=sm.GLM(yy,x,offset=oo,family=sm.families.NegativeBinomial(alpha=alpha)).fit(maxiter=300,tol=1e-11)
            summary['referenceWarnings'] += [dict(gene=f['featureID'],message=str(w.message)) for w in caught]
            assert fit.converged
            summary['independentFitsChecked']+=1
            effect=float(contrast@fit.params);se=float(np.sqrt(contrast@fit.cov_params()@contrast));z=effect/se
            w=fit.mu/(1+alpha*fit.mu);info=x.T@(w[:,None]*x)
            h=w*np.einsum('ij,jk,ik->i',x,np.linalg.inv(info),x)
            pearson=(yy-fit.mu)/np.sqrt(fit.mu+alpha*fit.mu**2)
            cooks=pearson**2*h/(x.shape[1]*(1-h)**2)
            errors=dict(effect=abs(effect-native['effect']),standardError=abs(se-native['standardError']),
                relativeMeans=float(np.max(abs(fit.mu-np.array(native['means']))/np.maximum(1,fit.mu))),
                leverage=float(np.max(abs(h-native['leverage']))),
                relativeCooks=float(np.max(abs(cooks-native['cooksDistances'])/np.maximum(1,cooks))))
            if f['status']=='tested':
                errors['z']=abs(z-f['zStatistic']);errors['p']=abs(float(erfc(abs(z)/np.sqrt(2)))-f['pValue'])
            for key,value in errors.items():summary['maxErrors'][key]=max(summary['maxErrors'][key],value)
            assert all(v<2e-5 for v in errors.values()),(f['featureID'],errors)
            entry['referenceErrors']=errors
            # Six source-order examples per active-donor stratum. Do not exclude
            # difficult optimizer outcomes silently; a failed assertion retains the log.
            if not diag['dispersionOutlier'] and profile_groups[len(active)]<6:
                target=diag['trendDispersion'];prior=nb['trend']['priorLogVariance']
                def objective(t):
                    alpha=float(np.exp(t))
                    ref=sm.GLM(yy,x,offset=oo,family=sm.families.NegativeBinomial(alpha=alpha)).fit(maxiter=300,tol=1e-11)
                    assert ref.converged
                    weights=ref.mu/(1+alpha*ref.mu)
                    sign,ld=np.linalg.slogdet(x.T@(weights[:,None]*x));assert sign==1
                    return -(ref.llf-.5*ld-(t-np.log(target))**2/(2*prior))
                opt=minimize_scalar(objective,bounds=(np.log(1e-8),np.log(100)),method='bounded',options={'xatol':1e-7})
                assert opt.success
                gap=abs(objective(np.log(alpha))-opt.fun)
                assert gap<1e-4,(f['featureID'],gap)
                summary['profileChecks'].append(dict(featureID=f['featureID'],activeDonors=len(active),objectiveGap=gap))
                profile_groups[len(active)]+=1
    if f['status']=='tested':tested.append(f)
    else:assert 'pValue' not in f and 'adjustedPValue' not in f
assert len(tested)==c['testedFeatures']
q=multipletests([f['pValue'] for f in tested],method='fdr_bh')[1]
assert np.max(abs(q-[f['adjustedPValue'] for f in tested]))<1e-12
summary.update(status='passed-active-donor-numerical-checks',testedGenes=len(tested),exactBH=True,
    sourceHashes={key:hashlib.sha256(path.read_bytes()).hexdigest() for key,path in [('baseline',a.baseline),('active',a.active),('counts',a.counts)]},
    versions={key:version(key) for key in ['numpy','pandas','scipy','statsmodels']},
    boundary='Conditional numerical and bookkeeping checks on these studies; no selection-adjusted FDR, biological ground truth or generalization qualification')
a.out.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='genes'},indent=2))
