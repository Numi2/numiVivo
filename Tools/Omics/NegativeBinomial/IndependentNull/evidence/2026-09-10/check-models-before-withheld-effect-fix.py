#!/usr/bin/env python3
"""Independent count/design, score, likelihood, information, tail and BH checks."""
import argparse, collections, gzip, hashlib, json
from pathlib import Path
import numpy as np
from scipy.special import gammaln, erfc
from scipy.stats import chi2, f as f_distribution

p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
p.add_argument('--available',action='store_true');a=p.parse_args();root=a.root
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();checker_sha=sha(Path(__file__))
read=lambda p:json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
audit=read(root/'count-audit.json');assert sha(root/'reference-counts.npz')==audit['referenceCountsSHA256']
with np.load(root/'reference-counts.npz',allow_pickle=True) as archive:
    all_y=archive['counts'];donors=archive['donors'].tolist();feature_ids=archive['features'].tolist()
manifest=read(root/'requests/manifest.json');summaries=[]
def bh(values):
    values=np.asarray(values);n=len(values)
    if not n:return values
    order=np.argsort(values,kind='stable');adjusted=np.minimum.accumulate((values[order]*n/np.arange(1,n+1))[::-1])[::-1]
    result=np.empty(n);result[order]=np.minimum(1,adjusted);return result
for case in manifest:
    directory=root/'native'/(case['cohort']+'-'+case['method']);receipt_path=directory/'receipt.json'
    if not receipt_path.exists():
        assert a.available,(case,'not terminal');continue
    receipt=read(receipt_path);assert sha(directory/'output.json.gz')==receipt['outputSHA256']
    destination=directory/'model-check.json'
    if destination.exists():
        previous=read(destination)
        if previous['checkerSHA256']==checker_sha and previous['outputSHA256']==receipt['outputSHA256']:
            summaries.append(previous);continue
    out=read(directory/'output.json.gz');r=out.get('result');checks=[];per_gene={}
    def check(name,actual,expected=0,tolerance=0,relative=False):
        actual=np.asarray(actual,dtype=float);expected=np.asarray(expected,dtype=float)
        errors=np.abs(actual-expected)
        if relative:errors=errors/np.maximum(1,np.abs(expected))
        error=float(np.max(errors,initial=0));passed=bool(np.isfinite(error) and error<=tolerance)
        checks.append(dict(name=name,error=error,tolerance=tolerance,passed=passed))
    if r is not None:
        d=r['design'];meta=read(root/'reference-inputs'/case['cohort']/'input.json')
        assert [f['featureID'] for f in r['features']]==feature_ids
        for key in ['columnNames','contrast','sourcePseudobulkIndices','libraryCounts','referenceFeatureIndices']:assert d[key]==meta[key],key
        assert d['rows']==meta['design'] and d['excludedSmallPseudobulkIndices']==[]
        assert d['controlReplicates']==d['treatmentReplicates']==6
        assert [o['donorID'] for o in d['observations']]==meta['donorIDs']
        y=all_y[[donors.index(donor) for donor in meta['donorIDs']]].T.astype(float)
        x=np.asarray(d['rows']);contrast=np.asarray(d['contrast']);factors=np.asarray(d['sizeFactorValues'])
        check('independent-native-offsets',factors,meta['sizeFactors'],2e-13,True)
        check('feature-total-counts',[f['totalCounts'] for f in r['features']],y.sum(axis=1))
        check('expressing-observations',[f['expressingPseudobulks'] for f in r['features']],(y>0).sum(axis=1))
        check('normalized-means',[f['meanNormalizedCount'] for f in r['features']],(y/factors).mean(axis=1),2e-12,True)
        features=r['features'];nb=r['negativeBinomial'];diagnostics=nb['features']
        assert [f['featureIndex'] for f in features]==[f['featureIndex'] for f in diagnostics]==list(range(len(feature_ids)))
        eligible=(y.sum(axis=1)>=10)&((y>0).sum(axis=1)>=3)
        assert int(eligible.sum())==meta['eligibleFeatures']
        assert all((f['status']=='filteredLowExpression')==bool(not eligible[g]) for g,f in enumerate(features))
        tested=[i for i,f in enumerate(features) if f['status']=='tested']
        assert len(tested)==r['testedFeatures']
        assert all(f.get('pValue') is None and f.get('adjustedPValue') is None for f in features if f['status']!='tested')
        ql=nb.get('quasiLikelihood');ql_tests={}
        if ql and ql.get('inference'):
            inf=ql['inference'];ql_tests={g:t for g,t in zip(ql['featureIndices'],inf['tests']) if t is not None}
            assert sorted(ql_tests)==tested and inf['poissonBound']=='not-applicable-to-modern-adjusted-QL'
            check('ordinary-DF-cap',inf['ordinaryResidualDFCap'],len(ql['featureIndices'])*(x.shape[0]-x.shape[1]))
        # Check each retained converged fit, including boundary fits without tests.
        rows=[g for g,f in enumerate(diagnostics) if f.get('finalFit',{}).get('converged')]
        if rows:
            fits=[diagnostics[g]['finalFit'] for g in rows];mu=np.asarray([f['means'] for f in fits]);yy=y[rows]
            phi=np.asarray([f['dispersion'] for f in fits])[:,None];beta=np.asarray([f['coefficients'] for f in fits])
            check('linear-predictor-means',mu,np.exp(beta@x.T)*factors,2e-12,True)
            weight=mu/(1+phi*mu);scores=((yy-mu)/(1+phi*mu))@x
            scaled=np.max(np.abs(scores)/np.sqrt(weight@(x*x)),axis=1)
            check('full-scaled-score',scaled,0,1.01e-7)
            size=1/phi
            ll=(gammaln(yy+size)-gammaln(size)-gammaln(yy+1)+yy*np.log(phi*mu)-(yy+size)*np.log1p(phi*mu)).sum(axis=1)
            check('NB-log-likelihood',[f['logLikelihood'] for f in fits],ll,2e-6,True)
            information=np.einsum('gn,np,nq->gpq',weight,x,x);inverse=np.linalg.inv(information)
            se=np.sqrt(np.einsum('p,gpq,q->g',contrast,inverse,contrast));effect=beta@contrast
            check('conditional-effect',[f['effect'] for f in fits],effect,2e-12,True)
            check('conditional-information-SE',[f['standardError'] for f in fits],se,2e-10,True)
            check('reported-log2-effect',[features[g]['log2FoldChange'] for g in rows],effect/np.log(2),2e-12,True)
            per_gene.update(fittedFeatureIndices=rows,fullScaledScores=scaled.tolist(),referenceLogLikelihood=ll.tolist())
        null_rows=[g for g in tested if diagnostics[g].get('likelihoodRatioFit')]
        if null_rows:
            fits=[diagnostics[g]['finalFit'] for g in null_rows];lr=[diagnostics[g]['likelihoodRatioFit'] for g in null_rows]
            mu=np.asarray([f['means'] for f in fits]);null_mu=np.asarray([f['nullMeans'] for f in lr]);yy=y[null_rows]
            phi=np.asarray([f['dispersion'] for f in fits])[:,None];null_beta=np.asarray([f['nullCoefficients'] for f in lr])
            nx=x[:,np.arange(x.shape[1])!=1]
            check('null-constraint',null_beta@contrast,0,2e-12)
            check('null-linear-predictor',null_mu,np.exp(null_beta@x.T)*factors,2e-12,True)
            scaled=np.max(np.abs(((yy-null_mu)/(1+phi*null_mu))@nx)/np.sqrt((null_mu/(1+phi*null_mu))@(nx*nx)),axis=1)
            check('null-scaled-score',scaled,0,1.01e-7)
            delta=mu-null_mu
            reference_lr=2*(yy*np.log1p(delta/null_mu)-(yy+1/phi)*np.log1p(phi*delta/(1+phi*null_mu))).sum(axis=1)
            check('gamma-cancelled-LR',[v['rawStatistic'] for v in lr],reference_lr,2e-6,True)
            check('nonnegative-LR',[v['statistic'] for v in lr],np.maximum(0,reference_lr),2e-6,True)
            per_gene.update(nullFeatureIndices=null_rows,nullScaledScores=scaled.tolist(),referenceLR=reference_lr.tolist())
        if case['method']=='wald':
            expected=[float(erfc(abs(f['zStatistic'])/np.sqrt(2))) for f in features if f['status']=='tested']
            check('normal-tail',[features[g]['pValue'] for g in tested],expected,2e-12)
        elif case['method']=='lrt':
            expected=chi2.sf([diagnostics[g]['likelihoodRatioFit']['statistic'] for g in tested],1)
            check('chi-square-tail',[features[g]['pValue'] for g in tested],expected,2e-12)
        elif ql_tests:
            positions={g:i for i,g in enumerate(ql['featureIndices'])};mod=ql['moderation']
            post=np.asarray([mod['posteriorVariances'][positions[g]] for g in tested])
            actual=np.asarray([ql_tests[g]['fStatistic'] for g in tested]);lr=np.asarray([ql_tests[g]['likelihoodRatio']['statistic'] for g in tested])
            df=np.asarray([ql_tests[g]['denominatorDegreesOfFreedom'] for g in tested])
            assert np.all(df>0) and np.all(df<=ql['inference']['ordinaryResidualDFCap'])
            check('F-arithmetic',actual,lr/post,2e-12,True)
            expected=f_distribution.logsf(actual,1,df)
            check('conditional-log-F-tail',[ql_tests[g]['logPValue'] for g in tested],expected,2e-7)
            check('F-probability',[features[g]['pValue'] for g in tested],np.exp(expected),2e-10)
            assert all(features[g].get('standardError') is None and features[g].get('intervalLower') is None for g in tested)
            # Successful cohort output retains test DF but not individual moment
            # quadrature records. This checks tails conditional on that DF and its
            # family cap, not a new independent moment/DF requalification.
        probabilities=[features[g]['pValue'] for g in tested]
        check('BH',[features[g]['adjustedPValue'] for g in tested],bh(probabilities),2e-12)
        statuses=dict(collections.Counter(f['status'] for f in features))
    else:
        rows=[];tested=[];statuses={};checks=[]
    result=dict(cohort=case['cohort'],method=case['method'],checkerSHA256=checker_sha,outputSHA256=receipt['outputSHA256'],
        status='passed' if r is not None and all(c['passed'] for c in checks) else 'method-failed' if r is None else 'failed',
        methodError=out.get('error'),checks=checks,featureStatuses=statuses,checkedFits=len(rows),testedFeatures=len(tested),
        scope='Count/design identity and conditional numerical arithmetic; QL tails use reported denominator DF. Statistical calibration is a separate descriptive null experiment.')
    (directory/'model-gene-checks.json.gz').write_bytes(gzip.compress(json.dumps(per_gene,separators=(',',':')).encode(),mtime=0))
    destination.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');summaries.append(result)
    print(case['cohort'],case['method'],result['status'],[c for c in checks if not c['passed']],flush=True)
summary=dict(cases=len(summaries),status='passed' if len(summaries)==27 and all(s['status']=='passed' for s in summaries) else 'partial-or-failed',
    results=[{k:v for k,v in s.items() if k!='checks'} for s in summaries],checkerSHA256=checker_sha)
(root/('model-check-available.json' if a.available else 'model-check-summary.json')).write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
