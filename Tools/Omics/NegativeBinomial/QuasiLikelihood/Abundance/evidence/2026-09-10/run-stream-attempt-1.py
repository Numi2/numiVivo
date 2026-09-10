#!/usr/bin/env python3
"""Qualify native abundance and its complete global QL consumer on prior families."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read(path, expected=None):
    raw = path.read_bytes()
    if expected is not None:
        assert sha(raw) == expected, path
    return json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)


def row_relative(a, b):
    return np.max(np.abs(np.asarray(a)-np.asarray(b))/np.maximum(1,np.abs(b)),axis=1)


def glm_scores(y, mu, x, phi):
    denominator = 1+phi[:,None]*mu
    return np.max(np.abs(((y-mu)/denominator)@x)/np.sqrt((mu/denominator)@(x*x)),axis=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ql-root', type=Path, required=True)
    parser.add_argument('--prior-global-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--binary', required=True)
    parser.add_argument('--jobs', type=int, default=1)
    parser.add_argument('--case')
    args = parser.parse_args()
    assert 1 <= args.jobs <= 3
    old_runs = read(args.prior_global_root/'native/complete.json')
    assert len(old_runs) == 58 and all(r['status']=='passed-native-stage' for r in old_runs)
    binary_sha = subprocess.check_output(['ssh','macmini','shasum','-a','256',args.binary]).decode().split()[0]
    protocol = dict(baseRevision='d45f93bb2b84e51208a4dcc9dc0c8cd7f78fe6db',
                    protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md').read_bytes()),
                    binary=args.binary,binarySHA256=binary_sha,host='macmini',
                    priorCompleteSHA256=sha((args.prior_global_root/'native/complete.json').read_bytes()),
                    family=[dict(case=r['case'],method=r['method']) for r in old_runs])
    args.out.mkdir(parents=True,exist_ok=True)
    if (args.out/'protocol.json').exists():
        assert read(args.out/'protocol.json') == protocol
    else:
        (args.out/'protocol.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n')

    def execute(old_run):
        case, method = old_run['case'], old_run['method']
        dest = args.out/case/method
        dest.mkdir(parents=True,exist_ok=True)
        data = read(args.ql_root/case/'input.json.gz',old_run['inputSHA256'])
        ref = read(args.prior_global_root/'reference'/case/'reference.json.gz',old_run['referenceSHA256'])['results'][method]
        payload = json.dumps(dict(counts=data['counts'],design=data['design'],offsets=data['offsets'],
                                  contrast=data['contrast'],trendDispersions=ref['trendDispersions']),
                             separators=(',',':'),allow_nan=False).encode()
        if (dest/'run.json').exists():
            receipt = read(dest/'run.json')
            assert receipt['binarySHA256']==binary_sha and receipt['requestSHA256']==sha(payload)
            assert receipt['priorNativeSHA256']==old_run['outputSHA256']
            assert receipt['outputSHA256']==sha((dest/'native.json.gz').read_bytes())
            if 'comparisonSHA256' in receipt:
                assert receipt['comparisonSHA256']==sha((dest/'comparison.json.gz').read_bytes())
            return receipt
        print('start',case,method,flush=True)
        start = time.monotonic()
        with (dest/'native.log').open('wb') as log:
            process = subprocess.run(['ssh','macmini','/usr/bin/time','-l',args.binary],
                                     input=payload,stdout=subprocess.PIPE,stderr=log)
        raw = process.stdout
        packed = gzip.compress(raw,mtime=0)
        (dest/'native.json.gz').write_bytes(packed)
        receipt = dict(case=case,method=method,exitCode=process.returncode,
                       secondsIncludingTransport=time.monotonic()-start,binarySHA256=binary_sha,
                       inputSHA256=old_run['inputSHA256'],referenceSHA256=old_run['referenceSHA256'],
                       priorNativeSHA256=old_run['outputSHA256'],requestSHA256=sha(payload),
                       outputSHA256=sha(packed),logicalSHA256=sha(raw))
        if process.returncode != 0:
            receipt['status'] = 'failed-process'
        else:
            output = json.loads(raw)
            if 'error' in output:
                receipt.update(status='failed-input',error=output['error'])
            elif not output['fit']['completed']:
                receipt.update(status='incomplete-native-stage',failures=output['fit']['failures'])
            else:
                native = output['fit']
                fit = native['globalFit']
                previous = read(args.prior_global_root/'native'/case/method/'native.json.gz',old_run['outputSHA256'])['fit']
                assert not native['failures'] and not fit['failures']
                assert len(native['abundanceFits']) == len(data['counts'])
                assert len(fit['updates']) == 2
                assert fit['updates'][0]['inputScale'] == 1
                assert fit['updates'][1]['inputScale'] == fit['updates'][0]['outputScale']
                assert fit['averageQuasiDispersion'] == fit['updates'][1]['outputScale']
                assert fit['initialFits'] == previous['initialFits'], 'Abundance priors must not alter count GLMs'
                for update, old in zip(fit['updates'],previous['updates'],strict=True):
                    assert update['eligibleIndices'] == old['eligibleIndices']
                    assert update['omittedIndices'] == old['omittedIndices']
                y, x = np.asarray(data['counts'],dtype=float), np.asarray(data['design'])
                phi = np.asarray(ref['trendDispersions'])
                abundance = native['abundanceFits']
                covariates = np.asarray([a['log2CountsPerMillion'] for a in abundance])
                abundance_errors = np.abs(covariates-np.asarray(ref['abundanceCovariates']))
                # Reconstruct augmented data directly, independently of native log scaling.
                libraries = np.exp(data['offsets'])
                prior = 2*libraries/np.mean(libraries)
                z, adjusted = y+prior[None,:], libraries+2*prior
                beta = np.asarray([a['logProportion'] for a in abundance])
                def abundance_score(coefficient):
                    means = np.exp(coefficient[:,None])*adjusted[None,:]
                    denominator = 1+phi[:,None]*means
                    return np.sum((z-means)/denominator,axis=1),np.sum(means/denominator,axis=1)
                ab_score, ab_information = abundance_score(beta)
                ab_scaled = np.abs(ab_score)/np.sqrt(ab_information)
                lower_score,_ = abundance_score(beta-2e-8*np.log(2))
                upper_score,_ = abundance_score(beta+2e-8*np.log(2))
                brackets = (lower_score >= 0)&(upper_score <= 0)
                final_means = np.asarray([f['means'] for f in fit['refittedFits']])
                final_score = glm_scores(y,final_means,x,phi/fit['averageQuasiDispersion'])
                reported_score = np.asarray([f['maximumScaledScore'] for f in fit['refittedFits']])
                mean_errors = row_relative(final_means,[f['means'] for f in previous['refittedFits']])
                residual_errors = {
                    field:np.asarray([abs(new[field]-old[field])/max(1,abs(old[field]))
                                      for new,old in zip(fit['adjustedResiduals'],previous['adjustedResiduals'],strict=True)])
                    for field in ['deviance','degreesOfFreedom','quasiDispersion']}
                scale_error = abs(fit['averageQuasiDispersion']-previous['averageQuasiDispersion'])/max(1,abs(previous['averageQuasiDispersion']))
                failed = ((abundance_errors>2e-8)|~brackets|(ab_scaled>1e-7+1e-10)|
                          (mean_errors>2e-7)|(final_score>1e-7+1e-10)|(reported_score>1e-7)|
                          (np.abs(final_score-reported_score)>1e-10))
                for errors in residual_errors.values():
                    failed |= errors>2e-7
                finite = [covariates,abundance_errors,ab_scaled,mean_errors,final_score,*residual_errors.values()]
                assert all(np.all(np.isfinite(a)) for a in finite)
                table = [dict(featureIndex=feature,nativeAbundance=float(covariates[g]),
                              referenceAbundance=ref['abundanceCovariates'][g],abundanceAbsoluteError=float(abundance_errors[g]),
                              abundanceScaledScore=float(ab_scaled[g]),scoreBracketPassed=bool(brackets[g]),
                              refitMeanRelativeError=float(mean_errors[g]),finalScaledScore=float(final_score[g]),
                              residualRelativeErrors={k:float(v[g]) for k,v in residual_errors.items()},
                              passed=not bool(failed[g])) for g,feature in enumerate(data['featureIndices'])]
                comparisons = gzip.compress(json.dumps(table,separators=(',',':'),allow_nan=False).encode(),mtime=0)
                (dest/'comparison.json.gz').write_bytes(comparisons)
                receipt.update(status='passed' if not np.any(failed) and scale_error<=2e-7 else 'failed-comparison',
                               genes=len(table),failedGenes=int(np.sum(failed)),comparisonSHA256=sha(comparisons),
                               maximumAbundanceAbsoluteError=float(np.max(abundance_errors)),
                               maximumAbundanceScaledScore=float(np.max(ab_scaled)),
                               maximumFinalScaledScore=float(np.max(final_score)),
                               maximumRefitMeanRelativeError=float(np.max(mean_errors)),
                               maximumResidualRelativeError={k:float(np.max(v)) for k,v in residual_errors.items()},
                               scale=fit['averageQuasiDispersion'],priorScale=previous['averageQuasiDispersion'],scaleRelativeError=scale_error,
                               abundanceScoreEvaluations=sum(a['scoreEvaluations'] for a in abundance),
                               maximumAbundanceScoreEvaluations=max(a['scoreEvaluations'] for a in abundance),
                               maximumAbundanceBracketWidth=max(a['logProportionBracketWidth'] for a in abundance),
                               omittedIndices=[u['omittedIndices'] for u in fit['updates']])
        (dest/'run.json').write_text(json.dumps(receipt,sort_keys=True,indent=2,allow_nan=False)+'\n')
        print(json.dumps(receipt),flush=True)
        return receipt

    selected = [r for r in old_runs if args.case is None or r['case']==args.case]
    assert selected
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(execute,selected))
    (args.out/('pilot-complete.json' if args.case else 'complete.json')).write_text(json.dumps(results,sort_keys=True,indent=2)+'\n')
    return 0 if all(r['status']=='passed' for r in results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
