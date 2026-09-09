#!/usr/bin/env python3
"""Pinned GO retrieval, annotation-kernel transfer, nested selection and sealed scoring."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
from sklearn.kernel_ridge import KernelRidge
from threadpoolctl import threadpool_limits
from combinations import REFERENCE_SHA256, counts, load, metrics, sha, write
from unseen_targets import select

METHODS = ('noChange', 'meanSingleResponse', 'meanSupportedResponse',
           'goRidgeFixed', 'goRidgeNested', 'shuffledGoRidgeNested')
LAMBDAS = (0.01, 0.1, 1.0, 10.0, 100.0)
ROOT_TERMS = {'GO:0008150', 'GO:0003674', 'GO:0005575'}


def identity():
    root = Path(__file__).parent
    return {n: sha(root / n) for n in ('go_transfer.py', 'GO_TRANSFER_PROTOCOL.md',
                                      'combinations.py', 'unseen_targets.py')}


def reference(path):
    if sha(path) != REFERENCE_SHA256:
        raise ValueError('reference fingerprint mismatch')
    return load(path)


def get(url):
    with urllib.request.urlopen(url, timeout=45) as response:
        raw = response.read(4_000_001)
    if len(raw) > 4_000_000:
        raise ValueError('annotation response limit')
    json.loads(raw)
    return raw


def fetch(ref_path, out):
    ref = reference(ref_path)
    out.mkdir(parents=True, exist_ok=False)
    rawdir = out / 'responses'
    rawdir.mkdir()
    before = get('https://mygene.info/v3/metadata')
    (out / 'metadata-before.json').write_bytes(before)
    targets = sorted(t for t in ref['conditions'] if t != 'control' and '_' not in t)
    symbols = ref['gene_symbols'].tolist()
    def one(target):
        record = dict(target=str(target), retrievedUTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
        if symbols.count(target) != 1:
            return dict(record, status='unresolved-source-symbol')
        gene = str(ref['feature_ids'][symbols.index(target)])
        url = 'https://mygene.info/v3/gene/' + urllib.parse.quote(gene) + '?fields=symbol,taxid,ensembl.gene,entrezgene,go'
        record.update(sourceFeatureID=gene, requestURL=url)
        try:
            raw = get(url)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return dict(record, status='annotation-not-found')
            raise
        path = rawdir / (gene + '.json')
        path.write_bytes(raw)
        return dict(record, status='retrieved', response=path.relative_to(out).as_posix(),
                    responseSHA256=sha(path), bytes=len(raw))
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(one, targets))
    write(out / 'requests.json', records)
    after = get('https://mygene.info/v3/metadata')
    (out / 'metadata-after.json').write_bytes(after)
    a, b = json.loads(before), json.loads(after)
    if a['build_version'] != b['build_version']:
        raise ValueError('annotation service build changed during panel retrieval')
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, implementation=identity(),
          requestsSHA256=sha(out / 'requests.json'), metadataBeforeSHA256=sha(out / 'metadata-before.json'),
          metadataAfterSHA256=sha(out / 'metadata-after.json'), buildVersion=a['build_version'],
          targets=len(records), retrieved=sum(r['status'] == 'retrieved' for r in records)))


def prepare(ref_path, source, out):
    ref = reference(ref_path)
    receipt = json.loads((source / 'receipt.json').read_text())
    if sha(source / 'requests.json') != receipt['requestsSHA256']:
        raise ValueError('annotation request fingerprint mismatch')
    for name, key in [('metadata-before.json','metadataBeforeSHA256'), ('metadata-after.json','metadataAfterSHA256')]:
        if sha(source / name) != receipt[key]:
            raise ValueError('annotation metadata fingerprint mismatch')
    requests = json.loads((source / 'requests.json').read_text())
    targets = sorted(t for t in ref['conditions'] if t != 'control' and '_' not in t)
    if [r['target'] for r in requests] != targets:
        raise ValueError('annotation target coverage differs')
    symbols = ref['gene_symbols'].tolist()
    coverage, annotations = [], {}
    for request in requests:
        target = request['target']
        item = dict(target=target, status=request['status'])
        if request['status'] == 'retrieved':
            path = source / request['response']
            if sha(path) != request['responseSHA256']:
                raise ValueError('annotation response fingerprint mismatch')
            gene = str(ref['feature_ids'][symbols.index(target)])
            value = json.loads(path.read_text())
            ensembles = value.get('ensembl', [])
            if isinstance(ensembles, dict):
                ensembles = [ensembles]
            ids = {x['gene'] for x in ensembles if isinstance(x, dict) and isinstance(x.get('gene'), str)}
            if value.get('taxid') != 9606 or gene not in ids or gene != request['sourceFeatureID']:
                item['status'] = 'identity-mismatch'
            else:
                terms, records = set(), []
                for namespace in ('BP', 'MF', 'CC'):
                    entries = value.get('go', {}).get(namespace, [])
                    if isinstance(entries, dict):
                        entries = [entries]
                    for entry in entries:
                        qualifier = entry.get('qualifier', '')
                        qualifiers = qualifier if isinstance(qualifier, list) else qualifier.split('|')
                        term = entry.get('id', '')
                        if 'NOT' in qualifiers or term in ROOT_TERMS:
                            continue
                        if len(term) != 10 or not term.startswith('GO:') or not term[3:].isdigit():
                            raise ValueError('malformed GO identifier')
                        terms.add(term)
                        records.append(dict(entry, namespace=namespace))
                item.update(sourceFeatureID=gene, currentSymbol=value.get('symbol'),
                            entrezID=value.get('entrezgene'), terms=len(terms),
                            status='supported' if terms else 'no-usable-GO-annotations')
                if terms:
                    annotations[target] = dict(terms=sorted(terms), annotations=records)
        coverage.append(item)
    supported = sorted(annotations)
    sets = [set(annotations[t]['terms']) for t in supported]
    kernel = np.array([[len(a & b) / len(a | b) for b in sets] for a in sets])
    vocabulary = sorted(set.union(*sets))
    binary = np.array([[term in s for term in vocabulary] for s in sets], dtype=np.int64)
    intersection = binary @ binary.T
    sizes = binary.sum(axis=1)
    independent = intersection / (sizes[:, None] + sizes[None, :] - intersection)
    np.testing.assert_array_equal(kernel, independent)
    if np.linalg.eigvalsh(kernel).min() < -1e-10:
        raise ValueError('Jaccard kernel is not positive semidefinite')
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out / 'kernel.npz', targetIDs=np.array(supported), kernel=kernel)
    write(out / 'annotations.json', annotations)
    write(out / 'coverage.json', coverage)
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, implementation=identity(),
          sourceReceiptSHA256=sha(source / 'receipt.json'), buildVersion=receipt['buildVersion'],
          kernelSHA256=sha(out / 'kernel.npz'), annotationsSHA256=sha(out / 'annotations.json'),
          coverageSHA256=sha(out / 'coverage.json'), supportedTargets=len(supported), terms=len(vocabulary),
          independentKernelExact=True, minimumKernelEigenvalue=float(np.linalg.eigvalsh(kernel).min())))


def solve_weights(kernel, query, lam):
    n = len(kernel)
    inverse = np.linalg.solve(kernel + np.eye(n) * lam, np.eye(n))
    u = inverse.sum(axis=1)
    precision = inverse - np.outer(u, u) / u.sum()
    weights = query @ precision + u / u.sum()
    loo = np.eye(n) - precision / np.diag(precision)[:, None]
    return weights, loo


def independent_prediction(kernel, query, y, lam):
    center = np.eye(len(kernel)) - np.ones_like(kernel) / len(kernel)
    centered = center @ kernel @ center
    centered_query = query - kernel.mean(axis=0) - query.mean() + kernel.mean()
    model = KernelRidge(alpha=lam, kernel='precomputed').fit(centered, y-y.mean(axis=0))
    return model.predict(centered_query[None])[0] + y.mean(axis=0)


def fit(data, descriptor):
    train = data['trainingTargets']
    raw, ctrl = counts(data['trainingCounts']), counts(data['control'])
    cpm = raw / raw.sum(axis=1, keepdims=True) * 1e6
    baseline = np.log1p(ctrl / ctrl.sum() * 1e6)
    delta = np.log1p(cpm) - baseline
    mean_cpm = (cpm.sum(axis=0) + np.expm1(baseline)) / (len(train)+1)
    top = np.array(sorted(range(len(ctrl)), key=lambda i: (-mean_cpm[i], data['featureIDs'][i]))[:1000])
    predictions = dict(noChange=baseline.copy(), meanSingleResponse=baseline+delta.mean(axis=0))
    info = dict(target=data['target'], trainingTargets=train, excludedConditions=data['excluded'],
                supported=data['target'] in descriptor['targetIDs'])
    if info['supported']:
        ids = descriptor['targetIDs'].tolist()
        usable = [i for i, t in enumerate(train) if t in ids]
        rows = [ids.index(train[i]) for i in usable]
        kernel = descriptor['kernel'][np.ix_(rows, rows)]
        query = descriptor['kernel'][ids.index(data['target']), rows]
        y = delta[usable]
        predictions['meanSupportedResponse'] = baseline+y.mean(axis=0)
        weights, losses, chosen, errors = {}, {}, {}, []
        for relation, response in [('ordinary', y), ('shuffled', np.roll(y, -1, axis=0))]:
            scores, candidates = [], []
            for lam in LAMBDAS:
                weight, loo = solve_weights(kernel, query, lam)
                cv = np.maximum(baseline+loo@response, 0)-baseline
                loss = float(np.sqrt(np.mean((cv-response)**2, axis=1)).mean())
                scores.append(dict(lambdaValue=lam, meanInnerRMSE=loss))
                candidates.append((loss, -lam, weight, lam))
            _, _, weight, lam = min(candidates, key=lambda v: (v[0],v[1]))
            method = 'goRidgeNested' if relation == 'ordinary' else 'shuffledGoRidgeNested'
            predicted = weight@response
            oracle = independent_prediction(kernel, query, response, lam)
            np.testing.assert_allclose(predicted, oracle, rtol=1e-9, atol=1e-10)
            errors.append(float(np.max(np.abs(predicted-oracle))))
            predictions[method] = baseline+predicted
            weights[method], chosen[method], losses[method] = weight.tolist(), lam, scores
        fixed, _ = solve_weights(kernel, query, 1)
        fixed_response = fixed@y
        oracle = independent_prediction(kernel, query, y, 1)
        np.testing.assert_allclose(fixed_response, oracle, rtol=1e-9, atol=1e-10)
        errors.append(float(np.max(np.abs(fixed_response-oracle))))
        predictions['goRidgeFixed'] = baseline+fixed_response
        weights['goRidgeFixed'] = fixed.tolist()
        info.update(descriptorTrainingTargets=[train[i] for i in usable], weights=weights,
                    selectedLambda=chosen, innerLosses=losses, independentMaxAbsoluteError=max(errors))
    result = {k: np.maximum(v, 0) for k,v in predictions.items()}
    info['diagnostics'] = {k: dict(clippedFeatures=int(np.count_nonzero(predictions[k] < 0)),
                                  impliedCPMSum=float(np.expm1(v).sum())) for k,v in result.items()}
    return result, baseline, top, info


def predict(ref_path, prepared, out):
    ref = reference(ref_path)
    receipt = json.loads((prepared / 'receipt.json').read_text())
    for name,key in [('kernel.npz','kernelSHA256'), ('annotations.json','annotationsSHA256'), ('coverage.json','coverageSHA256')]:
        if sha(prepared/name) != receipt[key]:
            raise ValueError('descriptor fingerprint differs')
    desc = load(prepared / 'kernel.npz')
    targets = sorted(t for t in ref['conditions'] if t != 'control' and '_' not in t)
    out.mkdir(parents=True, exist_ok=False)
    arrays = {m: [] for m in METHODS}
    infos, panels, baselines = [], [], []
    for target in targets:
        data = select(ref, target)
        changed = dict(ref, counts=ref['counts'].copy())
        names = ref['conditions'].tolist()
        changed['counts'][[names.index(n) for n in data['excluded']]] = 1
        mutated = select(changed, target)
        for k in ('control', 'trainingCounts', 'featureIDs', 'symbols'):
            np.testing.assert_array_equal(data[k], mutated[k])
        result, baseline, top, info = fit(data, desc)
        info['excludedOutcomeMutationInputsExact'] = True
        for method in METHODS:
            arrays[method].append(result.get(method, np.full(len(baseline), np.nan)))
        infos.append(info); panels.append(top); baselines.append(baseline)
        print(target+' frozen', flush=True)
    np.savez_compressed(out / 'predictions.npz', **{m:np.array(v) for m,v in arrays.items()},
                        targetIDs=np.array(targets), featureIDs=ref['feature_ids'],
                        baseline=np.array(baselines), trainingTop1000=np.array(panels),
                        supported=np.array([i['supported'] for i in infos]))
    write(out / 'folds.json', infos)
    write(out / 'receipt.json', dict(referenceSHA256=REFERENCE_SHA256, implementation=identity(),
          descriptorReceiptSHA256=sha(prepared/'receipt.json'), predictionsSHA256=sha(out/'predictions.npz'),
          foldsSHA256=sha(out/'folds.json'), targets=len(targets), supported=sum(i['supported'] for i in infos),
          outcomesReadByFitter=False, outerOutcomesUsedForLambda=False,
          independentMaxAbsoluteError=max(i.get('independentMaxAbsoluteError',0) for i in infos)))


def score(ref_path, predictions, out):
    ref = reference(ref_path)
    receipt = json.loads((predictions/'receipt.json').read_text())
    for name,key in [('predictions.npz','predictionsSHA256'), ('folds.json','foldsSHA256')]:
        if sha(predictions/name) != receipt[key]:
            raise ValueError('prediction fingerprint differs')
    pred = load(predictions/'predictions.npz')
    np.testing.assert_array_equal(pred['featureIDs'],ref['feature_ids'])
    names = ref['conditions'].tolist()
    results = []
    for i,target in enumerate(pred['targetIDs']):
        raw = counts(ref['counts'][names.index(target)])
        truth = np.log1p(raw/raw.sum()*1e6)-pred['baseline'][i]
        for method in METHODS:
            if method in METHODS[2:] and not pred['supported'][i]:
                continue
            for panel,ix in [('allGenes',np.arange(len(truth))),('trainingTop1000',pred['trainingTop1000'][i])]:
                results.append(dict(target=str(target),supported=bool(pred['supported'][i]),method=method,panel=panel,
                                    **metrics((pred[method][i]-pred['baseline'][i])[ix],truth[ix])))
    summary=[]
    for cohort in ('allAvailable','matchedSupported'):
        for method in METHODS:
            for panel in ('allGenes','trainingTop1000'):
                rows=[r for r in results if r['method']==method and r['panel']==panel and (cohort=='allAvailable' or r['supported'])]
                zero={r['target']:r['rmse'] for r in results if r['method']=='noChange' and r['panel']==panel}
                summary.append(dict(cohort=cohort,method=method,panel=panel,targets=len(rows),
                                    meanResponseRMSE=float(np.mean([r['rmse'] for r in rows])),
                                    worseThanNoChange=[r['target'] for r in rows if r['rmse']>zero[r['target']]]))
    out.mkdir(parents=True,exist_ok=False)
    write(out/'results.json',results);write(out/'summary.json',summary)
    write(out/'receipt.json',dict(referenceSHA256=REFERENCE_SHA256,predictionReceiptSHA256=sha(predictions/'receipt.json'),implementation=identity()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['fetch','prepare','predict','score'])
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--source',type=Path);p.add_argument('--prepared',type=Path);p.add_argument('--predictions',type=Path)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    with threadpool_limits(limits=1):
        if a.mode=='fetch':fetch(a.reference,a.out)
        elif a.mode=='prepare':
            if a.source is None:p.error('prepare requires --source')
            prepare(a.reference,a.source,a.out)
        elif a.mode=='predict':
            if a.prepared is None:p.error('predict requires --prepared')
            predict(a.reference,a.prepared,a.out)
        else:
            if a.predictions is None:p.error('score requires --predictions')
            score(a.reference,a.predictions,a.out)

if __name__=='__main__':main()
