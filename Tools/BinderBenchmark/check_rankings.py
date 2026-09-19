#!/usr/bin/env python3
"""Independent standard-library checks of native fixed rankings and outcome joins."""
from pathlib import Path
import csv
import json
import math
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
FEATURES = ['ipsae_min_boltz2', 'ipsae_min_ef2fast', 'ipsae_min_ef2full', 'ipsae_min_ptxv2']


def close(a, b):
    if not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10):
        raise AssertionError((a, b))


def command(binary: Path, *args, expected=0):
    result = subprocess.run([str(binary), *map(str, args)], text=True, capture_output=True)
    if result.returncode != expected:
        raise AssertionError((args, result.returncode, result.stdout, result.stderr))
    return result


def verify_ranking(query: dict, plan: dict, ranking: dict):
    assert set(query) == {'schemaVersion', 'sourceSHA256', 'candidates'}
    assert all(set(c) == {'id', 'target', 'scores'} for c in query['candidates'])
    targets = sorted({c['target'] for c in query['candidates']})
    assert [t['target'] for t in ranking['targets']] == targets
    features = [c['feature'] for c in plan['components']]
    total_weight = math.fsum(c['weight'] for c in plan['components'])
    for result in ranking['targets']:
        pool = [c for c in query['candidates'] if c['target'] == result['target']]
        complete = [c for c in pool if set(features) <= set(c['scores'])]
        excluded = {c['id']: sorted(set(features) - set(c['scores'])) for c in pool if c not in complete}
        assert result['inputCount'] == len(pool) and result['excluded'] == excluded
        expected = {}
        moments = {}
        if complete and plan['method'] == 'targetStandardizedMean':
            for feature in features:
                xs = [c['scores'][feature] for c in complete]
                mean = math.fsum(xs) / len(xs)
                sd = math.sqrt(math.fsum((x - mean) ** 2 for x in xs) / len(xs))
                moments[feature] = mean, sd
            for norm in result['normalization']:
                mean, sd = moments[norm['feature']]
                close(norm['mean'], mean); close(norm['populationSD'], sd)
        for c in complete:
            if plan['method'] == 'fixedScore':
                expected[c['id']] = c['scores'][features[0]]
            else:
                terms = []
                for component in plan['components']:
                    feature = component['feature']; mean, sd = moments[feature]
                    terms.append(component['weight'] * ((c['scores'][feature] - mean) / sd if sd else 0))
                expected[c['id']] = math.fsum(terms) / total_weight
        rows = result['candidates']
        assert {r['id'] for r in rows} == set(expected)
        assert rows == sorted(rows, key=lambda r: (-r['score'], r['id']))
        for row in rows:
            close(row['score'], expected[row['id']])
        k = min(plan['topK'], len(rows)); assert result['effectiveK'] == k
        if k:
            # Tie arithmetic is checked against the actual score words, after independently
            # checking the scores: arithmetic rounding must not redefine exact-score ties.
            cutoff = sorted([r['score'] for r in rows], reverse=True)[k-1]
            above = sum(r['score'] > cutoff for r in rows)
            equal = sum(r['score'] == cutoff for r in rows)
            for r in rows:
                p = 1 if r['score'] > cutoff else (k-above)/equal if r['score'] == cutoff else 0
                close(r['selectionWeight'], p)
        close(math.fsum(r['selectionWeight'] for r in rows), k)


def verify_assessment(records: list[dict], ranking: dict, assessment: dict):
    lookup = {r['id']: r for r in records}
    assert [r['target'] for r in ranking['targets']] == [r['target'] for r in assessment['targets']]
    for ranked, result in zip(ranking['targets'], assessment['targets']):
        weights = {}
        pairs = []
        for candidate in ranked['candidates']:
            row = lookup[candidate['id']]
            assert row['target'] == result['target']
            weights[row['outcome']] = weights.get(row['outcome'], 0) + candidate['selectionWeight']
            if row['outcome'] in ('binder', 'nonBinder'):
                pairs.append((int(row['outcome'] == 'binder'), candidate['score']))
        assert set(weights) == set(result['selectedWeightByOutcome'])
        for key in weights: close(weights[key], result['selectedWeightByOutcome'][key])
        close(result['measuredBinderHits'], weights.get('binder', 0))
        close(result['selectedNonbinaryWeight'], sum(v for k,v in weights.items() if k not in ('binder','nonBinder')))
        pos = [s for y,s in pairs if y]; neg = [s for y,s in pairs if not y]
        assert result['binaryOutcomeCount'] == len(pairs)
        if pos and neg:
            close(result['binaryAUROC'], sum((a>b) + .5*(a==b) for a in pos for b in neg)/(len(pos)*len(neg)))
        else: assert result.get('binaryAUROC') is None
        if pos:
            ap = 0
            for threshold in set(s for _,s in pairs):
                newpos = sum(y for y,s in pairs if s == threshold)
                selected = [y for y,s in pairs if s >= threshold]
                ap += newpos/len(pos) * sum(selected)/len(selected)
            close(result['binaryAveragePrecision'], ap)
        else: assert result.get('binaryAveragePrecision') is None


def check(binary: Path, root: Path):
    source = root/'source.csv'; config = root/'import.json'
    with source.open('w', newline='') as h:
        w = csv.writer(h); w.writerow(['uuid','target','sequence','adaptyv_binding'] + FEATURES)
        outcomes = ['not_tested', 'binder', 'non_binder', 'expression_failure', 'binder', 'inconclusive']
        for i in range(12):
            # Includes tied values and the same top untested candidate in both targets.
            w.writerow([f'c{i:02}', 'T' if i<6 else 'U', 'ACDEFGHIKLMN'[i]*4, outcomes[i%6]] + [1-(i%6)/10, (i%4)/4, .5, (i%3)/3])
    config.write_text(json.dumps({'schemaVersion':1,'assay':'adaptyv','targets':['T','U'],'features':FEATURES}))
    imported = root/'input'; query = root/'query'
    command(binary,'binder-import',source,config,imported)
    command(binary,'binder-ranking-query',imported,query)
    q = json.loads((query/'query.json').read_text())
    records = json.loads((imported/'imported.json').read_text())['dataset']['records']
    for name in ('boltz2','confidence3'):
        p = ROOT/f'Tools/BinderBenchmark/plans/ranking-{name}.json'
        ranked = root/f'rank-{name}'; assessed = root/f'assess-{name}'
        command(binary,'binder-rank',query,p,ranked)
        command(binary,'binder-assess-ranking',imported,ranked,assessed)
        for path in (query,ranked,assessed): command(binary,'binder-verify',path)
        r = json.loads((ranked/'ranking.json').read_text())
        verify_ranking(q,json.loads(p.read_text()),r)
        verify_assessment(records,r,json.loads((assessed/'assessment.json').read_text()))
    badplan = root/'bad-plan.json'
    badplan.write_text(json.dumps({'schemaVersion':1,'method':'fixedScore','components':[{'feature':FEATURES[0],'weight':1}],'topK':1,'outcomes':[]}))
    command(binary,'binder-rank',query,badplan,root/'bad',expected=65)
    assert not (root/'bad').exists()
    command(binary,'binder-rank',query,expected=64)
    print('PASS: outcome-blind CLI query/rank/assess/replay, independent scores, ties, missing outcomes and strict plan rejection')


if __name__ == '__main__':
    import sys
    with tempfile.TemporaryDirectory(prefix='binder-ranking-check-') as tmp:
        check(Path(sys.argv[1]).resolve(), Path(tmp))
