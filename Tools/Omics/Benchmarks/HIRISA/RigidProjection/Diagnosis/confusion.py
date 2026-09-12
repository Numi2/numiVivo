"""Hash-bound marginal confusion changes; not paired cell transitions or refitting."""
from pathlib import Path
import collections
import hashlib
import json
import math
import tarfile

root = Path(__file__).resolve().parent.parent
manifest = json.loads((root / 'manifest.json').read_text())
archive = root / 'evidence.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archiveSHA256']
with tarfile.open(archive) as source:
    def read(name):
        payload = source.extractfile(name).read()
        assert hashlib.sha256(payload).hexdigest() == manifest['members'][name]
        return json.loads(payload)
    baseline = read('evaluation/baseline.json')
    candidate = read('evaluation/rigid.json')
assert baseline['matrix'] == 'baseline' and candidate['matrix'] == 'rigid'
assert baseline['freezeSHA256'] == candidate['freezeSHA256']
labels = {}
for c in candidate['comparisons']:
    assert labels.setdefault(c['label'], c['labelName']) == c['labelName']
assert set(labels) == set(range(31))

def index(result):
    out = {}
    for f in result['folds']:
        key = (f['stratum'], f['donor'])
        assert key not in out
        out[key] = f
        matrix = f['confusion']
        assert len(matrix) == len(labels)
        for i, row in enumerate(matrix):
            assert len(row) == len(labels)
            assert all(math.isfinite(v) and v >= 0 and int(v) == v for v in row)
            assert sum(row) == f['queryCounts'][i]
            if sum(row):
                assert abs(row[i] / sum(row) - f['metrics'][i]['recall']) < 1e-12
        assert sum(map(sum, matrix)) == f['cells']
    return out
b, c = index(baseline), index(candidate)
assert b.keys() == c.keys() and len(b) == 114
for key in b:
    assert b[key]['id'] == c[key]['id']
    assert b[key]['queryCounts'] == c[key]['queryCounts']
    assert b[key]['trainingCounts'] == c[key]['trainingCounts']
rows = []
for item in candidate['comparisons']:
    label = item['label']
    fold_rows = []
    for f in item['folds']:
        key = (item['stratum'], f['donor'])
        before, after = b[key]['confusion'][label], c[key]['confusion'][label]
        assert sum(before) == f['queryCells']
        assert abs(before[label] / sum(before) - f['baselineRecall']) < 1e-12
        assert abs(after[label] / sum(after) - f['candidateRecall']) < 1e-12
        changes = [int(y-x) for x,y in zip(before, after)]
        assert sum(changes) == 0
        fold_rows.append(dict(foldID=c[key]['id'], donor=f['donor'], sufficient=f['sufficient'],
            queryCells=f['queryCells'], baselineCorrect=int(before[label]), candidateCorrect=int(after[label]),
            destinations=[dict(label=labels[j], baseline=int(before[j]), candidate=int(after[j]), delta=changes[j])
                          for j in range(len(labels)) if j != label and (before[j] or after[j])]))
    rows.append(dict(label=item['labelName'], stratum=item['stratum'], rare=item['rare'],
        sufficientSupport=item['sufficientSupport'], controlSensitive=item['controlSensitive'],
        marginFailure=item['marginFailure'], folds=fold_rows))
assert len(rows) == 471
summaries = []
for label in labels.values():
    eligible = [x for x in rows if x['label'] == label and x['sufficientSupport'] and x['controlSensitive']]
    fs = [f for x in eligible for f in x['folds'] if f['sufficient']]
    changes = collections.Counter()
    for f in fs:
        for d in f['destinations']:
            changes[d['label']] += d['delta']
    correct_delta = sum(f['candidateCorrect'] - f['baselineCorrect'] for f in fs)
    assert correct_delta + sum(changes.values()) == 0
    summaries.append(dict(label=label, eligibleComparisons=len(eligible), sufficientFolds=len(fs),
        queryAppearances=sum(f['queryCells'] for f in fs),
        baselineCorrect=sum(f['baselineCorrect'] for f in fs),
        candidateCorrect=sum(f['candidateCorrect'] for f in fs),
        destinationDeltas=dict(sorted(changes.items())),
        donorRecallLosses={str(donor):sum(f['candidateCorrect'] < f['baselineCorrect'] for f in fs if f['donor'] == donor)
                           for donor in sorted({f['donor'] for f in fs})}))
result = dict(schemaVersion=1, eligibleLabelSummaries=summaries, scope='Marginal destination counts, not paired cell transitions; all retained comparisons',
    archiveSHA256=manifest['archiveSHA256'],
    members={n:manifest['members'][n] for n in ['evaluation/baseline.json','evaluation/rigid.json']},
    foldCount=114, comparisons=rows)
print(json.dumps(result, separators=(',', ':'), sort_keys=True))
