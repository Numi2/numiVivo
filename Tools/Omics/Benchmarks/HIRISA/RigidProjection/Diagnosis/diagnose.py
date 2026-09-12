"""Recount retained annotation failures; no model fitting or new qualification."""
from pathlib import Path
import collections
import hashlib
import json
import tarfile

root = Path(__file__).resolve().parent
study = root.parent
manifest = json.loads((study / 'manifest.json').read_text())
archive = study / 'evidence.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == manifest['archiveSHA256']
with tarfile.open(archive) as source:
    payload = source.extractfile('evaluation/rigid.json').read()
assert hashlib.sha256(payload).hexdigest() == manifest['members']['evaluation/rigid.json']
comparisons = json.loads(payload)['comparisons']
by_label = collections.defaultdict(list)
for item in comparisons:
    by_label[item['labelName']].append(item)
rows = []
for label, items in sorted(by_label.items()):
    eligible = [x for x in items if x['sufficientSupport'] and x['controlSensitive']]
    failures = [x for x in eligible if x['marginFailure']]
    rows.append(dict(label=label, comparisons=len(items), eligible=len(eligible),
                     insufficientSupport=sum(not x['sufficientSupport'] for x in items),
                     supportedButInsensitive=sum(x['sufficientSupport'] and not x['controlSensitive'] for x in items),
                     failures=len(failures), rareFailures=sum(x['rare'] for x in failures),
                     failedStrata=[x['stratum'] for x in failures],
                     maximumMeanRecallLoss=max((x['meanRecallLoss'] for x in eligible), default=None),
                     maximumFoldRecallLoss=max((x['maximumFoldRecallLoss'] for x in eligible), default=None)))
assert sum(x['comparisons'] for x in rows) == 471
assert sum(x['eligible'] for x in rows) == 146
assert sum(x['failures'] for x in rows) == 31
assert sum(x['rareFailures'] for x in rows) == 9
result = dict(schemaVersion=1, scope='Retrospective accounting of existing full-cohort annotation evaluation; no refit',
              archiveSHA256=manifest['archiveSHA256'], evaluationSHA256=hashlib.sha256(payload).hexdigest(),
              comparisons=471, eligible=146, failures=31, rareFailures=9, labels=rows)
print(json.dumps(result, indent=2, sort_keys=True))
