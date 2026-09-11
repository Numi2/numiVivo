#!/usr/bin/env python3
"""Audit donor-score sensitivity of complete, already frozen HIRISA predictions.

This is retrospective score omission, with no refitting or interval estimation.
The primary endpoint remains response RMSE across all source genes. Both frozen
feature families and every method comparison are retained in the JSON output.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCES = (
    ('known-treatment', '2026-09-10-prediction', 'prediction-scores.json',
     '8d49cd4eb3ce2ced43f006cec931d848903251cfff0d496621004cc3637ce44f', 79),
    ('preparation-transfer', '2026-09-11-context-transfer', 'study/scores.json',
     '659181702f496192c678520e91157b7406a8a22a3197aaa94b843851ae7b4734', 120),
)
METHODS = ('noChange', 'meanResponse', 'medianResponse', 'contextRidge')
FAMILIES = ('all-source-genes', 'training-selected-context-genes')
PAIRS = (('meanResponse', 'noChange'), ('medianResponse', 'noChange'),
         ('contextRidge', 'noChange'), ('contextRidge', 'meanResponse'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_source(archive, name, expected):
    root = ROOT / 'evidence' / archive
    manifest_bytes = (root / 'manifest.json').read_bytes()
    manifest = json.loads(manifest_bytes)
    parts = [r for r in manifest['records'] if r['sourcePath'] == name]
    require(bool(parts), 'Missing source records')
    payload = bytearray()
    for part in sorted(parts, key=lambda r: r.get('sourceOffset', 0)):
        require(part.get('sourceOffset', 0) == len(payload), 'Noncontiguous source')
        path = root / part['storedPath']
        require(path.resolve().is_relative_to(root.resolve()) and not path.is_symlink(),
                'Invalid archived path')
        stored = path.read_bytes()
        require((len(stored), sha(stored)) ==
                (part['storedBytes'], part['storedSHA256']), 'Stored hash mismatch')
        decoded = gzip.decompress(stored) if part['gzipEncoded'] else stored
        require((len(decoded), sha(decoded)) ==
                (part['sourceBytes'], part['sourceSHA256']), 'Decoded hash mismatch')
        payload.extend(decoded)
    require(sha(payload) == expected, 'Frozen score source changed')
    return json.loads(payload), dict(archive=archive, sourcePath=name,
                                   SHA256=expected, bytes=len(payload),
                                   manifestSHA256=sha(manifest_bytes))


def comparison(values, candidate, baseline):
    """values maps donor to a paired difference: candidate error minus baseline."""
    donors = sorted(values)
    require(len(donors) >= 3, 'Insufficient donors for omission diagnostic')
    require(all(math.isfinite(values[d]) for d in donors), 'Nonfinite difference')
    mean = math.fsum(values.values()) / len(donors)
    omitted = {d: math.fsum(values[x] for x in donors if x != d) / (len(donors) - 1)
               for d in donors}
    low, high = min(omitted.values()), max(omitted.values())
    status = ('stable-improvement' if mean < 0 and high < 0 else
              'stable-worse' if mean > 0 and low > 0 else
              'all-tied' if mean == low == high == 0 else 'donor-sensitive')
    return dict(candidate=candidate, baseline=baseline, donorCount=len(donors),
                pairedDifferences=values, meanDifference=mean,
                donorWins=sum(v < 0 for v in values.values()),
                donorTies=sum(v == 0 for v in values.values()),
                omittedDonorMeans=omitted, omittedMinimum=low, omittedMaximum=high,
                status=status)


def audit():
    results, identities, transfer = [], [], {}
    total_folds = 0
    for study, archive, name, digest, count in SOURCES:
        data, identity = read_source(archive, name, digest)
        identities.append(identity)
        require(data['status'] == 'passed', 'Upstream scoring failed')
        require(data['completedFolds'] == data['foldCount'] == len(data['folds']) == count,
                'Incomplete fold inventory')
        require(len({f['id'] for f in data['folds']}) == count, 'Duplicate fold')
        groups = {}
        for fold in data['folds']:
            key = ((fold['population'], fold['treatment']) if study == 'known-treatment'
                   else (fold['kind'], fold['queryPreparation'], fold['lineage']))
            scores = {(s['family'], s['method']): s for s in fold['scores']}
            require(len(scores) == len(fold['scores']) == 8, 'Duplicate/missing score')
            require(set(scores) == {(f, m) for f in FAMILIES for m in METHODS},
                    'Unexpected method or feature family')
            for (family, method), score in scores.items():
                require(math.isfinite(score['responseRMSE']) and score['responseRMSE'] >= 0,
                        'Invalid response RMSE')
                if family == FAMILIES[0]:
                    require(score['features'] == 18082, 'Incomplete primary feature axis')
            group = groups.setdefault(key, {})
            require(fold['heldOutDonor'] not in group, 'Duplicate donor within contrast')
            group[fold['heldOutDonor']] = scores
        require(len(groups) == (16 if count == 79 else 24), 'Missing contrasts')
        for key, donors in sorted(groups.items()):
            expected_count = 4 if study == 'known-treatment' and key == ('Bcell', 'IFNg') else 5
            require(len(donors) == expected_count, 'Wrong contrast donor coverage')
            for family in FAMILIES:
                # Reconcile the source's published equal-donor summary, independently
                # of the paired-difference calculation used below.
                for method in METHODS:
                    summaries = [s for s in data['summaries'] if s['family'] == family
                                 and s['method'] == method and
                                 ((s['population'], s['treatment']) if study == 'known-treatment'
                                  else (s['kind'], s['queryPreparation'], s['lineage'])) == key]
                    require(len(summaries) == 1, 'Ambiguous archived contrast summary')
                    field = 'equalFoldMeans' if study == 'known-treatment' else 'equalDonorMeans'
                    observed = math.fsum(v[family, method]['responseRMSE']
                                         for v in donors.values()) / len(donors)
                    require(math.isclose(observed, summaries[0][field]['responseRMSE'],
                                         rel_tol=0, abs_tol=1e-12), 'Summary does not reconstruct')
                for candidate, baseline in PAIRS:
                    diff = {d: s[family, candidate]['responseRMSE'] - s[family, baseline]['responseRMSE']
                            for d, s in donors.items()}
                    results.append(dict(study=study, contrast=list(key), family=family,
                                        **comparison(diff, candidate, baseline)))
        if study == 'preparation-transfer':
            transfer = groups
        total_folds += count
    # Across-preparation versus within-preparation comparisons use identical query
    # donors and the common full-gene axis; training-selected panels are not paired.
    penalties = []
    for key, cross in sorted(transfer.items()):
        if key[0] != 'cross':
            continue
        within = transfer[('within', *key[1:])]
        require(set(cross) == set(within), 'Unmatched transfer donors')
        for method in METHODS:
            diff = {d: cross[d][FAMILIES[0], method]['responseRMSE'] -
                    within[d][FAMILIES[0], method]['responseRMSE'] for d in cross}
            penalties.append(dict(queryPreparation=key[1], lineage=key[2],
                                  **comparison(diff, 'cross-' + method, 'within-' + method)))
    return dict(schemaVersion=1, status='passed', sourceFolds=total_folds,
                inputs=identities, comparisons=results, transferPenalties=penalties,
                scope='Retrospective donor-score omission; no refitting, new predictions, '
                      'confidence intervals, independent-study validation or model promotion',
                scoringUnit='Paired held-out donor RMSE difference in natural-log(1+CPM) units',
                negativeDifference='Candidate has lower error',
                independence='Overlapping training donors and reused contexts prevent treating '
                             'folds or contrasts as independent experiments')


def render(data):
    primary = [r for r in data['comparisons'] if r['family'] == FAMILIES[0]
               and r['candidate'] == 'contextRidge']
    summary = []
    for study, kind, label in [('known-treatment', None, 'Known treatment'),
                               ('preparation-transfer', 'cross', 'Cross preparation'),
                               ('preparation-transfer', 'within', 'Matched within preparation')]:
        for base, base_label in [('noChange', 'No change'), ('meanResponse', 'Mean response')]:
            selected = [r for r in primary if r['study'] == study and r['baseline'] == base
                        and (kind is None or r['contrast'][0] == kind)]
            summary.append('| ' + ' | '.join([label, base_label,
                str(sum(r['meanDifference'] < 0 for r in selected)) + '/' + str(len(selected)),
                str(sum(r['status'] == 'stable-improvement' for r in selected)),
                str(sum(r['status'] == 'stable-worse' for r in selected)),
                str(sum(r['status'] == 'donor-sensitive' for r in selected)),
                str(sum(r['donorWins'] for r in selected)) + '/' +
                str(sum(r['donorCount'] for r in selected))]) + ' |')
    lines = ['# Sensitivity of completed prediction scores to individual donors', '',
             'This retrospective diagnostic uses all **199 completed HIRISA folds**: '
             '79 known-treatment folds and 120 cross/within-preparation folds. It '
             'reconstructs the published summaries from hash-verified score archives. '
             'No counts, models, folds or predictions change.', '',
             'For each comparison, subtract baseline RMSE from candidate RMSE for '
             'each donor, then omit each donor score once and average the remaining '
             'differences. A negative difference favors the candidate. A stable '
             'improvement stays strictly negative in the full mean and every omission; '
             'stable worse stays strictly positive. Otherwise the result is donor-sensitive '
             '(or all-tied). These are score-aggregation sensitivities, not retrained '
             'leave-one-out experiments or confidence intervals.', '',
             '| Ridge comparison setting | Baseline | Full-mean wins | Stable improvement | Stable worse | Donor-sensitive | Individual donor wins |',
             '| --- | --- | ---: | ---: | ---: | ---: | ---: |', *summary, '',
             'For known treatment, the three stable ridge improvements over the mean '
             'are Monocyte IFNa, IFNb and IFNg. Its additional full-mean win, NK IFNa, '
             'changes sign under a donor omission. NK IFN-L1 is similarly sensitive '
             'against no-change. All twelve cross-preparation mean wins over no-change '
             'survive every omission, although four individual ridge predictions are '
             'worse. All twelve cross-versus-within ridge penalties also survive. '
             'This supports conditional average-response signal and a limited benefit '
             'from donor-context ridge, with the original failures retained.', '',
             'All rows below use **18,082 source genes** and equally weighted donors. '
             'The JSON retains both original feature families, mean/median/ridge '
             'comparisons, every paired difference and every omitted-donor result.', '',
             '| Setting / contrast | Donors | Ridge − no change | Omission range | Status | Ridge − mean | Omission range | Status |',
             '| --- | ---: | ---: | --- | --- | ---: | --- | --- |']
    rows = [r for r in data['comparisons'] if r['family'] == FAMILIES[0]
            and r['candidate'] == 'contextRidge']
    keys = sorted({(r['study'], tuple(r['contrast'])) for r in rows})
    for study, contrast in keys:
        pair = {r['baseline']: r for r in rows if r['study'] == study and tuple(r['contrast']) == contrast}
        first = pair['noChange']
        cells = [' / '.join(contrast), str(first['donorCount'])]
        for base in ('noChange', 'meanResponse'):
            row = pair[base]
            cells += [f'{row["meanDifference"]:+.6f}',
                      f'[{row["omittedMinimum"]:+.6f}, {row["omittedMaximum"]:+.6f}]', row['status']]
        lines.append('| ' + ' | '.join(cells) + ' |')
    lines += ['', '## Transfer penalty', '',
              '| Query preparation / lineage | Cross ridge − within ridge | Omission range | Status |',
              '| --- | ---: | --- | --- |']
    for row in data['transferPenalties']:
        if row['candidate'] == 'cross-contextRidge':
            lines.append(f'| {row["queryPreparation"]} / {row["lineage"]} | '
                         f'{row["meanDifference"]:+.6f} | '
                         f'[{row["omittedMinimum"]:+.6f}, {row["omittedMaximum"]:+.6f}] | {row["status"]} |')
    lines += ['', '## Interpretation and reproduction', '',
              'A stable mean improvement can coexist with individual donor failures. '
              'Models share training donors; lineages and treatments reuse donors. '
              'Do not count 199 folds, genes or cells as independent biological replicates. '
              'Omitting a scoring donor does not remove that donor from other fitted '
              'models. This post hoc diagnostic cannot establish predictive interval '
              'coverage, repair preparation/batch confounding or validate a new study.', '',
              'The [original known-treatment results](PREDICTION_RESULTS.md), '
              '[preparation-transfer results](CONTEXT_TRANSFER.md) and '
              '[biological assessment](../../../../Documentation/BiologicalPrediction.md) '
              'retain their scientific limits and negative findings.', '',
              'From the repository root, using standard-library Python and new output paths:', '',
              '```sh',
              'python3 Tools/Omics/Benchmarks/HIRISA/prediction_sensitivity.py --out /tmp/hirisa-sensitivity.json --report /tmp/hirisa-sensitivity.md',
              '```', '',
              'The [evidence manifest](evidence/2026-09-11-prediction-sensitivity/manifest.json) '
              'binds the complete diagnostic and independent arithmetic check. '
              'Source archive identities are embedded in the diagnostic.']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    require(not args.out.exists() and not args.report.exists(), 'Outputs must be new')
    data = audit()
    args.out.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n')
    args.report.write_text(render(data))
    print(json.dumps(dict(status=data['status'], folds=data['sourceFolds'],
                          comparisons=len(data['comparisons']), transferPenalties=len(data['transferPenalties']))))


if __name__ == '__main__':
    main()
