#!/usr/bin/env python3
"""Render all frozen HIRISA contrast results without selecting a winning method."""
import argparse,hashlib,json,re
from pathlib import Path

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();root=a.root;source=root/'prediction-scores.json';data=json.loads(source.read_text())
 assert data['status']=='passed' and data['completedFolds']==data['foldCount']==79
 pub=json.loads((root/'prediction-batch-publish-status.json').read_text());verify=json.loads((root/'prediction-batch-verify-status.json').read_text())
 assert pub['returnCode']==verify['returnCode']==0
 methods=['noChange','meanResponse','medianResponse','contextRidge'];cohorts=sorted({s['cohort'] for s in data['summaries']})
 def rows(family,c):return {s['method']:s for s in data['summaries'] if s['family']==family and s['cohort']==c}
 wins={m:sum(rows('all-source-genes',c)[m]['equalFoldMeans']['responseRMSE']<rows('all-source-genes',c)['noChange']['equalFoldMeans']['responseRMSE'] for c in cohorts) for m in methods[1:]}
 ridge_mean=sum(rows('all-source-genes',c)['contextRidge']['equalFoldMeans']['responseRMSE']<rows('all-source-genes',c)['meanResponse']['equalFoldMeans']['responseRMSE'] for c in cohorts)
 error=max(f['numericalMaximumAbsoluteErrors']['contextRidge-treated'] for f in data['folds'])
 rss={name:int(re.search(r'(\d+)\s+maximum resident set size',(root/f'prediction-batch-{name}.log').read_text()).group(1)) for name in ('publish','verify')}
 lines=['# HIRISA frozen donor-held-out prediction results','',
  'All 79 prespecified folds completed with native source reconstruction and byte-exact replay. Every model and query output was hashed before opening held-out treated counts for scoring. Independent NumPy reconstruction checks all selected features, control moments, coefficients and all four prediction vectors.',
  '',f'The maximum ridge treated-expression difference was {error:.3e}. Publication (including one full-source reconstruction and all fits) took {pub["seconds"]:.2f} seconds; replay took {verify["seconds"]:.2f} seconds on the physical M4 Pro Mac mini. Maximum RSS was {rss["publish"]:,}/{rss["verify"]:,} bytes. These are observed CPU timings, not a controlled scverse speed comparison or a Metal claim.',
  '',f'Mean, median and ridge responses improve the full-gene contrast RMSE over no-change in {wins["meanResponse"]}, {wins["medianResponse"]} and {wins["contextRidge"]} of 16 contrasts respectively. Ridge improves on the training-mean response in only {ridge_mean} of 16. Monocyte IFN-L1 and NK IFNg retain worse errors than no-change for every learned baseline. No method is selected or promoted from these outcomes.',
  '', 'The tables report equally weighted held-out donor means within each contrast, in natural-log(1+CPM) units. Training has four donor pairs, except Bcell IFNg which has three. All 18,082 source genes and the separately frozen training-selected context family are shown. Both treated-expression and response RMSE are retained in the full results; subtracting the same control makes them equal apart from rounding. All fold-level MAE, response Pearson correlations and implied CPM sums remain in the archived JSON. Constant-response correlation, including no-change, remains null.', '']
 for family,title in [('all-source-genes','All 18,082 source genes'),('training-selected-context-genes','Training-selected context genes')]:
  lines += ['## '+title,'','| Population | Treatment | Donor folds | No change | Mean response | Median response | Context ridge |','| --- | --- | ---: | ---: | ---: | ---: | ---: |']
  for c in cohorts:
   s=rows(family,c);v=s['noChange'];assert all(t['status']=='complete' for t in s.values())
   lines.append('| '+ ' | '.join([v['population'],v['treatment'],str(v['expectedFolds']),*[f'{s[m]["equalFoldMeans"]["responseRMSE"]:.6f}' for m in methods]])+' |')
  lines.append('')
 lines += ['## Evidence and limits','',
  'The [prediction archive](evidence/2026-09-10-prediction) retains every native model/prediction, source and implementation bindings, original protocols/folds, both output freezes, numerical checks, all scores, product tests and failed attempts. The first transport preparation missed the pool component of the batch identity and was repaired before publication. A later transport description contained the wrong GEO accession; that native attempt was explicitly interrupted during source reconstruction before fold fitting. Its records remain alongside the corrected GSE306664 run. Neither repair changed the frozen scientific protocol, folds or numerical model.',
  '', 'This is empirical prediction of a known perturbation in an unseen donor from its observed control, using experimental enrichment populations rather than verified single-cell annotations. The separate PBMC transfer split is still unfrozen. The results do not establish unseen perturbation identity, calibrated uncertainty, causal mechanisms or general biological validity. Metadata/report residency remains a limitation. Subsequent complete PCA/graph and seed-7 integration results are linked from the benchmark README; full-cohort clustering and broader biological preservation remain open.',
  '', 'Score JSON SHA256: `'+hashlib.sha256(source.read_bytes()).hexdigest()+'`.',
  '', 'Reproduce the numerical checks and scores with `python score_prediction_batch.py --root /absolute/path/hirisa`, then render this report with `python report_prediction.py --root /absolute/path/hirisa --output PREDICTION_RESULTS.md`. Restore the complete native source/report and frozen bundle according to the archive links; verification of a historical native receipt requires its original executable identity.']
 a.output.write_text('\n'.join(lines)+'\n')
 print(json.dumps({'contrasts':len(cohorts),'folds':79,'ridgeBetterThanMeanContrasts':ridge_mean,'output':str(a.output)}))

if __name__=='__main__':main()
