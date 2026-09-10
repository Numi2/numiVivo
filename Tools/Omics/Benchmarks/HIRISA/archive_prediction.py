#!/usr/bin/env python3
"""Archive HIRISA frozen prediction outputs, all scores and retained failed attempts."""
import argparse,gzip,json,shutil
from pathlib import Path
from acquire import digest

def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();root=a.root;out=a.output;scores=read(root/'prediction-scores.json')
 assert scores['foldCount']==79 and len(scores['summaries'])==128
 assert read(root/'prediction-batch-publish-status.json')['returnCode']==0
 assert read(root/'prediction-batch-verify-status.json')['returnCode']==0
 assert '7 tests in 2 suites passed' in (root/'aggregate-prediction-tests.log').read_text()
 for name in ('prediction-product-check','prediction-product-check-release'):
  r=read(root/name/'checks.json');assert r['status']=='passed' and len(r['checks'])==10
 freeze=read(root/'prediction-native-output-freeze.json')
 for name,expected in freeze['files'].items():assert digest(root/'prediction-batch'/name)==expected
 names=['prediction-folds.json','prediction-batch-plan.json','prediction-scoring-targets.json','prediction-transport-receipt.json',
  'prediction-execution-environment.json','prediction-batch-publish-status.json','prediction-batch-verify-status.json',
  'prediction-batch-publish.log','prediction-batch-verify.log','prediction-batch-driver.log','prediction-native-output-freeze.json',
  'prediction-output-freeze.json','prediction-scores.json','prediction-score.log','prepare-prediction-batch.log',
  'prepare-prediction-batch-attempt-1.log','prepare-prediction-batch-attempt-1.py','prediction-release-build.log',
  'prediction-runtime-debug-manifest.json','prediction-runtime-release-manifest.json',
  'aggregate-prediction-tests.log','aggregate-prediction-tests-attempt-1.log',
  'prediction-product-check.log','prediction-product-check-release.log','run_prediction_batch.py',
  'remote_prediction_product.py','remote_prediction_product_release.py','remote_product_release.py',
  'native-independent-verification.json','prediction-host.json']
 for folder in ('prediction-batch','prediction-product-check','prediction-product-check-release','prediction-transport-attempt-1'):
  names += [str(path.relative_to(root)) for path in sorted((root/folder).rglob('*')) if path.is_file()
            and str(path.relative_to(root)) not in ('prediction-batch/source/original.h5ad','prediction-batch/source/report.json')]
 names=sorted(set(names));assert all((root/name).is_file() and not (root/name).is_symlink() for name in names)
 out.mkdir(parents=True,exist_ok=False);records=[]
 for name in names:
  source=root/name;encoded=source.suffix not in ('.gz','.npz','.h5ad');destination=out/(name+('.gz' if encoded else ''))
  destination.parent.mkdir(parents=True,exist_ok=True);before=digest(source)
  with source.open('rb') as src,destination.open('xb') as dst:
   if encoded:
    with gzip.GzipFile(fileobj=dst,mode='wb',filename='',mtime=0) as zipped:shutil.copyfileobj(src,zipped,1048576)
   else:shutil.copyfileobj(src,dst,1048576)
  assert digest(source)==before
  records.append(dict(sourcePath=name,storedPath=str(destination.relative_to(out)),sourceBytes=source.stat().st_size,
      sourceSHA256=before,storedBytes=destination.stat().st_size,storedSHA256=digest(destination),gzipEncoded=encoded))
 manifest=dict(schemaVersion=1,records=records,folds=79,completedFolds=scores['completedFolds'],numericalAndScoringStatus=scores['status'],
  sourceSHA256=read(root/'native-independent-verification.json')['sourceSHA256'],archiveScriptSHA256=digest(Path(__file__)),
  sourceIngestionArchive='2026-09-10-native',nativeReplayPassed=True,
  scope='All frozen known-perturbation/new-donor predictions, independent numerical reconstruction and empirical errors. No unseen-perturbation identity, uncertainty or biological calibration. Exact original source and native executables remain external; the full source report and library count references are in the ingestion archive.')
 (out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
 print(json.dumps(dict(members=len(records),storedBytes=sum(r['storedBytes'] for r in records),manifestSHA256=digest(out/'manifest.json'))))

if __name__=='__main__':main()
