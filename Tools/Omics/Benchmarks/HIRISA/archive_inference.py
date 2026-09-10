#!/usr/bin/env python3
"""Archive all frozen HIRISA native and reference DE families with diagnostics."""
import argparse,gzip,json,shutil
from pathlib import Path
from acquire import digest

def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();root=a.root;inputs=root/'inference-inputs';out=a.output
 summary=read(inputs/'inference-summary.json');checks=read(inputs/'native-model-check-summary.json')
 assert summary['allCasesTerminal'] and summary['nativeRawIngestionLinked'] and not summary['missingCases']
 assert len(summary['methods'])==96 and len(summary['comparisons'])==144 and len(summary['sixMethodJointFamilies'])==16
 assert checks['status']=='passed' and checks['cases']==48
 assert read(root/'native-full-verify-status.json')['returnCode']==0
 names=['native-full-verify-status.json','native-full-verify.log','native-inference-driver.log','native-independent-verification.json',
        'model-check-final.log','inference-summary-final.log','inference-inputs/receipt.json','inference-inputs/counts.npz',
        'inference-inputs/inference-summary.json','inference-inputs/native-model-check-summary.json',
        'inference-inputs/native/environment.json','inference-inputs/native/driver-complete.json']
 requests=read(inputs/'requests/manifest.json');assert len(requests)==48
 for case in requests:
  d=inputs/'native'/(case['cohort']+'-'+case['method']);r=read(d/'receipt.json')
  assert r['returnCode']==0 and r['outputSHA256']==digest(d/'output.json.gz')
  c=read(d/'model-check.json');assert c['status']=='passed' and c['outputSHA256']==r['outputSHA256']
  names += [str((d/name).relative_to(root)) for name in ('receipt.json','output.json.gz','stderr.log','model-check.json','model-gene-checks.json.gz')]
 for cohort in range(1,17):
  r=read(inputs/'reference'/f'{cohort:02d}-receipt.json');assert r['returnCode']==0
  for name,expected in r['files'].items():assert digest(inputs/'reference'/f'{cohort:02d}'/name)==expected
 for folder in ('cohorts','requests','reference-inputs','reference'):
  names += [str(path.relative_to(root)) for path in sorted((inputs/folder).rglob('*')) if path.is_file()]
 names=sorted(set(names));assert all((root/name).is_file() and not (root/name).is_symlink() for name in names)
 out.mkdir(parents=True,exist_ok=False);records=[]
 for name in names:
  source=root/name;encoded=source.suffix not in ('.gz','.npz');destination=out/(name+('.gz' if encoded else ''))
  destination.parent.mkdir(parents=True,exist_ok=True);before=digest(source)
  with source.open('rb') as src,destination.open('xb') as dst:
   if encoded:
    with gzip.GzipFile(fileobj=dst,mode='wb',filename='',mtime=0) as zipped:shutil.copyfileobj(src,zipped,1048576)
   else:shutil.copyfileobj(src,dst,1048576)
  assert digest(source)==before
  records.append(dict(sourcePath=name,storedPath=str(destination.relative_to(out)),sourceBytes=source.stat().st_size,
      sourceSHA256=before,storedBytes=destination.stat().st_size,storedSHA256=digest(destination),gzipEncoded=encoded))
 manifest=dict(schemaVersion=1,records=records,nativeCases=48,referenceCases=48,cohorts=16,numericalChecks=checks['status'],
   sourceSHA256=read(root/'native-independent-verification.json')['sourceSHA256'],archiveScriptSHA256=digest(Path(__file__)),
   sourceIngestionArchive='2026-09-10-native',predictionBenchmarkIncluded=False,
   scope='All frozen paired DE families and full diagnostics on deposited experimental counts, with independent conditional numerical checks. Concordance does not establish FDR, power, interval coverage or biological accuracy. Original source data and executable binaries remain external.')
 (out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
 print(json.dumps(dict(members=len(records),storedBytes=sum(r['storedBytes'] for r in records),manifestSHA256=digest(out/'manifest.json'))))

if __name__=='__main__':main()
