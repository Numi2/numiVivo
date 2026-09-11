"""Validate pooled results from per-gene rows and reject altered evidence."""
from pathlib import Path
import gzip,json,math,shutil,subprocess,sys
study,root=map(Path,sys.argv[1:3]);summary=json.loads((root/'summary.json').read_text());out=root/'summary-contract';out.mkdir(exist_ok=False);recipe=Path(__file__).with_name('summarize.py');errors={k:[] for k in ['joint','noChange','trainingMeanResponse']}
for tag in summary['bindings']:
 rows=json.loads(gzip.decompress((root/tag/'scores-by-gene.json.gz').read_bytes()))
 for row in rows:
  if row['predictionLog1pCPM'] is None:continue
  y=row['observedLog1pCPM']
  for model,key in [('joint','predictionLog1pCPM'),('noChange','controlLog1pCPM'),('trainingMeanResponse','trainingMeanResponseLog1pCPM')]:errors[model].append(row[key]-y)
for name,values in errors.items():
 metrics=summary['all']['metrics'][name];assert len(values)==metrics['genes'];assert math.isclose(math.sqrt(math.fsum(v*v for v in values)/len(values)),metrics['RMSE'],rel_tol=1e-12);assert math.isclose(math.fsum(abs(v) for v in values)/len(values),metrics['MAE'],rel_tol=1e-12)
assert summary['status']=='partial-completed-folds' and summary['all']['frozenFivePercentDevelopmentCriterionPassed'] is None
cases={};first=next(iter(summary['bindings']))
for name in ['no-completed-folds','changed-score','changed-plan']:
 d=out/name;d.mkdir();shutil.copyfile(root/'scoring-plan.json',d/'scoring-plan.json')
 if name!='no-completed-folds':
  target=d/first;target.mkdir()
  for file in ['scores-verification.json','scores.json','verification.json']:shutil.copyfile(root/first/file,target/file)
  if name=='changed-score':
   p=target/'scores.json';v=json.loads(p.read_text());v['metrics']['joint']['RMSE']+=1;p.write_text(json.dumps(v))
  else:
   p=d/'scoring-plan.json';v=json.loads(p.read_text());v['primaryEndpoint']='altered endpoint';p.write_text(json.dumps(v))
 p=subprocess.run([sys.executable,str(recipe),str(study),str(d)],capture_output=True,text=True);(d/'stdout.log').write_text(p.stdout);(d/'stderr.log').write_text(p.stderr)
 if name=='no-completed-folds':
  assert p.returncode==0;v=json.loads((d/'summary.json').read_text());assert v['all']['scoredFolds']==0 and v['all']['frozenFivePercentDevelopmentCriterionPassed'] is None and v['all']['metrics']=={};cases[name]='partial, no pass'
 else:assert p.returncode!=0 and not (d/'summary.json').exists();cases[name]='rejected before summary publication'
report={'status':'passed-real-summary-and-failure-contracts','geneFoldsCompared':len(errors['joint']),'cases':cases};(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
