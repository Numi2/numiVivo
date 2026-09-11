"""Aggregate verified donor scores without turning partial runs into passes."""
from pathlib import Path
import collections,hashlib,json,math,sys,time
study,root=map(Path,sys.argv[1:3]);destination=Path(sys.argv[3]) if len(sys.argv)>3 else root/'summary.json'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=read(study/'protocol.json');folds=protocol['folds'];assert len({f['tag'] for f in folds})==len(folds);plan=read(root/'scoring-plan.json');assert plan['folds']==len(folds) and '5 percent pooled RMSE' in plan['threshold'];plan_hash=sha(root/'scoring-plan.json');completed=[];pending=[];bindings={}
for fold in folds:
 tag=fold['tag'];out=root/tag;checked=out/'scores-verification.json'
 if not checked.exists():pending.append(tag);continue
 check=read(checked);score=read(out/'scores.json');pred=read(out/'verification.json')
 assert check['status']=='passed-independent-sparse-scoring' and check['tag']==tag
 assert sha(out/'scores.json')==check['scoreReportSHA256']
 assert score['tag']==tag and score['queryDonorID']==fold['excludedDonorID'] and score['sourceGenes']==fold['featureCount']
 assert score['scoringPlanSHA256']==plan_hash and sha(out/'verification.json')==score['verifiedPredictionManifestSHA256']
 assert pred['status']=='passed-control-only-predictions' and not pred['errors'] and pred['genes']==fold['featureCount']
 assert sha(out/'scores-by-gene.json.gz')==score['scoresSHA256']
 for shard,b in pred['shardBindings'].items():
  assert sha(out/(shard+'-input.jsonl.gz'))==b['inputSHA256'] and sha(out/(shard+'-output.jsonl.gz'))==b['outputSHA256']
 assert score['scoredGenes']==pred['states']['conditionalPrediction'] and score['scoredGenes']+sum(score['excludedStates'].values())==fold['featureCount']
 for v in score['metrics'].values():
  assert v['genes']==score['scoredGenes'] and v['genes']>0 and v['sumSquaredError']>=0 and math.isfinite(v['sumSquaredError'])
  assert math.isclose(v['RMSE'],math.sqrt(v['sumSquaredError']/v['genes']),rel_tol=1e-12)
 completed.append(dict(score,origin=fold['origin']));bindings[tag]={'scoreReportSHA256':sha(out/'scores.json'),'scoringVerificationSHA256':sha(checked),'predictionVerificationSHA256':sha(out/'verification.json')}
def aggregate(rows,declared):
 complete=len(rows)==len(declared);n=sum(r['scoredGenes'] for r in rows);metrics={};gains={};worse={};excluded=collections.Counter()
 for row in rows:excluded.update(row['excludedStates'])
 if n:
  for name in ['joint','noChange','trainingMeanResponse']:
   ss=math.fsum(r['metrics'][name]['sumSquaredError'] for r in rows);metrics[name]={'genes':n,'sumSquaredError':ss,'RMSE':math.sqrt(ss/n),'MAE':math.fsum(r['metrics'][name]['MAE']*r['scoredGenes'] for r in rows)/n}
  for name in ['noChange','trainingMeanResponse']:
   gains[name]=1-metrics['joint']['RMSE']/metrics[name]['RMSE'] if metrics[name]['RMSE']>0 else None
   worse[name]=[r['tag'] for r in rows if r['metrics']['joint']['RMSE']>r['metrics'][name]['RMSE']]
 passed=all(g is not None and g>=0.05 for g in gains.values()) if complete and n else None
 return {'declaredFolds':len(declared),'scoredFolds':len(rows),'allDeclaredFoldsScored':complete,'declaredSourceGeneFolds':sum(f['featureCount'] for f in declared),'completedSourceGeneFolds':sum(r['sourceGenes'] for r in rows),'scoredGeneFolds':n,'excludedStates':dict(excluded),'metrics':metrics,'relativeRMSEGain':gains,'worseDonors':worse,'pointMassPredictions':sum(r['pointMassPredictions'] for r in rows),'frozenFivePercentDevelopmentCriterionPassed':passed}
report={'status':'complete-development-evaluation' if not pending else 'partial-completed-folds','recordedUnix':time.time(),'scoringPlanSHA256':plan_hash,'pendingFolds':pending,'all':aggregate(completed,folds),'byOrigin':{origin:aggregate([r for r in completed if r['origin']==origin],[f for f in folds if f['origin']==origin]) for origin in sorted({f['origin'] for f in folds})},'bindings':bindings,'qualification':'Pooled matched-gene development RMSE; incomplete groups have no pass/fail result. Missing and excluded predictions remain explicit. Does not establish calibrated uncertainty, independent biological validation, or general phenotype prediction.'}
temp=destination.with_name(destination.name+'.tmp');temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(destination);print({k:report[k] for k in ['status','pendingFolds','all']})
