"""Retrospective complete-target diagnostic; no model fitting or new validation."""
from pathlib import Path
import subprocess,gzip,json,hashlib,statistics,argparse

def main():
 a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);args=a.parse_args()
 root=Path(__file__).resolve().parents[5]
 base='Tools/Omics/PerturbationPrediction/Replogle2020/evidence/2026-09-11-target-kernel/'
 def blob(name):return subprocess.check_output(['git','-C',str(root),'show','HEAD:'+base+name])
 manifest=json.loads(blob('manifest.json'));record=next(x for x in manifest['records'] if x['sourcePath']=='scoring.json')
 stored=blob(record['storedPath']);assert len(stored)==record['storedBytes'] and hashlib.sha256(stored).hexdigest()==record['storedSHA256']
 raw=gzip.decompress(stored);assert len(raw)==record['sourceBytes'] and hashlib.sha256(raw).hexdigest()==record['sourceSHA256']
 files={}
 for f in json.loads(raw)['sourceFiles']:
  b=f['rawUTF8'].encode();assert len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['SHA256']
  if f['path'].endswith('.json'):files[f['path']]=json.loads(b)
 path='study/target-kernel-check/'
 scores=files[path+'scores.json'];summary=files[path+'summary.json']
 rows=[x for x in scores if x['panel']=='allGenes'];groups=sorted({x['gemgroup'] for x in rows});targets=sorted({x['target'] for x in rows})
 assert len(groups)==5 and len(targets)==30 and len(rows)==750
 index={(x['gemgroup'],x['target'],x['method']):x for x in rows};assert len(index)==750
 methods=['noChange','meanSingleResponse','meanSupportedResponse','goRidgeFixed','shuffledGoRidgeFixed']
 result=[]
 for target in targets:
  folds=[]
  for group in groups:
   values={m:index[group,target,m]['rmse'] for m in methods}
   assert all(index[group,target,m]['supported'] for m in methods)
   assert values['meanSingleResponse']==values['meanSupportedResponse']
   folds.append(dict(gemgroup=group,rmse=values,worseThanMean=values['goRidgeFixed']>values['meanSingleResponse'],worseThanNoChange=values['goRidgeFixed']>values['noChange'],worseThanShuffle=values['goRidgeFixed']>values['shuffledGoRidgeFixed']))
  result.append(dict(target=target,folds=folds,worseThanMeanGroups=sum(f['worseThanMean'] for f in folds),worseThanNoChangeGroups=sum(f['worseThanNoChange'] for f in folds),worseThanShuffleGroups=sum(f['worseThanShuffle'] for f in folds),meanRMSE={m:statistics.mean(f['rmse'][m] for f in folds) for m in methods}))
 # Reproduce original complete failure lists, independently of target aggregation.
 for s in summary:
  if s['panel']!='allGenes':continue
  group=s['gemgroup']
  for baseline,key in [('meanSingleResponse','kernelWorseThanMean'),('noChange','kernelWorseThanNoChange')]:
   assert set(s[key])=={t for t in targets if index[group,t,'goRidgeFixed']['rmse']>index[group,t,baseline]['rmse']}
 report=dict(schemaVersion=1,scope='Retrospective all-target consistency across technical groups; no new fit or independent validation',sourceScoringSHA256=record['sourceSHA256'],groups=groups,targets=30,folds=150,results=result,histograms={field:{str(n):sum(t[field]==n for t in result) for n in range(6)} for field in ['worseThanMeanGroups','worseThanNoChangeGroups','worseThanShuffleGroups']})
 with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
 print(json.dumps(report['histograms'],indent=2))
 print(json.dumps({field:[t['target'] for t in result if t[field]==5] for field in report['histograms']},indent=2))
if __name__=='__main__':main()
