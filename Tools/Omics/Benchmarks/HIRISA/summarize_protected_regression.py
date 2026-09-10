#!/usr/bin/env python3
"""Compare the single frozen development candidate with the published native fit."""
import argparse,statistics
from pathlib import Path
from prepare_annotation_retention import read,write
from check_integration_response import sha


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ('study','root'):p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();s=a.study;r=a.root;e=r/'evaluation';complete=read(r/'complete.json')
 assert complete['allNumericalChecksPassed'] and complete['allOriginalCellsPreserved']
 for name,h in complete['results'].items():assert sha(e/name)==h
 paths={'nativeAnnotation':s/'annotation-retention/native.json','protectedAnnotation':e/'protected.json','nativeResponse':s/'integration-response-v2/native-retry.json','baselineResponse':s/'integration-response-v2/baseline.json','protectedResponse':e/'response.json','nativeWithin':s/'program-calibration/results/native.json','protectedWithin':e/'within-programs.json','nativeMixed':s/'integration-programs/native.json','protectedMixed':e/'mixed-programs.json'}
 x={k:read(v) for k,v in paths.items()};old=x['nativeAnnotation'];new=x['protectedAnnotation'];legacy=read(e/'legacy-reconstructed.json')
 assert len(old['folds'])==len(new['folds'])==len(legacy['folds'])==114
 assert legacy['comparisons']==old['comparisons']
 for f,g in zip(old['folds'],legacy['folds']):
  for k in ('id','confusion','metrics','accuracy','balancedRecall','trainingCounts','queryCounts'):assert f[k]==g[k]
 aa=old['comparisons'];bb=new['comparisons'];assert len(aa)==len(bb)==471
 eligible=[]
 for i,(f,g) in enumerate(zip(aa,bb)):
  for k in ('stratum','label','labelName','rare','sufficientSupport','controlSensitive','baselineMeanRecall'):assert f[k]==g[k]
  if f['sufficientSupport'] and f['controlSensitive']:eligible.append(i)
 delta=[bb[i]['candidateMeanRecall']-aa[i]['candidateMeanRecall'] for i in eligible]
 strata=read(r/'execution-freeze.json')['strata']
 def desc(v):return dict(**strata[v['stratum']],label=v['labelName'],rare=v['rare'])
 annotation={name:dict(marginFailures=sum(v['marginFailure'] for v in c),rareMarginFailures=sum(v['marginFailure'] and v['rare'] for v in c),eligibleMeanRecall=statistics.mean(c[i]['candidateMeanRecall'] for i in eligible)) for name,c in [('native',aa),('protected',bb)]}
 annotation.update(eligible=len(eligible),eligibleRare=sum(aa[i]['rare'] for i in eligible),insufficient=471-len(eligible),resolved=[desc(f) for f,g in zip(aa,bb) if f['marginFailure'] and not g['marginFailure']],introduced=[desc(g) for f,g in zip(aa,bb) if not f['marginFailure'] and g['marginFailure']],meanRecallChange=statistics.mean(delta),minimumRecallChange=min(delta),maximumRecallChange=max(delta),improved=sum(v>0 for v in delta),worsened=sum(v<0 for v in delta),unchanged=sum(v==0 for v in delta))
 response={}
 for name in ('baseline','native','protected'):
  v=x[name+'Response'];scatter=[z['fraction'] for z in v['conditionalDonorScatter'] if z['status']=='measured']
  response[name]=dict(conditionalGroups=len(scatter),unweightedMeanConditionalDonorScatterFraction=statistics.mean(scatter))
  if name!='baseline':response[name].update(contrasts=len(v['comparisons']),allResponseGatesPassed=v['allResponseDiagnosticGatesPassed'],unweightedMeanRelativeResponseDrift=statistics.mean(z['meanRelativeResponseDrift'] for z in v['comparisons']))
 programs={}
 for name in ('nativeWithin','protectedWithin','nativeMixed','protectedMixed'):
  v=x[name];sens=[z for z in v['comparisons'] if z['controlSensitive']]
  programs[name]=dict(comparisons=len(v['comparisons']),sensitive=len(sens),insufficient=len(v['comparisons'])-len(sens),sensitiveMarginFailures=sum(not z['gates']['meanPreserved'] or not z['gates']['everyFoldPreserved'] for z in sens),completeQualified=v['allProgramGradientGatesPassed'])
 write(r/'summary.json',dict(status='completed-development-hypothesis-not-promoted',annotation=annotation,response=response,programs=programs,inputSHA256={k:sha(v) for k,v in paths.items()},completeSHA256=sha(r/'complete.json'),summarizerSHA256=sha(Path(__file__)),scope='One post-result fixed-membership regression experiment; no native iterative implementation or independent biological validation. Descriptive means are unweighted across declared comparisons, not cells.'))
 print(annotation);print(response);print(programs)
if __name__=='__main__':main()
