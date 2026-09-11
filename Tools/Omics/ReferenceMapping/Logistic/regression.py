#!/usr/bin/env python3
"""Exercise native logistic lifecycle and exact historical default kNN compatibility."""
import argparse,copy,hashlib,json,math,shutil,subprocess
from pathlib import Path
from prepare import sha,write
from bundles import pack

def main():
 p=argparse.ArgumentParser();p.add_argument('--restored',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
 ref=a.restored/'reference';query=a.restored/'inputs/query.h5ad';base=json.loads((a.restored/'mapped/plan.json').read_text())
 def run(label,args,success=True):
  r=subprocess.run([str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr);commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=success));write(a.out/'commands.json',commands);assert (r.returncode==0)==success,(label,r.stderr[-1500:])
 def mapping(label,source,plan,success=True,reference=ref):
  write(a.out/(label+'-plan.json'),plan);dest=a.out/label
  run(label,['singlecell-reference-map',source,'--plan',a.out/(label+'-plan.json'),'--reference',reference,'--output',dest],success);assert dest.exists()==success;return dest
 run('restored-model-verify',['singlecell-reference-verify',ref]);run('restored-prediction-verify',['singlecell-reference-map-verify',a.restored/'mapped'])
 for label,edit in [('label-leak',lambda p:p['mapping'].update(groupColumn='cell_type')),('classifier-budget',lambda p:p.update(maximumClassifierOperations=1)),('projection-budget',lambda p:p.update(maximumProjectionUpdates=1)),('namespace',lambda p:p.update(featureNamespace='unmatched'))]:
  plan=copy.deepcopy(base);edit(plan);mapping(label,query,plan,False)
 mapping('missing-gene',a.fixtures/'missing-gene.h5ad',base,False)
 empty=mapping('empty',a.fixtures/'empty-library.h5ad',base);er=json.loads((empty/'report.json').read_text());assert len(er['cells'])==3 and er['classifierOperations']==2*21*len(er['classes']);zero=[c for c in er['cells'] if c['totalCounts']==0];assert len(zero)==1;row=zero[0];assert row.get('candidateLabel') is None and row.get('classProbabilities') is None and row.get('scores') is None;assert all(c.get('classProbabilities') is not None for c in er['cells'] if c['totalCounts']>0)
 reordered=mapping('reordered',a.fixtures/'reordered-unlabelled.h5ad',base);old=json.loads((a.restored/'mapped/report.json').read_text());new=json.loads((reordered/'report.json').read_text())
 assert [x['candidateLabel'] for x in old['cells']]==[x['candidateLabel'] for x in new['cells']]
 for x,y in zip(old['cells'],new['cells']):
  assert all(math.isclose(u,v,rel_tol=1e-11,abs_tol=1e-11) for u,v in zip(x['classProbabilities'],y['classProbabilities']))
 run('existing-output',['singlecell-reference-map',query,'--plan',a.restored/'mapped/plan.json','--reference',ref,'--output',empty],False)
 overlap=copy.deepcopy(base);overlap['mapping']['samples']=json.loads((ref/'plan.json').read_text())['mapping']['samples'];mapping('overlap',a.restored/'inputs/train.h5ad',overlap,False)
 def clone(source,target):
  subprocess.run(['/bin/cp','-c',str(source),str(target)],check=True);return str(target)
 tamper=a.out/'tampered-reference';shutil.copytree(ref,tamper,copy_function=clone)
 model=json.loads((tamper/'model.json').read_text());model['logistic']['parameters'][0][0]+=0.1;write(tamper/'model.json',model)
 run('tampered-model',['singlecell-reference-verify',tamper],False)
 receipt=json.loads((tamper/'receipt.json').read_text());receipt['result']['bytes']=list(hashlib.sha256((tamper/'model.json').read_bytes()).digest());write(tamper/'receipt.json',receipt)
 run('tampered-rehashed-model',['singlecell-reference-verify',tamper],False)
 # The default classifier still emits exactly the historical model/report.
 fit=json.loads((ref/'plan.json').read_text());fit.pop('logistic');write(a.out/'default-fit.json',fit)
 default=a.out/'default-reference';run('default-fit',['singlecell-reference-fit',a.restored/'inputs/train.h5ad','--plan',a.out/'default-fit.json','--output',default])
 defaultMap=mapping('default-mapped',query,base,reference=default)
 expected=json.loads((a.fixtures/'default-hashes.json').read_text())
 assert sha(default/'model.json')==expected['model.json'] and sha(defaultMap/'report.json')==expected['report.json']
 assert 'logistic' not in json.loads((default/'model.json').read_text()) and 'classifierOperations' not in json.loads((defaultMap/'report.json').read_text())
 assert not list(a.out.glob('.numivivo-reference-*'))
 write(a.out/'checks.json',dict(status='passed',commands=len(commands),restoredNativeModelAndPredictionVerified=True,emptyQueriesUnmapped=True,reversedFeaturesExactLabels=True,reversedFeatureProbabilityTolerance=1e-11,workBudgetsAndLabelLeakRejected=True,sourceOverlapRejected=True,tamperedModelAndRehashedModelRejected=True,defaultHistoricalModelAndReportBytesExact=True,binarySHA256=sha(a.binary),checkerSHA256=sha(__file__)))
 print((a.out/'checks.json').read_text())
if __name__=='__main__':main()
