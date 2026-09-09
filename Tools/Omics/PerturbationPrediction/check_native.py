#!/usr/bin/env python3
"""Real donor-held-out product qualification against the frozen external baselines."""
import argparse,copy,hashlib,json,subprocess,shutil
from pathlib import Path
import anndata as ad
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--study',choices=['kang','hagai'],required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--source-bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--negative-controls',action='store_true');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
original=a.source_bundle/'original.h5ad';mapping=json.loads((a.source_bundle/'plan.json').read_text())['mapping'];obj=ad.read_h5ad(original)
with original.open('rb') as f:sourceSHA=hashlib.file_digest(f,'sha256').hexdigest()
assert sourceSHA==json.loads((a.reference/'inputs.json').read_text())['source']['sourceH5ADSHA256']
control,treated=('ctrl','stim') if a.study=='kang' else ('unstimulated','LPS6')
spec={s['id']:s for s in mapping['samples']};cellSamples=obj.obs[mapping['sampleColumn']].astype(str).tolist();cellDonors=np.array([spec[s]['donorID'] for s in cellSamples]);conditions=np.array([spec[s]['condition'] for s in cellSamples]);donors=sorted(set(cellDonors));commands=[];results=[]
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def run(command,log,error=None):
 r=subprocess.run([str(a.binary.resolve()),*map(str,command)],capture_output=True,text=True);log.write_text(r.stdout+r.stderr)
 commands.append({'command':list(map(str,command)),'returncode':r.returncode,'expectedFailure':error});write(a.out/'commands.json',commands)
 if error:assert r.returncode!=0 and error in r.stdout+r.stderr,(command,r.stdout,r.stderr)
 else:assert r.returncode==0,(command,r.stdout,r.stderr)
for index,donor in enumerate(donors):
 root=a.out/donor;root.mkdir()
 for role,mask in [('training',cellDonors!=donor),('control',(cellDonors==donor)&(conditions==control))]:
  rows=np.flatnonzero(mask).tolist();write(root/(role+'-projection.json'),{'schemaVersion':1,'source':{'bytes':list(bytes.fromhex(sourceSHA))},'provenance':'Real '+a.study+' held-out '+donor+' '+role+'; all original source features retained','observationIndices':rows})
  run(['singlecell-h5ad-project',original,'--plan',root/(role+'-projection.json'),'--output',root/role],root/(role+'-projection.log'))
 trainMapping=copy.deepcopy(mapping);trainMapping['samples']=[s for s in mapping['samples'] if s['donorID']!=donor]
 queryMapping=copy.deepcopy(mapping);queryMapping['samples']=[s for s in mapping['samples'] if s['donorID']==donor and s['condition']==control]
 namespace='HUMAN_GENE_SYMBOL' if a.study=='kang' else 'MOUSE_ENSEMBL_GENE_ID'
 fit={'schemaVersion':1,'mapping':trainMapping,'featureNamespace':namespace,'perturbationID':a.study+'-'+treated,'controlCondition':control,'treatmentCondition':treated,'provenance':'Source-qualified '+a.study+' paired donor study. Training donors exclude '+donor+'. Fixed predeclared four-baseline protocol.'}
 query={'schemaVersion':1,'mapping':queryMapping,'featureNamespace':namespace,'perturbationID':fit['perturbationID']}
 write(root/'fit.json',fit);write(root/'query.json',query)
 run(['singlecell-perturbation-fit',root/'training/projected.h5ad','--plan',root/'fit.json','--output',root/'model'],root/'fit.log')
 run(['singlecell-perturbation-verify',root/'model'],root/'verify-fit.log')
 args=['singlecell-perturbation-predict',root/'control/projected.h5ad','--plan',root/'query.json','--reference',root/'model','--output']
 run(args+[root/'prediction'],root/'predict.log');run(['singlecell-perturbation-prediction-verify',root/'prediction'],root/'verify-prediction.log');run(args+[root/'repeat'],root/'repeat.log')
 assert (root/'prediction/report.json').read_bytes()==(root/'repeat/report.json').read_bytes()
 model=json.loads((root/'model/model.json').read_text());report=json.loads((root/'prediction/report.json').read_text());prediction=report['predictions'][0]
 assert len(report['predictions'])==1 and prediction['group']['donorID']==donor
 with np.load(a.reference/donor/'model.npz',allow_pickle=False) as ref:
  assert model['trainingDonors']==ref['trainingDonors'].tolist();assert model['featureIDs']==ref['featureIDs'].tolist();assert model['selectedFeatureIndices']==ref['selected'].tolist()
  for native,external in [('contextCenters','center'),('contextScales','scale'),('contexts','context'),('meanResponse','mean'),('medianResponse','median'),('dualCoefficients','dual')]:
   np.testing.assert_allclose(model[native],ref[external],rtol=1e-8,atol=1e-10,err_msg=native)
 errors={}
 for estimate in prediction['estimates']:
  name=estimate['baseline']
  with np.load(a.reference/donor/(name+'.npz'),allow_pickle=False) as ref:
   np.testing.assert_allclose(prediction['control'],ref['control'],rtol=1e-10,atol=1e-10)
   for key in ['unclippedResponse','predictedTreated','predictedResponse']:
    np.testing.assert_allclose(estimate[key],ref[key],rtol=1e-8,atol=1e-9,err_msg=name+'/'+key)
   implied=float(np.expm1(ref['predictedTreated']).sum());np.testing.assert_allclose(estimate['impliedCPMSum'],implied,rtol=1e-10,atol=1e-6)
   errors[name]=float(np.max(np.abs(np.asarray(estimate['predictedTreated'])-ref['predictedTreated'])))
 results.append({'donor':donor,'features':len(model['featureIDs']),'trainingDonors':model['trainingDonors'],'maximumPredictionErrors':errors,'maximumSolveResidual':model['maximumSolveResidual'],'exactReplay':True});print(json.dumps(results[-1]),flush=True)
 if index!=0 or not a.negative_controls:continue
 for name,change,error in [
  ('wrong-perturbation',lambda q:q.update(perturbationID='unseen-perturbation'),'identity, namespace or unit mismatch'),
  ('wrong-namespace',lambda q:q.update(featureNamespace='other'),'identity, namespace or unit mismatch'),
  ('wrong-unit',lambda q:q['mapping'].update(countUnit='readCount'),'identity, namespace or unit mismatch'),
  ('treated-query',lambda q:[s.update(condition=treated) for s in q['mapping']['samples']],'control-only query samples'),
  ('wrong-organism',lambda q:[s.update(organism='unknown-organism') for s in q['mapping']['samples']],'feature universe, organism, cell group or donor multiplicity mismatch'),
  ('donor-overlap',lambda q:[s.update(donorID=model['trainingDonors'][0]) for s in q['mapping']['samples']],'query donor overlaps training')]:
  bad=copy.deepcopy(query);change(bad);write(root/(name+'.json'),bad)
  run(['singlecell-perturbation-predict',root/'control/projected.h5ad','--plan',root/(name+'.json'),'--reference',root/'model','--output',root/name],root/(name+'.log'),error);assert not (root/name).exists()
 run(args+[root/'prediction'],root/'overwrite.log','output already exists')
 q=ad.read_h5ad(root/'control/projected.h5ad');q=q[:,::-1].copy();q.write_h5ad(root/'reordered.h5ad')
 run(['singlecell-perturbation-predict',root/'reordered.h5ad','--plan',root/'query.json','--reference',root/'model','--output',root/'reordered'],root/'reordered.log')
 reordered=json.loads((root/'reordered/report.json').read_text())
 for lhs,rhs in zip(reordered['predictions'][0]['estimates'],prediction['estimates']):
  np.testing.assert_array_equal(lhs['predictedTreated'],rhs['predictedTreated'])
 q=q[:,1:].copy();q.write_h5ad(root/'missing-gene.h5ad')
 run(['singlecell-perturbation-predict',root/'missing-gene.h5ad','--plan',root/'query.json','--reference',root/'model','--output',root/'missing'],root/'missing.log','feature universe, organism, cell group or donor multiplicity mismatch')
 shutil.copytree(root/'model',root/'altered-model');path=root/'altered-model/model.json';value=json.loads(path.read_text());value['meanResponse'][0]+=1;write(path,value)
 receiptPath=root/'altered-model/receipt.json';receipt=json.loads(receiptPath.read_text());receipt['result']={'bytes':list(hashlib.sha256(path.read_bytes()).digest())};write(receiptPath,receipt)
 run(['singlecell-perturbation-verify',root/'altered-model'],root/'altered-model.log','perturbation model does not reconstruct')
 assert not list(root.glob('.numivivo-perturbation-*'))
write(a.out/'summary.json',{'study':a.study,'sourceSHA256':sourceSHA,'binarySHA256':hashlib.file_digest(open(a.binary,'rb'),'sha256').hexdigest(),'folds':results,'commands':len(commands)})
