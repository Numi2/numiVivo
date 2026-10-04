"""Thin, portable laboratory adapter over frozen native cellular operations.

No fitting, neural implementation or alternative scientific scoring lives here.
The initial artifact exposes historical population means, not single-cell data.
"""
import shutil,tempfile,uuid
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,inventory,timestamp
from train_intervention_design import native
from design_objective import objective,utility,selection_score

SCOPE='Previously exposed technical development data; independent biological validation and population-mean uncertainty remain unqualified'

def load(config,verify=False):
 config=Path(config);a=read(config);root=config.parent
 if verify:
  for p,h in a['files'].items():require(sha(root/p)==h,'Frozen cellular artifact changed: '+p)
  for p,h in a['scientificOwners'].items():require(sha(Path(__file__).with_name(p))==h,'Scientific owner changed: '+p)
 return a,root

def selected(root,s):
 return [c for c in read(root/'cases.json') if c['specimenID']==s['specimen'] and c['populationID']==s['population'] and c['conditionID']==s['condition'] and (not s.get('targets') or c['target'] in s['targets'])]

def spec(genes):return {'genes':genes,'preserveGenes':[],'penalty':0}

def check(run):
 run=Path(run);reg=read(run/'registration.json')
 for p,h in read(run/'seal.json')['files'].items():require(sha(run/p)==h,'Sealed prediction changed: '+p)
 require(reg['ownerSHA256']==sha(__file__),'Use retained cellular adapter version')
 load(run/'artifact/assay.json',True)
 if (run/'comparison.json').exists():
  for p,h in read(run/'observation-seal.json')['files'].items():require(sha(run/p)==h,'Revealed observation changed: '+p)
 return reg

class CellularResponseAdapter:
 family='learned-cellular-response';adapter_id='numivivo-frozen-cellular/v1'
 def catalog(self,config):
  a,root=load(config);return {k:a[k] for k in ('id','title','family','modelVersion','presentation','geometry','specimens','populations','conditions','targets','features','evidence','limits','biologicalPromotion')}|{'capabilities':self.capabilities(config)}
 def capabilities(self,config):
  a,root=load(config);return {k:a[k] for k in ('presentation','geometry','specimens','populations','conditions','targets','features')}|{'prediction':True,'observationComparison':True,'readouts':['gene','population-summary','objective-coverage'],'measuredFeatures':a['features'],'observationAccess':a['observationAccess'],'preservationSupported':False}
 def validate_selection(self,config,selection):
  a,root=load(config);s=dict(selection)
  require(set(s)<= {'specimen','population','condition','targets','objective','preregistration'},'Unknown cellular selection field')
  require(s.get('condition') in {c['id'] for c in a['conditions']},'Unsupported condition: no registered matched control input')
  require(s.get('specimen') in {c['id'] for c in a['specimens']},'Unknown cellular specimen')
  rows=selected(root,{**s,'targets':[]});require(rows,'Population and condition have no matched source input')
  require(isinstance(s.get('targets'),list) and len(s['targets'])==len(set(s['targets'])),'Distinct intervention list required')
  require(set(s['targets'])<={c['target'] for c in rows},'Unsupported intervention in selected population and condition')
  obj=s.get('objective');require(isinstance(obj,dict),'Preregistered molecular objective required');objective(a['features'],obj)
  require(not obj['preserveGenes'] and obj['penalty']==0,'Receiving populations unavailable; preservation disabled')
  return s
 def support(self,config,s,genes):
  a,root=load(config);rows=selected(root,{**s,'targets':[]});controls={x for c in rows for x in c['controlCellIDs']}
  return {'population':s['population'],'referenceCells':len(controls),'selectedPopulationReferenceCells':len(controls),'reference':{'referenceCells':len(controls),'uniqueCells':len(controls),'safeHarbourCells':None,'barcodeNegativeCells':None},'candidates':[{'target':c['target'],'role':'direct','controlCells':len(set(c['controlCellIDs'])),'outcomeCells':len(set(c['outcomeCellIDs'])),'canExecute':True,'recommendationSupported':False,'biologicalUnits':None,'reason':SCOPE} for c in rows],'genes':[{'gene':g,'measured':g in a['features'],'modeled':g in a['features']} for g in genes],'recommendationSupported':False,'reliableWinner':False,'scope':SCOPE,'reason':SCOPE,'preservationSupported':False}
 def inspect_objective(self,config,s,obj):
  a,root=load(config);genes=obj['genes'];rows=selected(root,{**s,'targets':[]});control={x for c in rows for x in c['controlCellIDs']};eligible=bool(genes) and len(genes)==len(set(genes)) and set(genes)<=set(a['features']) and not obj['preserveGenes'] and obj['penalty']==0
  effects=[]
  if eligible:
   q=load_file(str(root/'query.safetensors'));p=load_file(str(root/'preview.safetensors'));ix=objective(a['features'],obj)
   for c in rows:effects.append({'target':c['target'],'role':'direct','population':s['population'],'predictedUtility':float(utility(p['mean'][c['index']]-q['context'][c['index']],None,ix)[0]),'evidence':'MODEL INFERENCE'})
  return {'genes':[{'gene':g,'eligible':g in a['features'],'referenceQuantity':'MEASURED' if g in a['features'] else 'UNAVAILABLE','predictionQuantity':'MODEL INFERENCE' if g in a['features'] else 'UNAVAILABLE','selectedPopulationReferenceCells':len(control),'referenceCells':len(control),'measurementUnits':'log1p(CPM)','referenceScope':'selected source/context unique control cells; raw counts not exposed by this frozen artifact'} for g in genes],'canExecute':eligible and bool(rows),'candidateEffects':effects,'previewReference':'PREVIEW from retained frozen checkpoint, not this draft’s sealed prediction','objective':'reduce mean log1p(CPM) across exactly these genes','preservationSupported':False,'preservationReason':'No receiving-population model in cellular artifact','biologicalPromotion':False,'qualification':SCOPE}
 def predict(self,config,runtime,workspace,selection):
  a,root=load(config,True);s=self.validate_selection(config,selection);require(s['targets'],'Select explicit supported interventions');cases=selected(root,s);require(cases,'Select supported interventions')
  run=Path(workspace)/uuid.uuid4().hex;run.mkdir();shutil.copytree(root,run/'artifact');ix=[c['index'] for c in cases];query=load_file(str(root/'query.safetensors'));save_file({k:v[ix] for k,v in query.items()},str(run/'query.safetensors'))
  reg={'format':'numivivo-cellular-experiment/v1','family':self.family,'adapterID':self.adapter_id,'createdAt':timestamp(),'assay':a,'selection':s,'specimen':s['specimen'],'plan':s,'cases':cases,'ownerSHA256':sha(__file__),'featureAxisSHA256':sha(root/'features.json'),'querySHA256':sha(run/'query.safetensors'),'objective':s['objective'],'scoringVersion':a['scientificOwners'],'biologicalPromotion':False,'scope':SCOPE}
  # Registration is durable before the native owner can execute.
  write(run/'registration.json',reg)
  native(run/'artifact'/a['runtime']['binary'],'predict',run/'artifact/plan.json',run/'query.safetensors',run/'prediction',run/'artifact/weights.safetensors')
  write(run/'seal.json',{'createdAt':timestamp(),'files':inventory(run)})
  return run
 def readout(self,config,record,query,selection=None):
  if record:
   run=Path(record);reg=check(run);root=run/'artifact';a=reg['assay'];cases=reg['cases'];q=load_file(str(run/'query.safetensors'));p=load_file(str(run/'prediction/prediction.safetensors'));revealed=(run/'comparison.json').exists();y=load_file(str(run/'observations.safetensors'))['mean'] if revealed else None
  else:
   a,root=load(config);cases=selected(root,{**selection,'targets':[]});q0=load_file(str(root/'query.safetensors'));p0=load_file(str(root/'preview.safetensors'));indices=[c['index'] for c in cases];q={k:v[indices] for k,v in q0.items()};p={k:v[indices] for k,v in p0.items()};y=None;revealed=False
  require(query.get('mode')!='measured' or revealed,'Measured endpoints remain closed for this experiment; authorize reveal explicitly')
  require(query.get('kind','gene') in ('gene','population-summary','objective-coverage'),'Unsupported readout');gene=query.get('gene',a['features'][0]);require(gene in a['features'],'Feature not measured on this artifact axis');j=a['features'].index(gene);limit=query.get('limit',128);offset=query.get('offset',0);require(type(limit)==int and 1<=limit<=256 and type(offset)==int and offset>=0,'Bounded pagination required');arms=[]
  for i,c in enumerate(cases):
   if query.get('target') and query['target']!=c['target']:continue
   observed=float(y[i,j]) if y is not None else None
   arms.append({'target':c['target'],'role':'direct','population':c['populationID'],'control':float(q['context'][i,j]),'predicted':float(p['mean'][i,j]),'variance':float(p['variance'][i,j]),'observed':observed,'residual':float(p['mean'][i,j])-observed if observed is not None else None,'controlCells':len(set(c['controlCellIDs'])),'outcomeCells':len(set(c['outcomeCellIDs'])),'observedValues':[],'evidence':'MODEL INFERENCE','biologicalUnits':None})
  baseline=load_file(str(root/'baseline.safetensors'))['delta']
  for arm in arms:
   c=next(c for c in cases if c['target']==arm['target']);arm['matchedBaseline']=arm['control']+float(baseline[c['index'],j]);arm['noChangeResidual']=arm['control']-arm['observed'] if arm['observed'] is not None else None;arm['matchedResidual']=arm['matchedBaseline']-arm['observed'] if arm['observed'] is not None else None
  return {'gene':gene,'units':'log1p(CPM)','presentation':'population','modeled':True,'arms':arms[offset:offset+limit],'total':len(arms),'revealed':revealed,'access':'REVEALED OBSERVATION' if revealed else 'SEALED PREDICTION' if record else 'PREVIEW','scope':SCOPE,'distributionStatus':'Population mean only; individual-cell distributions unavailable in this frozen artifact'}
 def evaluate_objective(self,config,record,obj):
  run=Path(record);reg=check(run);require(obj==reg['objective'],'Objective differs from preregistration');a=reg['assay'];ix=objective(a['features'],obj);q=load_file(str(run/'query.safetensors'));p=load_file(str(run/'prediction/prediction.safetensors'));y=load_file(str(run/'observations.safetensors'))['mean'] if (run/'comparison.json').exists() else None;b=load_file(str(run/'artifact/baseline.safetensors'))['delta'];arms=[];pred=[];obs=[];base=[]
  for i,c in enumerate(reg['cases']):
   pu=float(utility(p['mean'][i]-q['context'][i],None,ix)[0]);ou=float(utility(y[i]-q['context'][i],None,ix)[0]) if y is not None else None;bu=float(utility(b[c['index']],None,ix)[0]);pred.append(pu);obs.append(ou);base.append(bu)
   arms.append({'target':c['target'],'role':'direct','population':c['populationID'],'predictedUtility':pu,'observedUtility':ou,'trainingMeanUtility':bu,'baselineScope':c['baselineScope'],'controlCells':len(set(c['controlCellIDs'])),'outcomeCells':len(set(c['outcomeCellIDs'])),'genes':[{'gene':g,'control':float(q['context'][i,j]),'predicted':float(p['mean'][i,j]),'observed':float(y[i,j]) if y is not None else None,'residual':float(p['mean'][i,j]-y[i,j]) if y is not None else None,'trainingMeanDelta':float(b[c['index'],j])} for g,j in zip(obj['genes'],ix[0])]})
  candidates=[c['target'] for c in reg['cases']]+['no-intervention'];scores=selection_score(candidates,pred+[0.],obs+[0.],base+[0.]) if y is not None else None
  choice=sorted(zip(candidates,pred+[0.]),key=lambda x:(-x[1],x[0]))[0][0]
  return {'objectiveGenes':obj['genes'],'arms':arms,'numericalChoice':choice,'scores':scores,'reliableWinner':False,'scope':SCOPE,'preservation':'UNAVAILABLE'}
 def reveal(self,run,runtime):
  run=Path(run);reg=check(run);require(not (run/'comparison.json').exists(),'Already revealed');source=load_file(str(run/'artifact/observations.safetensors'));ix=[c['index'] for c in reg['cases']];save_file({k:v[ix] for k,v in source.items()},str(run/'observations.safetensors'))
  comparison=self._measure(run);write(run/'comparison.json',comparison);write(run/'observation-seal.json',{'files':{p:sha(run/p) for p in ('observations.safetensors','comparison.json')}});return comparison
 def _measure(self,run):
  from train_cohort_recovery import measure
  reg=read(run/'registration.json');q=load_file(str(run/'query.safetensors'));y=load_file(str(run/'observations.safetensors'));p=load_file(str(run/'prediction/prediction.safetensors'));rows=[{'source':c['specimenID'],'context':c['conditionID'],'target':c['target']} for c in reg['cases']]
  return {'metrics':measure(p['mean'],{'context':q['context'],'observed':y['mean'],'mask':y['mask']},rows),'scope':SCOPE,'biologicalPromotion':False}
 def verify(self,run,runtime):
  run=Path(run);reg=check(run)
  with tempfile.TemporaryDirectory(prefix='cellular-replay-') as d:
   out=Path(d)/'prediction';native(run/'artifact'/reg['assay']['runtime']['binary'],'predict',run/'artifact/plan.json',run/'query.safetensors',out,run/'artifact/weights.safetensors');old=load_file(str(run/'prediction/prediction.safetensors'));new=load_file(str(out/'prediction.safetensors'));require(old.keys()==new.keys() and all(np.array_equal(old[k],new[k]) for k in old),'Native replay differs')
  if (run/'comparison.json').exists():require(self._measure(run)==read(run/'comparison.json'),'Scientific comparison differs')
  return {'status':'passed','exactPredictionReplay':True,'observationComparisonVerified':(run/'comparison.json').exists(),'biologicalPromotion':False}
 def summary(self,run):
  run=Path(run);reg=check(run);return {'id':run.name,'registration':{k:v for k,v in reg.items() if k not in ('cases',)},'revealed':(run/'comparison.json').exists(),'comparison':read(run/'comparison.json') if (run/'comparison.json').exists() else None,'evidence':'MODEL INFERENCE','scope':SCOPE}
