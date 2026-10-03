"""Investigable native spatial-response experiments through the existing lifecycle."""
from pathlib import Path
from functools import lru_cache
import uuid,shutil,subprocess
import numpy as np,anndata as ad
from safetensors.numpy import load_file,save_file
from adapters import WetLabExperimentAdapter
from wetlab import read,write,sha,require,timestamp,inventory
from prepare_spatial_learning import logcounts
from evaluate_spatial_response import metrics

FORMAT='numivivo-learned-spatial-assay/v1'

def load(config):
 p=Path(config).resolve();a=read(p);require(a['format']==FORMAT,'Learned spatial assay format')
 for name,digest in a['files'].items():require(sha(p.parent/name)==digest,'Assay artifact changed: '+name)
 return a,p.parent

@lru_cache(maxsize=4)
def source(path,digest):
 require(sha(path)==digest,'Measured source changed');return ad.read_h5ad(path)

def check(run):
 run=Path(run);s=read(run/'seal.json')
 for f,d in s['files'].items():
  p=run/f;require(p.resolve().is_relative_to(run.resolve()) and not p.is_symlink() and sha(p)==d,'Changed sealed artifact: '+f)
 return read(run/'registration.json')

def check_observation(run):
 run=Path(run);s=read(run/'observation-seal.json');require(s['comparisonSHA256']==sha(run/'comparison.json') and s['predictionSealSHA256']==sha(run/'seal.json'),'Changed observation/evaluation artifact')

def build(inputs,binary,output):
 inputs=Path(inputs).resolve();output=Path(output);output.mkdir(exist_ok=False);prepared=read(inputs/'prepared.json');cohort=Path(prepared['cohort']);receipt=read(cohort/'cohort.json');rows=read(inputs/'chip1-rows.json');query=load_file(str(inputs/'neighborhood/chip1.safetensors'))
 populations=[];common={};meta=[]
 for celltype,hood in sorted(set((r['cellType'],r['neighborhood']) for r in rows)):
  idx=[i for i,r in enumerate(rows) if r['cellType']==celltype and r['neighborhood']==hood and r['role']=='direct'];pid=celltype+'|'+str(hood)
  populations.append({'id':pid,'cellType':celltype,'neighborhood':hood,'anchors':len(idx),'coveredAnchors':sum(rows[i]['referenceCoverage'] for i in idx),'anchorRows':[rows[i]['anchor'] for i in idx],'support':'experimental; four isolated controls, biological qualification blocked'})
  common[pid+'-context']=query['context'][idx].mean(0);common[pid+'-descriptor']=query['descriptor'][idx].mean(0)
 save_file(common,str(output/'common-reference.safetensors'));np.save(output/'model-feature-mask.npy',load_file(str(inputs/'neighborhood/chip2.safetensors'))['mask'][0].astype(bool));shutil.copy2(inputs/'features.json',output/'features.json');shutil.copy2(inputs/'chip1-rows.json',output/'query-rows.json');shutil.copy2(inputs/'preregistration.json',output/'protocol.json');shutil.copy2(inputs/'baselines.safetensors',output/'baselines.safetensors');shutil.copy2(inputs/'baseline-types.json',output/'baseline-types.json')
 models={}
 for variant in ('neighborhood','no-neighborhood','shuffled-neighborhood'):
  folder=inputs/variant;selection=read(folder/'selection.json');weights=folder/'training'/('weights-'+str(selection['selected']['step'])+'.safetensors');destination=output/(variant+'.safetensors');shutil.copy2(weights,destination);models[variant]={'weights':destination.name,'sha256':sha(destination),'selection':selection}
 # A source-bound campaign effect landscape; no outcome RNA is read here.
 prediction=load_file(str(inputs/'neighborhood/test-prediction/prediction.safetensors'));landmarks=[];effects=[]
 for target,celltype,hood,role in sorted(set((r['target'],r['cellType'],r['neighborhood'],r['role']) for r in rows)):
  ix=[i for i,r in enumerate(rows) if (r['target'],r['cellType'],r['neighborhood'],r['role'])==(target,celltype,hood,role)]
  landmarks.append({'target':target,'population':celltype+'|'+str(hood),'role':role,'anchors':len(ix),'evidence':'MODEL INFERENCE','reference':'anchor-matched campaign references; comparison workspace uses one pooled reference'})
  effects.append((prediction['mean'][ix]-query['context'][ix]).mean(0))
 write(output/'landscape.json',landmarks);save_file({'effect':np.asarray(effects,np.float32)},str(output/'landscape.safetensors'))
 shutil.copy2(inputs/'neighborhood/plan.json',output/'plan.json')
 specimens=[]
 for chip in ('chip1','chip2'):
  s=next(x for x in receipt['sources'] if x['chip']==chip);specimens.append({'id':chip,'title':chip+' · '+('held-out mouse 1' if chip=='chip1' else 'training pool, mice 2 + 3'),'cells':s['cells'],'animals':s['animals'],'source':str(cohort/(chip+'.h5ad')),'sourceSHA256':receipt['files'][chip+'.h5ad'],'inferenceSupported':chip=='chip1'})
 artifact={'format':FORMAT,'family':'learned-spatial-response','id':'spatial-response-v03','title':'Which populations respond differently?','runtime':{'binary':str(Path(binary).resolve()),'sha256':sha(binary)},'models':models,'populations':populations,'targets':prepared['targets'],'excludedTargets':[x for x in receipt['targets'] if not x['eligible']],'specimens':specimens,'inputs':str(inputs),'metricPanel':prepared['metricPanelIndices'],'evidence':'MODEL INFERENCE','biologicalPromotion':False,'limits':['Fixed endpoint geometry; no cell trajectories or tissue remodeling','Control-projected annotations are uncertain, not measured cell types','Only four isolated mSafe controls in held-out chip; qualification requires five','One training chip pools two animals; neighborhoods and cells are not independent animals','Native diagonal Gaussian marginals are uncalibrated','Observed comparisons are retrospective after initial campaign reveal; no new independent validation'], 'files':{p.name:sha(p) for p in output.iterdir() if p.is_file()},'sourceCitation':'https://www.nature.com/articles/s41467-026-69677-6'}
 write(output/'assay.json',artifact);return output/'assay.json'

class LearnedSpatialResponseAdapter(WetLabExperimentAdapter):
 adapter_id='learned-spatial-response/v1';family='learned-spatial-response'
 def catalog(self,config):
  a,root=load(config);return {k:a[k] for k in ('id','family','title','populations','targets','excludedTargets','limits','sourceCitation','biologicalPromotion')}|{'adapterID':self.adapter_id,'specimens':[{k:v for k,v in s.items() if k not in ('source','sourceSHA256')} for s in a['specimens']],'features':read(root/'features.json'),'model':'Native MLX cell / reference-neighborhood / perturbation response model','evidence':'MODEL INFERENCE','completedCampaign':read(root/'campaign-summary.json') if (root/'campaign-summary.json').exists() else None}
 def compile(self,config,selection):
  a,root=load(config);require(set(selection)=={'specimen','population','targets'},'Typed specimen/population/targets required');require(selection['specimen']=='chip1','Training specimen is inspectable, not a held-out prediction');require(selection['population'] in {p['id'] for p in a['populations']},'Unsupported reference population');targets=selection['targets'];require(isinstance(targets,list) and 1<=len(targets)<=7 and len(targets)==len(set(targets)) and set(targets)<={t['target'] for t in a['targets']},'Unsupported intervention')
  return {'format':'numivivo-experiment-plan/v2',**selection,'targets':sorted(targets),'time':'measured study endpoint','geometry':'fixed measured chip coordinates','predictionLevel':'population distribution','evidence':'MODEL INFERENCE','support':'experimental; insufficient isolated controls','biologicalPromotion':False}
 def predict(self,config,runtime,workspace,selection):
  a,root=load(config);plan=self.compile(config,selection);require(sha(a['runtime']['binary'])==a['runtime']['sha256'],'Native binary changed');run=Path(workspace)/uuid.uuid4().hex;run.mkdir(parents=True)
  write(run/'registration.json',{'format':FORMAT,'adapterID':self.adapter_id,'family':self.family,'assay':a,'specimen':plan['specimen'],'plan':plan,'createdAt':timestamp(),'ownerSHA256':sha(__file__)})
  for f in a['files']:shutil.copy2(root/f,run/f)
  shutil.copy2(__file__,run/'record-owner.py')
  ref=load_file(str(run/'common-reference.safetensors'));tid={t['target']:i for i,t in enumerate(a['targets'])};arms=[{'target':t,'role':role,'population':plan['population'],'id':t+'|'+role} for t in plan['targets'] for role in ('direct','neighbor')];write(run/'arms.json',arms)
  data={'context':np.asarray([ref[plan['population']+'-context'] for _ in arms],np.float32),'descriptor':np.asarray([ref[plan['population']+'-descriptor'] for _ in arms],np.float32),'target':np.asarray([tid[x['target']] for x in arms],np.int32),'known':np.ones((len(arms),1),np.float32)};data['descriptor'][:,-1]=[float(x['role']=='neighbor') for x in arms]
  for variant in a['models']:
   d={k:v.copy() for k,v in data.items()}
   if variant=='no-neighborhood':d['descriptor'][:,16:32]=0
   # The shuffled model is a training ablation; query context is held fixed
   # for this investigative comparison, explicitly distinct from campaign shuffle.
   save_file(d,str(run/(variant+'-query.safetensors')))
   subprocess.run([a['runtime']['binary'],'spatial-response','predict',str(run/'plan.json'),str(run/(variant+'-query.safetensors')),str(run/variant),str(run/a['models'][variant]['weights'])],check=True,capture_output=True)
  write(run/'seal.json',{'format':FORMAT,'createdAt':timestamp(),'files':inventory(run),'allSelectedArmsSealedTogether':True});return run
 def _observations(self,run):
  reg=check(run);a=reg['assay'];s=next(x for x in a['specimens'] if x['id']=='chip1');return source(s['source'],s['sourceSHA256'])
 def comparison(self,run):
  run=Path(run);reg=check(run);a=reg['assay'];measured=self._observations(run);arms=read(run/'arms.json');rows=read(run/'query-rows.json');query=load_file(str(run/'neighborhood-query.safetensors'));pred={m:load_file(str(run/m/'prediction.safetensors')) for m in a['models']};baseline=load_file(str(run/'baselines.safetensors'));types=[tuple(k) for k in read(run/'baseline-types.json')];mask=measured.var.measured_in_source.to_numpy()&np.load(run/'model-feature-mask.npy');out=[]
  for i,arm in enumerate(arms):
   matched=[r for r in rows if r['target']==arm['target'] and r['role']==arm['role'] and r['cellType']+'|'+str(r['neighborhood'])==arm['population']];ids=sorted(set(j for r in matched for j in r['outcomeRows']));entry={**arm,'outcomeRows':ids,'count':len(ids),'controlReferenceCells':len(set(j for r in matched for j in r['referenceRows'])),'isolatedSafeControls':4,'biologicalUnits':1,'pooledTrainingAnimals':2,'referenceCoveredAnchors':sum(r['referenceCoverage'] for r in matched),'anchorCount':len(matched),'support':'insufficient control coverage','metrics':{},'diagnostics':[]}
   if ids:
    obs=logcounts(measured,ids);control=query['context'][i]
    for name,p in pred.items():entry['metrics'][name]=metrics(p['mean'][i],p['variance'][i],obs,control,a['metricPanel'],mask)
    entry['metrics']['no-change']=metrics(control,None,obs,control,a['metricPanel'],mask)
    key=(arm['target'],arm['role'],arm['population'].split('|')[0]);typeprior=baseline['typePrior'][types.index(key)] if key in types else baseline['prior'][query['target'][i]];entry['metrics']['cell-type-matched']=metrics(control+typeprior,None,obs,control,a['metricPanel'],mask)
    rmse={k:v['rmse'] for k,v in entry['metrics'].items()}
    for name in ('no-neighborhood','no-change','cell-type-matched'):
     if rmse[name]<rmse['neighborhood']:entry['diagnostics'].append(name+' has lower whole-gene error')
   else:entry['diagnostics'].append('No measured target/role population in this specimen')
   if len(ids)<5:entry['diagnostics'].append('Fewer than five measured outcome cells')
   if entry['referenceCoveredAnchors']<entry['anchorCount']:entry['diagnostics'].append('Reference geometry coverage is incomplete')
   entry['diagnostics'].append('Four isolated mSafe controls; minimum five not met');out.append(entry)
  return {'format':FORMAT,'arms':out,'biologicalPromotion':False,'calibration':'Uncalibrated marginal Gaussian; no animal-level confidence interval','predictionSealSHA256':sha(run/'seal.json')}
 def reveal(self,run,runtime):
  run=Path(run);require(not (run/'comparison.json').exists(),'Already revealed');self.verify(run,runtime);result=self.comparison(run);write(run/'comparison.json',result);write(run/'observation-seal.json',{'comparisonSHA256':sha(run/'comparison.json'),'predictionSealSHA256':sha(run/'seal.json'),'revealedAt':timestamp()});return result
 def verify(self,run,runtime):
  run=Path(run);reg=check(run);a=reg['assay'];require(sha(a['runtime']['binary'])==a['runtime']['sha256'],'Exact native runtime unavailable');require(sha(__file__)==reg['ownerSHA256'],'Exact record owner changed');checks={}
  for variant in a['models']:
   out=run/('replay-'+variant+'-'+uuid.uuid4().hex);subprocess.run([a['runtime']['binary'],'spatial-response','predict',str(run/'plan.json'),str(run/(variant+'-query.safetensors')),str(out),str(run/a['models'][variant]['weights'])],check=True,capture_output=True);original=load_file(str(run/variant/'prediction.safetensors'));replay=load_file(str(out/'prediction.safetensors'));require(all(np.array_equal(v,replay[k]) for k,v in original.items()),'Native replay differs');checks[variant]='bit-exact mean and variance tensors'
  if (run/'comparison.json').exists():require(read(run/'comparison.json')==self.comparison(run),'Evaluation replay differs')
  return {'status':'verified','run':run.name,'nativeReplay':checks,'biologicalPromotion':False}
 def summary(self,run):
  run=Path(run);reg=check(run);result={'id':run.name,'adapterID':self.adapter_id,'family':self.family,'registration':reg,'recordDirectory':str(run.resolve()),'arms':read(run/'arms.json'),'revealed':(run/'comparison.json').exists(),'observationsAvailable':True}
  if result['revealed']:
   check_observation(run);result['comparison']=read(run/'comparison.json')
  return result
 def geometry(self,config,specimen):
  a,root=load(config);s=next(x for x in a['specimens'] if x['id']==specimen);data=source(s['source'],s['sourceSHA256']);types=sorted(data.obs.projected_cell_type.astype(str).unique());coords=data.obsm['spatial'];rows=read(root/'query-rows.json') if specimen=='chip1' else [];mapping={r['anchor']:r['cellType']+'|'+str(r['neighborhood']) for r in rows};cells=[[str(data.obs.source_cell_id.iloc[i]),float(x),float(y),types.index(str(data.obs.projected_cell_type.iloc[i])),bool(data.obs.reference_coverage.iloc[i]),mapping.get(i),str(data.obs.guide_assignment.iloc[i]),bool(data.obs.reference_admitted.iloc[i])] for i,(x,y) in enumerate(coords)]
  return {'specimen':specimen,'cells':cells,'types':types,'units':'original GEF chip coordinates','evidence':'MEASURED geometry; control-projected types MODEL INFERENCE','sourceSHA256':s['sourceSHA256'],'measuredGeometryTime':'endpoint, not pretreatment','populations':a['populations'] if specimen=='chip1' else []}
 def feature(self,config,specimen,gene,run=None):
  a,root=load(config);s=next(x for x in a['specimens'] if x['id']==specimen);data=source(s['source'],s['sourceSHA256']);features=read(root/'features.json');require(gene in features,'Gene not in measured union');j=features.index(gene);available=bool(data.var.measured_in_source.iloc[j]);modeled=bool(np.load(root/'model-feature-mask.npy')[j]) and available;ref=np.flatnonzero(data.obs.reference_admitted.to_numpy());raw=data.X[ref,j].toarray().ravel();totals=np.asarray(data.X[ref].sum(1)).ravel();values=np.log1p(raw/np.maximum(1,totals)*1e6)
  result={'gene':gene,'units':'log1p(CPM)','rawUnits':'UMI counts','measuredInSource':bool(data.var.measured_in_source.iloc[j]),'referenceRows':ref.tolist(),'referenceRaw':raw.tolist() if available else [None]*len(raw),'referenceNormalized':values.tolist() if available else [None]*len(raw),'modeled':modeled,'arms':[]}
  landscape=read(root/'landscape.json');effects=load_file(str(root/'landscape.safetensors'))['effect']
  result['populationEffects']=[{**row,'effect':float(effects[i,j]) if modeled else None} for i,row in enumerate(landscape)]
  if run is not None:
   run=Path(run);reg=check(run);require(reg['assay']['id']==a['id'] and reg['specimen']==specimen,'Record/specimen mismatch');arms=read(run/'arms.json');q=load_file(str(run/'neighborhood-query.safetensors'));pred=load_file(str(run/'neighborhood/prediction.safetensors'));rows=read(run/'query-rows.json');revealed=(run/'comparison.json').exists()
   for i,arm in enumerate(arms):
    y={'id':arm['id'],'target':arm['target'],'role':arm['role'],'control':float(q['context'][i,j]),'predicted':float(pred['mean'][i,j]),'variance':float(pred['variance'][i,j]),'evidence':'MODEL INFERENCE','uncertainty':'uncalibrated distribution; not a cell-specific interval','observed':None,'residual':None,'observedValues':[],'observedRaw':[],'observedRows':[]}
    if revealed:
     check_observation(run)
     ids=sorted(set(k for r in rows if r['target']==arm['target'] and r['role']==arm['role'] and r['cellType']+'|'+str(r['neighborhood'])==reg['plan']['population'] for k in r['outcomeRows']))
     if ids:
      counts=data.X[ids,j].toarray().ravel();tot=np.asarray(data.X[ids].sum(1)).ravel();v=np.log1p(counts/np.maximum(1,tot)*1e6);y.update(observed=float(v.mean()),residual=float(pred['mean'][i,j]-v.mean()),observedValues=v.tolist(),observedRaw=counts.tolist(),observedRows=ids)
    if not modeled:
     y.update(predicted=None,variance=None,residual=None,evidence='UNAVAILABLE',uncertainty='gene unmeasured in query or training source')
    if not available:y.update(control=None,observed=None,observedValues=[],observedRaw=[],observedRows=[])
    result['arms'].append(y)
  return result
