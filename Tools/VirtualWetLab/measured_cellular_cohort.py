"""Bounded measured-evidence access to existing native count-store artifacts.

No model or second scientific scorer. Access is checked per source row before
reading expression; reserved development responses are never previewed.
"""
from functools import lru_cache
from pathlib import Path
import math
import numpy as np
from wetlab import read,require,sha

BLOCK='No compatible trained model is registered for this cohort. Register an artifact with this source, context and feature axis; measured exploration remains available.'
SCOPE='Measured technical development data. Biological-unit independence unresolved; cell counts do not establish biological replication.'
DTYPE=np.dtype([('row','<u4'),('feature','<u4'),('count','<u8')])

@lru_cache(maxsize=8)
def dataset(path,fingerprint,artifact_stats):
 root=Path(path).parent;a=read(path);require(a['datasetFingerprint']==fingerprint,'Cohort identity changed')
 require(sha(root/'registration.json')==fingerprint,'Frozen cohort registration differs')
 metadata_files={'registration.json','observation-access.json'}
 for source in a['sources']:
  metadata_files.update(source['directory']+'/'+name for name in ('rows.json','source.json','admission.json','row-offsets.npy','store/metadata.json','store/quality.json','store/receipt.json'))
 for name in metadata_files:require(sha(root/name)==a['files'][name],'Admitted metadata changed: '+name)
 sources={}
 for source in a['sources']:
  declared=next(s for s in read(root/'registration.json')['sources'] if s['id']==source['id'])
  require(all(source[k]==declared[k] for k in ('context','study','conditionID','species','interventionType')),'Source context binding differs')
  d=root/source['directory'];metadata=read(d/'store/metadata.json');rows=read(d/'rows.json');quality=read(d/'store/quality.json');records=None;offsets=np.load(d/'row-offsets.npy',mmap_mode='r');lookup={};groups={}
  for j,f in enumerate(metadata['features']):
   lookup.setdefault(f['name'],[]).append(j);lookup.setdefault(f['id'],[]).append(j) if f['id']!=f['name'] else None
  for i,row in enumerate(rows):groups.setdefault(row['target'] or 'no-intervention',[]).append(i)
  sources[source['id']]={'source':source,'metadata':metadata,'rows':rows,'totals':quality['rowTotals'],'records':records,'offsets':offsets,'features':lookup,'groups':groups,'countPath':str(d/'store/counts.bin'),'expressionFiles':tuple((str(root/name),h,(root/name).stat().st_ino,(root/name).stat().st_size,(root/name).stat().st_mtime_ns) for name,h in a['files'].items() if name.startswith(source['directory']+'/') or name=='source-metadata/'+source['id']+'.json')}
 return a,sources

def load(config):
 config=Path(config).resolve();a=read(config)
 stats=tuple((name,(config.parent/name).stat().st_ino,(config.parent/name).stat().st_size,(config.parent/name).stat().st_mtime_ns) for name in a['files'])
 return dataset(str(config),a['datasetFingerprint'],stats)
def selection(config,s):
 a,sources=load(config);require(s.get('specimen') in sources,'Unknown source specimen');d=sources[s['specimen']];src=d['source'];require(s.get('population')==src['populationID'],'Population does not match source');require(s.get('condition')==src['conditionID'],'Unsupported condition: select the recorded assay endpoint');return a,d

def page(q,maximum=128):
 limit=q.get('limit',128);offset=q.get('offset',0);require(type(limit)==int and 1<=limit<=maximum and type(offset)==int and offset>=0,'Bounded offset/limit required');return offset,limit

def accessible(d,target):return [i for i in d['groups'].get(target,[]) if d['rows'][i]['observationAccess']=='development-open']
@lru_cache(maxsize=16)
def verify_expression(files):
 for path,expected,*identity in files:require(sha(path)==expected,'Admitted expression artifact changed: '+path)
 return True

def values(d,rows,gene,normalization):
 verify_expression(d['expressionFiles'])
 if d['records'] is None:d['records']=np.memmap(d['countPath'],dtype=DTYPE,mode='r')
 ix=d['features'].get(gene,[]);require(len(ix)==1,'Feature absent or ambiguous on this source measurement axis');j=ix[0];out=[]
 for i in rows:
  require(d['rows'][i]['observationAccess']=='development-open','Reserved response access denied')
  v=d['records'][d['offsets'][i]:d['offsets'][i+1]];matches=v['count'][v['feature']==j];count=int(matches.sum()) if len(matches) else 0
  out.append(count if normalization=='raw_counts' else math.log1p(count/d['totals'][i]*10000))
 return out

def candidate(d,t):
 controls=accessible(d,'no-intervention');allrows=d['groups'][t];opened=accessible(d,t)
 return {'target':t,'role':'direct','controlCells':len(controls),'outcomeCells':len(allrows),'accessibleOutcomeCells':len(opened),'reservedOutcomeCells':len(allrows)-len(opened),'canExecute':False,'recommendationSupported':False,'biologicalUnits':None,'supportStatus':'MEASURED' if opened else 'RESERVED','access':'development-open' if len(opened)==len(allrows) else 'partially-reserved' if opened else 'reserved-development','reason':BLOCK}

class MeasuredCellularCohort:
 family='measured-cellular-cohort';adapter_id='numivivo-measured-cellular-cohort/v1'
 def catalog(self,config):
  a=read(config)
  return {k:a[k] for k in ('id','title','family','presentation','geometry','specimens','populations','conditions','targets','features','featuresTotal','evidence','limits','biologicalPromotion','datasetFingerprint','observationAccessFingerprint')}|{'capabilities':self.capabilities(config)}
 def capabilities(self,config):
  a=read(config)
  return {k:a[k] for k in ('presentation','geometry','specimens','populations','conditions','targets','features')}|{'prediction':False,'observationComparison':False,'measuredExploration':True,'readouts':['gene','interventions','features','populations','expression-tile','objective-coverage'],'preservationSupported':False,'predictionUnavailableReason':BLOCK}
 def validate_selection(self,config,s):
  selection(config,s);raise ValueError(BLOCK)
 def support(self,config,s,genes):
  a,d=selection(config,s);controls=accessible(d,'no-intervention')
  return {'population':s['population'],'populationAnchors':len(d['rows']),'referenceCells':len(controls),'reference':{'referenceCells':len(controls),'uniqueCells':len(controls)},'candidates':[candidate(d,t) for t in sorted(d['groups']) if t!='no-intervention'],'genes':[{'gene':g,'measured':len(d['features'].get(g,[]))==1,'modeled':False} for g in genes],'recommendationSupported':False,'reliableWinner':False,'reason':BLOCK,'scope':SCOPE,'preservationSupported':False}
 def inspect_objective(self,config,s,obj):
  a,d=selection(config,s)
  return {'genes':[{'gene':g,'eligible':len(d['features'].get(g,[]))==1,'measured':len(d['features'].get(g,[]))==1,'modeled':False,'referenceQuantity':'MEASURED' if len(d['features'].get(g,[]))==1 else 'UNAVAILABLE','predictionQuantity':'UNAVAILABLE','missingFeatureReason':None if len(d['features'].get(g,[]))==1 else 'Absent or ambiguous source feature; not a zero measurement','selectedPopulationReferenceCells':len(accessible(d,'no-intervention'))} for g in obj['genes']],'canExecute':False,'candidateEffects':[],'executionReason':BLOCK,'objective':'Unexecuted proposal; intent retained','biologicalPromotion':False}
 def readout(self,config,record,q,s):
  require(record is None,'Measured cohort has no sealed prediction records');a,d=selection(config,s);require(q.get('mode','measured')=='measured',BLOCK);offset,limit=page(q);kind=q.get('kind','gene');search=q.get('search','');require(isinstance(search,str) and len(search)<=128,'Bounded search required')
  common={'datasetFingerprint':a['datasetFingerprint'],'observationAccessFingerprint':a['observationAccessFingerprint'],'offset':offset,'limit':limit,'scope':SCOPE,'access':{'status':'development-open-only','reservedResponsesReturned':False},'presentation':'population','modeled':False,'revealed':False}
  if kind in ('interventions','features','populations'):
   if kind=='interventions':items=[candidate(d,t) for t in sorted(d['groups']) if t!='no-intervention' and search.lower() in t.lower()]
   elif kind=='features':items=sorted([{'id':f['name'] if len(d['features'][f['name']])==1 else f['id'],'title':f['name']+' · '+f['id']} for f in d['metadata']['features'] if search.lower() in f['name'].lower() or search.lower() in f['id'].lower()],key=lambda x:x['id'])
   else:items=[p for p in a['populations'] if search.lower() in p['title'].lower()]
   return common|{'items':items[offset:offset+limit],'total':len(items)}
  if kind=='objective-coverage':return common|self.inspect_objective(config,s,s['objective'])
  normalization=q.get('normalization','raw_counts');require(normalization in ('raw_counts','log1p_10000'),'Unsupported normalization');units='raw UMI counts' if normalization=='raw_counts' else 'log1p(count / retained-panel cell total × 10000)'
  target=q.get('target');require(target in d['groups'],'Select an explicit supported intervention for a bounded measured readout')
  opened=accessible(d,target);controls=accessible(d,'no-intervention');allrows=d['groups'][target]
  if kind=='expression-tile':
   genes=q.get('genes');require(isinstance(genes,list) and 1<=len(genes)<=32 and limit<=64,'Tile bound is 64 cells × 32 genes');rows=opened[offset:offset+limit];mask=[len(d['features'].get(g,[]))==1 for g in genes];columns=[values(d,rows,g,normalization) if ok else [None]*len(rows) for g,ok in zip(genes,mask)]
   return common|{'units':units,'genes':genes,'featureMask':mask,'cellIDs':[d['rows'][i]['id'] for i in rows],'values':list(map(list,zip(*columns))),'total':len(opened)}
  require(kind=='gene','Unsupported readout');gene=q.get('gene');require(isinstance(gene,str),'Gene required');measured=len(d['features'].get(gene,[]))==1
  missing=None if measured else 'Feature absent or ambiguous on this source panel; unmeasured features are not zero counts'
  # Means are computed over accessible cells only, never over reserved outcomes.
  cv=values(d,controls,gene,normalization) if measured else [];ov=values(d,opened,gene,normalization) if measured else []
  arm={'target':target,'role':'direct','population':s['population'],'control':float(np.mean(cv)) if cv else None,'observed':float(np.mean(ov)) if ov else None,'predicted':None,'residual':None,'variance':None,'controlCells':len(controls),'outcomeCells':len(allrows),'accessibleOutcomeCells':len(opened),'reservedOutcomeCells':len(allrows)-len(opened),'controlValues':cv[offset:offset+limit],'observedValues':ov[offset:offset+limit],'controlCellIDs':[d['rows'][i]['id'] for i in controls[offset:offset+limit]] if measured else [],'observedCellIDs':[d['rows'][i]['id'] for i in opened[offset:offset+limit]] if measured else [],'evidence':'MEASURED' if ov else 'UNAVAILABLE','biologicalUnits':None,'access':candidate(d,target)['access'],'meanScope':'Accessible selected-population cells only; reserved responses excluded'}
  return common|{'gene':gene,'units':units,'normalization':normalization,'measured':measured,'missingFeatureReason':missing,'arms':[arm],'total':max(len(controls),len(opened)),'distributionStatus':'Individual measured cells; bounded page. Cell variability is not uncertainty about biological transfer.'}
 def predict(self,*args,**kwargs):raise ValueError(BLOCK)
 def reveal(self,*args,**kwargs):raise ValueError('Reserved cohort outcomes require a separately frozen evaluation campaign; exploration cannot reveal them')
 def verify(self,*args,**kwargs):raise ValueError('No prediction to replay; use native corpus verification')
 def summary(self,*args,**kwargs):raise ValueError('No experiment record for measured-only cohort')
 def evaluate_objective(self,*args,**kwargs):raise ValueError(BLOCK)
