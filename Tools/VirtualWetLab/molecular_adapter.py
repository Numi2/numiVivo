"""Gene perturbation on a real spatial specimen, using the existing native RNA owner.

This does not implement a second predictor or infer a binding/pathway mechanism.
The measured barcode-to-knockout transition remains a hypothesis.
"""
from pathlib import Path
import math, shutil, uuid
import wetlab as rna
from adapters import WetLabExperimentAdapter

FORMAT='numivivo-molecular-spatial-assay/v1'


def load(path):
    path=Path(path).resolve();a=rna.read(path)
    rna.require(a['format']==FORMAT and a['family']=='molecular-perturbation','Molecular assay schema')
    for name,digest in a['files'].items():
        p=path.parent/name
        rna.require(p.resolve().is_relative_to(path.parent) and rna.sha(p)==digest,'Changed source: '+name)
    return a,path.parent


def identity(a):
    return {**rna.runtime_identity(a['runtime']['binary']),'tissueBinary':rna.sha(a['runtime']['tissueBinary']),
            'molecularAdapter':rna.sha(__file__)}


def check(run,exact=False):
    run=Path(run);seal=rna.read(run/'seal.json');rna.require(seal['format']==FORMAT,'Molecular seal schema')
    for name,digest in seal['files'].items():
        p=run/name;rna.require(p.resolve().is_relative_to(run.resolve()) and not p.is_symlink() and rna.sha(p)==digest,'Changed sealed artifact: '+name)
    reg=rna.read(run/'registration.json')
    if exact:rna.require(identity(reg['assay'])==reg['runtimeIdentity'],'Exact molecular runtime changed')
    return reg


def predictions(run):
    run=Path(run);reg=check(run);rows=[];features=None
    for i,region in enumerate(reg['plan']['regionIDs']):
        p=rna.read(run/'prediction/folds'/f'{i:03d}'/'prediction.json');features=p['featureIDs']
        rows.append({'regionID':region,**p['predictions'][0]})
    return features,rows


def score(run):
    run=Path(run);reg=check(run);features,rows=predictions(run)
    bulk=rna.read(run/'prediction/source/report.json')['pseudobulk']
    m=bulk['matrix'];regions=[];geneSquared={k:[0.0]*len(features) for k in ('contextRidge','noChange','meanResponse','medianResponse')}
    definitions={r['id']:r for r in reg['assay']['regions']}
    for p in rows:
        d=definitions[p['regionID']];j=d['observedRow'];g=bulk['groups'][j]
        rna.require(g['donorID']==reg['specimen'] and g['condition']==reg['assay']['intervention']['barcode'],'Held-out observation identity')
        counts=[0]*len(features)
        for k in range(m['rowOffsets'][j],m['rowOffsets'][j+1]):counts[m['featureIndices'][k]]=m['counts'][k]
        total=sum(counts);rna.require(total>0,'Empty observation library');observed=[math.log1p(v/total*1e6) for v in counts]
        metrics=[];errors=None
        for e in p['estimates']:
            residual=[x-y for x,y in zip(e['predictedTreated'],observed)]
            metrics.append({'model':e['baseline'],'rmse':math.sqrt(math.fsum(x*x for x in residual)/len(residual)),'mae':math.fsum(map(abs,residual))/len(residual)})
            for z,v in enumerate(residual):geneSquared[e['baseline']][z]+=v*v
            if e['baseline']=='contextRidge':errors=residual
        by={m['model']:m['rmse'] for m in metrics}
        regions.append({'regionID':p['regionID'],'observed':observed,'residual':errors,'metrics':metrics,'cellCount':len(g['sourceCellIndices']),
                        'criterionPassed':all(by['contextRidge']<by[b] for b in ('noChange','meanResponse'))})
    aggregate=[{'model':k,'rmse':math.sqrt(math.fsum(v)/(len(features)*len(rows)))} for k,v in geneSquared.items()]
    errors=geneSquared['contextRidge'];order=sorted(range(len(features)),key=lambda i:(-errors[i],features[i]))
    return {'format':FORMAT,'predictionSealSHA256':rna.sha(run/'seal.json'),'units':'log1p(CPM over common gene universe)',
        'regions':regions,'metrics':aggregate,'geneRMSE':[math.sqrt(v/len(rows)) for v in errors],
        'largestGeneErrors':[{'gene':features[i],'rmse':math.sqrt(errors[i]/len(rows))} for i in order[:30]],
        'failureClusters':[{'regionID':r['regionID'],'reason':'does not improve over both baselines'} for r in regions if not r['criterionPassed']],
        'verdict':'criterion-met-on-development-specimen' if all(r['criterionPassed'] for r in regions) else 'criterion-contradicted',
        'calibration':{'status':'UNAVAILABLE','reason':'Point predictor has no calibrated predictive distribution'},
        'qualification':'retrospective held-out mouse, known perturbation; no zero-shot/causal/biological qualification',
        'biologicalUnit':reg['specimen'],'biologicalReplicates':1,'regionsAreReplicates':False}


class MolecularPerturbationAdapter(WetLabExperimentAdapter):
    adapter_id='molecular-perturbation/v1';family='molecular-perturbation'
    def catalog(self,config):
        a,root=load(config)
        return {'id':a['id'],'adapterID':self.adapter_id,'family':self.family,'title':a['title'],'intervention':a['intervention'],
            'specimens':[{'id':a['specimen'],'regions':a['regions']}],'timepoints':a['timepoints'],'limits':a['limits'],
            'sourceCitation':a['sourceCitation'],'provenance':'GSE274447 · three source-labelled mice · held-out mouse 3',
            'excludedRegions':a['excludedRegions'],'design':rna.read(root/'design-specimen.json')}
    def compile(self,config,selection):
        a,root=load(config)
        rna.require(set(selection)<= {'specimen','target','regionIDs','timepoint'},'Unsupported plan fields; no silent partial execution')
        rna.require(selection.get('specimen')==a['specimen'],'Unsupported specimen')
        rna.require(selection.get('target')==a['intervention']['target'],'Unsupported target; binding and pathway transfer not available')
        rna.require(selection.get('timepoint') in a['timepoints'],'Unsupported time; no temporal interpolation')
        ids=selection.get('regionIDs');supported={r['id'] for r in a['regions']}
        rna.require(isinstance(ids,list) and ids and len(set(ids))==len(ids) and set(ids)<=supported,'Select supported regions')
        return {'format':'numivivo-experiment-plan/v2','specimenID':a['specimen'],'assayID':a['id'],'adapterID':self.adapter_id,
            'intervention':a['intervention'],'regionIDs':sorted(ids),'timepoint':selection['timepoint'],'control':'safe-harbour barcode-positive cells in each selected region',
            'transitions':[{'from':'measured guide barcode','to':'Clu knockout','state':'HYPOTHESIS','provenance':a['sourceCitation'],'modelID':'unique-positive-guide-v1','uncertainty':'per-cell functional knockout unverified','validation':'not-established'},
                {'from':'training mouse control/knockout RNA','to':'held-out regional RNA','state':'MODEL INFERENCE','provenance':a['files']['source/report.json'],'modelID':'native-context-ridge-alpha1','uncertainty':'uncalibrated','validation':'development evaluation pending'},
                {'from':'gene knockout','to':'binding/reaction/pathway kinetics','state':'UNAVAILABLE','provenance':a['sourceCitation'],'modelID':None,'uncertainty':'unavailable','validation':'not-implemented'},
                {'from':'RNA change','to':'growth/force/phenotype','state':'UNAVAILABLE','provenance':a['sourceCitation'],'modelID':None,'uncertainty':'unavailable','validation':'not-implemented'}]}
    def predict(self,config,runtime,workspace,selection):
        a,root=load(config);plan=self.compile(config,selection);rid=identity(a)
        workspace=Path(workspace);workspace.mkdir(parents=True,exist_ok=True);run=workspace/uuid.uuid4().hex;run.mkdir()
        reg={'format':FORMAT,'adapterID':self.adapter_id,'family':self.family,'assay':a,'specimen':a['specimen'],'createdAt':rna.timestamp(),
             'plan':plan,'runtimeIdentity':rid,'purpose':'retrospective development qualification','limits':a['limits']}
        rna.write(run/'registration.json',reg)
        for name in ('specimen.json','design-specimen.json','protocol.json','lineage.json'):shutil.copy2(root/name,run/name)
        shutil.copy2(__file__,run/'molecular_adapter.py')
        definitions={r['id']:r for r in a['regions']}
        folds=[{'id':region,'perturbationID':'Clu-KO','controlCondition':a['intervention']['control'],'treatmentCondition':a['intervention']['barcode'],
            'cellGroup':'segmented-brain-cells','trainingGroupIndices':a['trainingRows'],'queryGroupIndices':[definitions[region]['controlRow']]} for region in plan['regionIDs']]
        rna.write(run/'native-plan.json',{'schemaVersion':1,'sourceReport':{'bytes':list(bytes.fromhex(a['files']['source/report.json']))},'featureNamespace':'MOUSE_GENE_SYMBOL','provenance':'Spatial regions within held-out mouse; not independent replicates. Native donor response applied to regional control reference.','folds':folds})
        try:
            native=a['runtime']['binary']
            rna.native(native,['singlecell-perturbation-batch',root/'source','--plan',run/'native-plan.json','--output',run/'prediction'],run/'prediction-log.json')
            receipt=rna.read(run/'prediction/receipt.json');rna.require(all(f['status']=='completed' for f in receipt['folds']),'Native fold failed')
            rna.write(run/'seal.json',{'format':FORMAT,'sealedAt':rna.timestamp(),'files':rna.inventory(run)})
        except Exception as e:rna.write(run/'failure.json',{'error':str(e),'at':rna.timestamp()});raise
        return run
    def reveal(self,run,runtime):
        run=Path(run);reg=check(run,True);rna.require(not (run/'comparison.json').exists(),'Already revealed')
        self.verify(run,runtime);result=score(run);rna.write(run/'comparison.json',result);return result
    def verify(self,run,runtime):
        run=Path(run);reg=check(run,True)
        rna.native(reg['assay']['runtime']['binary'],['singlecell-perturbation-batch-verify',run/'prediction'],run/('replay-'+uuid.uuid4().hex+'.json'))
        rna.native(reg['assay']['runtime']['tissueBinary'],[run/'specimen.json'],run/('tissue-replay-'+uuid.uuid4().hex+'.json'))
        if (run/'comparison.json').exists():rna.require(rna.read(run/'comparison.json')==score(run),'Comparison does not reconstruct')
        return {'status':'verified','run':run.name,'family':self.family}
    def summary(self,run):
        run=Path(run);reg=check(run);features,rows=predictions(run);revealed=(run/'comparison.json').exists()
        result={'id':run.name,'family':self.family,'adapterID':self.adapter_id,'registration':reg,'recordDirectory':str(run.resolve()),'featureIDs':features,
            'regionalPredictions':rows,'revealed':revealed,'observationsAvailable':True,'design':rna.read(run/('specimen.json' if revealed else 'design-specimen.json'))}
        if revealed:
            result['comparison']=rna.read(run/'comparison.json');rna.require(result['comparison']==score(run),'Comparison changed')
        return result
