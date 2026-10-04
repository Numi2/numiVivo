"""Read-only objective support on existing assays; never opens perturbed RNA."""
from pathlib import Path
import numpy as np
from wetlab import read,require


def inspect(config, selection, genes):
    from learned_spatial import load, source
    a,root=load(config)
    require(genes and len(genes)==len(set(genes)),'Distinct objective genes required')
    require(selection['population'] in {p['id'] for p in a['populations']},'Unknown population')
    s=next(x for x in a['specimens'] if x['id']==selection['specimen'])
    data=source(s['source'],s['sourceSHA256']);features=read(root/'features.json');mask=np.load(root/'model-feature-mask.npy')
    refs=np.flatnonzero(data.obs.reference_admitted.to_numpy());rows=[]
    query=read(root/'query-rows.json');matched=sorted({j for r in query if r['cellType']+'|'+str(r['neighborhood'])==selection['population'] for j in r['referenceRows']})
    for gene in genes:
        j=features.index(gene) if gene in features else None
        measured=j is not None and bool(data.var.measured_in_source.iloc[j]);modeled=measured and bool(mask[j])
        counts=data.X[refs,j].toarray().ravel() if measured else np.array([])
        rows.append({'gene':gene,'referenceQuantity':'MEASURED' if measured else 'UNAVAILABLE','predictionQuantity':'MODEL INFERENCE' if modeled else 'UNAVAILABLE','referenceScope':'all admitted specimen reference cells, not independent units','selectedPopulationReferenceCells':len(matched),'referenceCells':len(refs) if measured else 0,'detectedReferenceCells':int((counts>0).sum()),'measurementUnits':'raw UMI','populationReference':'control-projected population; annotation uncertain','eligible':bool(modeled)})
    from safetensors.numpy import load_file
    landscape=read(root/'landscape.json');effects=load_file(str(root/'landscape.safetensors'))['effect'];js=[features.index(g) for g in genes if g in features]
    comparisons=[{'target':r['target'],'role':r['role'],'population':r['population'],'predictedUtility':-float(effects[i,js].mean()),'reference':r['reference'],'evidence':'MODEL INFERENCE'} for i,r in enumerate(landscape) if r['population']==selection['population'] and len(js)==len(genes) and all(x['eligible'] for x in rows)]
    return {'candidateEffects':comparisons,'previewReference':'Prior campaign anchor-matched reference; execution recomputes the pooled-reference comparison. Not a sealed ranking.','genes':rows,'canExecute':all(r['eligible'] for r in rows),'objective':'reduce mean log1p(CPM) across exactly these genes','preservationSupported':False,'preservationReason':'Verified sections and receiving-population controls unavailable','biologicalPromotion':False,'qualification':'Development only; measured feature availability is not adequate biological validation','sourceSHA256':s['sourceSHA256']}


def evaluate(config,run,genes):
    """Score exact sealed objective using the existing owner's marginal readouts."""
    from learned_spatial import LearnedSpatialResponseAdapter,check
    from design_objective import selection_score
    run=Path(run);reg=check(run);owner=LearnedSpatialResponseAdapter();fields=[owner.feature(config,reg['specimen'],g,run) for g in genes]
    arms=[]
    for i,arm in enumerate(fields[0]['arms']):
        rows=[f['arms'][i] for f in fields];available=all(x['predicted'] is not None for x in rows)
        predicted=-float(np.mean([x['predicted']-x['control'] for x in rows])) if available else None
        observed=-float(np.mean([x['observed']-x['control'] for x in rows])) if all(x['observed'] is not None for x in rows) else None
        arms.append({'target':arm['target'],'role':arm['role'],'population':reg['plan']['population'],'predictedUtility':predicted,'observedUtility':observed,'genes':[{'gene':g,'control':x['control'],'predicted':x['predicted'],'observed':x['observed'],'residual':x['residual']} for g,x in zip(genes,rows)]})
    direct=[a for a in arms if a['role']=='direct' and a['predictedUtility'] is not None];candidates=[a['target'] for a in direct]+['no-intervention'];pred=[a['predictedUtility'] for a in direct]+[0.]
    scores=None
    if all(a['observedUtility'] is not None for a in direct) and direct:
        measured=[a['observedUtility'] for a in direct]+[0.]
        scores=selection_score(candidates,pred,measured,[0.]*len(candidates))
        scores['trainingMeanRanking']='UNAVAILABLE in this spatial card; zero-change ranking is not a fitted training baseline'
        scores.pop('gainOverTrainingMean',None)
    return {'objectiveGenes':genes,'arms':arms,'numericalChoice':candidates[int(np.argmax(pred))],'scores':scores,'reliableWinner':False,'scope':'Exposed spatial development; one mouse and inadequate isolated controls; no independent biological confidence','preservation':'UNAVAILABLE'}
