"""Intervention-choice operations on the existing experiment campaign records.

No second simulator: all candidate distributions are sealed native model outputs.
This module compiles a molecular objective and evaluates the choice itself.
"""
from pathlib import Path
import uuid, shutil, subprocess
import numpy as np
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp,inventory
from design_objective import objective,utility,selection_score,acquisition_backtest
from prepare_intervention_design import load_source,axis,orthology,normalized,V03


def bound(campaign):
    campaign=Path(campaign);seal=read(campaign/'prediction-seal.json');reg=read(campaign/'registration.json')
    for f,h in seal['files'].items():require(sha(campaign/f)==h,'Changed trained campaign artifact: '+f)
    inputs=Path(reg['inputs'])
    require(sha(inputs/'prepared.json')==reg['preparedSHA256'],'Prepared cohort changed')
    prepared=read(inputs/'prepared.json')
    for f,h in prepared['files'].items():require(sha(inputs/f)==h,'Changed input: '+f)
    return campaign,inputs,reg


def catalog(campaign):
    campaign,inputs,reg=bound(campaign);rows=read(inputs/'reserved-rows.json');contexts=[]
    for source,context in sorted({(r['source'],r['context']) for r in rows}):
        selected=[r for r in rows if r['source']==source and r['context']==context]
        contexts.append({'id':source+'|'+context,'source':source,'context':context,'targets':[r['target'] for r in selected],
          'controlCells':selected[0]['controlCellCount'],'unit':selected[0]['unit'],'modality':selected[0]['modality'],
          'receiverSupport':'UNAVAILABLE: dissociated population has no measured spatial relationships',
          'geometry':'UNAVAILABLE: no invented cell positions','biologicalUnits':1})
    return {'contexts':contexts,'features':read(inputs/'features.json'),'objectives':read(inputs/'protocol.json')['objectives'],
            'models':list(load_file(str(campaign/'all-predictions.safetensors'))),
            'evidence':'MODEL INFERENCE','reliableWinner':False,'biologicalPromotion':False,
            'limits':['RNA research utility, not treatment efficacy, safety or survival','No qualified transfer uncertainty; candidate ordering is experimental',
              'New mouse batch and Jurkat study reserved once; later workspace comparisons are retrospective',
              'GO-neighbor transfer is a transparent baseline, not the GEARS neural model; prior State comparator is PerturbMean, not State Transition']}


def preview(campaign,selection):
    campaign,inputs,reg=bound(campaign);require(set(selection)=={'context','objective'},'Context and objective required')
    features=read(inputs/'features.json');indices=objective(features,selection['objective']);require(not indices[1],'Preservation requires supported receiving-population measurements; unavailable for dissociated cohorts')
    rows=read(inputs/'reserved-rows.json');ix=[i for i,r in enumerate(rows) if r['source']+'|'+r['context']==selection['context']]
    require(len(ix)>=2,'Fewer than two supported interventions');query=load_file(str(inputs/'reserved.safetensors'));pred=load_file(str(campaign/'all-predictions.safetensors'))
    candidates=[rows[i]['target'] for i in ix]+['no-intervention'];require(len(set(candidates))==len(candidates),'Duplicate target in comparison context')
    utilities={};effects={}
    for model,value in pred.items():
        delta=value[ix]-query['context'][ix];u,benefit,outside=utility(delta,None,indices)
        utilities[model]=u.tolist()+[0.];effects[model]=delta[:,indices[0]].tolist()+[[0.]*len(indices[0])]
    order=sorted(range(len(candidates)),key=lambda i:(-utilities['target-descriptor'][i],candidates[i]))
    arms=[{'target':candidates[i],'predictedUtility':utilities['target-descriptor'][i],
      'programDelta':-utilities['target-descriptor'][i],'outsidePopulationChange':None,
      'trainingMeanUtility':utilities['training-mean'][i],'modelRange':[min(utilities[m][i] for m in ('target-descriptor','target-ID','shuffled-target-descriptor')),max(utilities[m][i] for m in ('target-descriptor','target-ID','shuffled-target-descriptor'))],
      'contributingControls':rows[ix[0]]['controlCellCount'],'contributingTreated':rows[ix[i]]['outcomeCellCount'] if i<len(ix) else rows[ix[0]]['controlCellCount'],
      'evidence':'HYPOTHESIS' if i<len(ix) else 'CONTROL REFERENCE','support':'unqualified transfer uncertainty; no reliable winner'} for i in order]
    native=load_file(str(campaign/'target-descriptor/reserved-prediction/prediction.safetensors'))
    distributions={rows[i]['target']:{'control':query['context'][i,indices[0]].tolist(),'predicted':native['mean'][i,indices[0]].tolist(),'cellVariance':native['variance'][i,indices[0]].tolist()} for i in ix}
    return {'selection':selection,'rows':ix,'candidates':candidates,'utilities':utilities,'geneEffects':effects,'rankedCandidates':arms,'distributions':distributions,
      'numericalChoice':candidates[order[0]],'reliableWinner':False,'conclusion':'The model does not support a reliable winner. This is an experimental numerical ranking.',
      'objectiveEquation':'utility = -mean(log1p(CPM) response over selected genes); no intervention = 0',
      'sourcePredictionSealSHA256':sha(campaign/'prediction-seal.json')}


def seal(campaign,workspace,selection):
    result=preview(campaign,selection);root=Path(workspace)/('campaign-'+uuid.uuid4().hex);root.mkdir(parents=True)
    campaign,inputs,reg=bound(campaign)
    write(root/'registration.json',{'format':'numivivo-experiment-campaign/v1','kind':'intervention-design','createdAt':timestamp(),
      'campaign':str(campaign.resolve()),'selection':selection,'ownerSHA256':sha(__file__),'objectiveOwnerSHA256':sha(Path(__file__).with_name('design_objective.py')),
      'biologicalUnitSemantics':'Source batch or study, not cells or captures','revealPolicy':'All candidate predictions and exact objective sealed together before observations'})
    write(root/'prediction.json',result)
    write(root/'prediction-seal.json',{'registrationSHA256':sha(root/'registration.json'),'predictionSHA256':sha(root/'prediction.json'),'sealedAt':timestamp()})
    return {'id':root.name,'prediction':result,'revealed':False}


def check(root):
    root=Path(root);s=read(root/'prediction-seal.json');reg=read(root/'registration.json')
    require(s['registrationSHA256']==sha(root/'registration.json') and s['predictionSHA256']==sha(root/'prediction.json'),'Changed sealed objective or candidate predictions')
    require(reg['ownerSHA256']==sha(__file__) and reg['objectiveOwnerSHA256']==sha(Path(__file__).with_name('design_objective.py')),'Exact campaign owner changed')
    require(read(root/'prediction.json')==preview(reg['campaign'],reg['selection']),'Candidate prediction reconstruction differs')
    return reg


def observations(campaign,indices):
    campaign,inputs,reg=bound(campaign);rows=read(inputs/'reserved-rows.json');specs=read(inputs/'prepared.json')['sources'];features=read(inputs/'features.json');mapping=orthology(V03/'research/MGI-HOM_MouseHumanSequence.rpt')
    values={};detail={}
    for sid in sorted({rows[i]['source'] for i in indices}):
        s=next(s for s in specs if s['id']==sid);require(sha(s['path'])==s['sha256'],'Raw source changed');a=load_source(s);ax=axis(a,mapping,s['species']);cols=[ax[g] for g in features]
        for i in indices:
            r=rows[i]
            if r['source']!=sid:continue
            y=normalized(a,r['rows'],cols);values[i]=y.mean(0);detail[i]={'mean':y.mean(0),'variance':y.var(0),'count':len(y)}
        del a
    return np.asarray([values[i] for i in indices]),detail


def reveal(root):
    root=Path(root);reg=check(root);require(not (root/'comparison.json').exists(),'Already revealed')
    p=read(root/'prediction.json');campaign,inputs,_=bound(reg['campaign']);indices=p['rows'];q=load_file(str(inputs/'reserved.safetensors'));features=read(inputs/'features.json');ox=objective(features,reg['selection']['objective'])
    # Durable observation snapshot permits interrupted analysis recovery without
    # re-extraction or silently changing measured inputs.
    snapshot=root/'observations.safetensors'
    if not snapshot.exists():
        y,detail=observations(campaign,indices);save_file({'mean':y,'variance':np.asarray([detail[i]['variance'] for i in indices]),'count':np.asarray([detail[i]['count'] for i in indices],np.int32)},str(snapshot))
        write(root/'observation-source.json',{'sourceSealSHA256':sha(campaign/'prediction-seal.json'),'observationSHA256':sha(snapshot),'openedAt':timestamp()})
    else:require(sha(snapshot)==read(root/'observation-source.json')['observationSHA256'],'Interrupted observation snapshot changed')
    snapshot_data=load_file(str(snapshot));y=snapshot_data['mean'];u,_,_=utility(y-q['context'][indices],None,ox);observed=u.tolist()+[0.]
    results={model:selection_score(p['candidates'],score,observed,p['utilities']['training-mean']) for model,score in p['utilities'].items()}
    pred=load_file(str(campaign/'all-predictions.safetensors'));all_metrics={}
    for model,value in pred.items():
        err=value[indices]-y;response=y-q['context'][indices];pd=value[indices]-q['context'][indices];mask=np.abs(response)>.1
        all_metrics[model]={'endpointRMSE':float(np.sqrt(np.mean(err**2))),'responseRMSE':float(np.sqrt(np.mean((pd-response)**2))),
           'responseDirection':float(np.mean(np.sign(pd[mask])==np.sign(response[mask]))) if mask.any() else None,
           'magnitudeMAE':float(np.mean(np.abs(np.abs(pd)-np.abs(response))))}
    acquisition=acquisition_backtest([p['utilities'][m][:-1] for m in ('target-descriptor','target-ID','shuffled-target-descriptor')],observed[:-1],[r['contributingTreated'] for r in sorted(p['rankedCandidates'],key=lambda a:p['candidates'].index(a['target']))[:-1]])
    measured={target:{'mean':y[j,ox[0]].tolist(),'cellVariance':snapshot_data['variance'][j,ox[0]].tolist(),'count':int(snapshot_data['count'][j])} for j,target in enumerate(p['candidates'][:-1])}
    result={'format':'numivivo-campaign-comparison/v1','observedUtilities':dict(zip(p['candidates'],observed)),'selectionResults':results,'metrics':all_metrics,'distributions':measured,
      'acquisition':acquisition,'biologicalPromotion':False,'reliableWinner':False,'verdict':'experimental-selection-evaluated; biological promotion unavailable',
      'uncertainty':'One batch/study; cell spread and technical captures cannot qualify transfer confidence',
      'predictionSealSHA256':sha(root/'prediction-seal.json'),'observationSHA256':sha(snapshot)}
    write(root/'comparison.json',result);write(root/'observation-seal.json',{'comparisonSHA256':sha(root/'comparison.json'),'predictionSealSHA256':sha(root/'prediction-seal.json')})
    return {'id':root.name,'prediction':p,'comparison':result,'revealed':True}


def verify(root):
    root=Path(root);reg=check(root);campaign,inputs,binding=bound(reg['campaign']);require(sha(binding['binary'])==binding['binarySHA256'],'Native runtime changed')
    from train_intervention_design import native,VARIANTS
    checks={}
    for variant in VARIANTS:
        folder=campaign/variant;step=read(folder/'selection.json')['selected']['step'];dest=root/('replay-'+variant+'-'+uuid.uuid4().hex)
        native(binding['binary'],'predict',folder/'plan.json',folder/'reserved.safetensors',dest,folder/'training'/f'weights-{step}.safetensors')
        original=load_file(str(folder/'reserved-prediction/prediction.safetensors'));actual=load_file(str(dest/'prediction.safetensors'))
        require(all(np.array_equal(v,actual[k]) for k,v in original.items()),'Native replay mismatch');checks[variant]='bit-exact mean and variance'
    if (root/'comparison.json').exists():
        s=read(root/'observation-seal.json');require(s['comparisonSHA256']==sha(root/'comparison.json'),'Changed selection evaluation')
        p=read(root/'prediction.json');c=read(root/'comparison.json');observed=[c['observedUtilities'][t] for t in p['candidates']]
        for model,values in p['utilities'].items():require(selection_score(p['candidates'],values,observed,p['utilities']['training-mean'])==c['selectionResults'][model],'Selection evaluation replay mismatch')
        require(c['observationSHA256']==sha(root/'observations.safetensors'),'Changed observations')
    return {'status':'verified','id':root.name,'nativeReplay':checks,'biologicalPromotion':False}
