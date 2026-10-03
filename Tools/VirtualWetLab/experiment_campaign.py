"""ExperimentCampaign recorder over assay owners; no second simulation scheduler."""
from pathlib import Path
import itertools, uuid
import laboratory
from wetlab import read,write,require,sha,timestamp


def expand(matrix):
    require(isinstance(matrix,dict) and matrix and all(isinstance(v,list) and v for v in matrix.values()),'Nonempty intervention matrix required')
    keys=sorted(matrix);size=1
    for k in keys:size*=len(matrix[k])
    require(size<=32,'Campaign limited to 32 declared arms')
    return [dict(zip(keys,values)) for values in itertools.product(*(matrix[k] for k in keys))]


def register_and_predict(config,runtime,workspace,request):
    require(set(request)=={'hypothesis','matrix','successCriterion','simulationReplicates'},'Campaign contract fields')
    require(isinstance(request['hypothesis'],str) and request['hypothesis'].strip(),'Hypothesis required')
    require(request['successCriterion']=='adapter-primary-criterion-in-every-arm','Unsupported success criterion')
    require(type(request['simulationReplicates']) is int and 1<=request['simulationReplicates']<=3,'Simulation replicate budget')
    selections=expand(request['matrix']);adapter=laboratory.adapter_for_config(config)
    require(hasattr(adapter,'compile'),'This campaign version requires a typed-plan adapter')
    plans=[adapter.compile(config,s) for s in selections]
    root=Path(workspace)/('campaign-'+uuid.uuid4().hex);root.mkdir(parents=True)
    write(root/'registration.json',{'format':'numivivo-experiment-campaign/v1','createdAt':timestamp(),'request':request,'plans':plans,
        'assaySHA256':sha(config),'biologicalUnits':sorted({p['specimenID'] for p in plans}),
        'replicateSemantics':'Deterministic simulation repetitions are replay checks, never biological replicates.',
        'revealPolicy':'All declared predictions must be sealed before any campaign reveal.'})
    arms=[]
    try:
        for i,s in enumerate(selections):
            for replicate in range(request['simulationReplicates']):
                run=adapter.predict(config,runtime,root/'arms',s)
                entry={'selection':s,'replicate':replicate,'run':run.name,'sealSHA256':sha(run/'seal.json')}
                arms.append(entry);write(root/f'arm-{len(arms):03d}.json',entry)
        write(root/'prediction-seal.json',{'registrationSHA256':sha(root/'registration.json'),'arms':arms,'sealedAt':timestamp()})
    except Exception as error:
        write(root/'failure.json',{'error':str(error),'retainedArms':arms,'at':timestamp()});raise
    return root


def verify_barrier(root):
    root=Path(root);seal=read(root/'prediction-seal.json');reg=read(root/'registration.json')
    require(sha(root/'registration.json')==seal['registrationSHA256'],'Campaign registration changed')
    require(len(seal['arms'])==len(reg['plans'])*reg['request']['simulationReplicates'],'Incomplete campaign')
    require(len({a['run'] for a in seal['arms']})==len(seal['arms']),'Duplicate campaign run')
    selections=expand(reg['request']['matrix'])
    for i,entry in enumerate(seal['arms']):
        require(entry['selection']==selections[i//reg['request']['simulationReplicates']] and entry['replicate']==i%reg['request']['simulationReplicates'],'Campaign arm differs from preregistered matrix')
        require(len(entry['run'])==32 and all(c in '0123456789abcdef' for c in entry['run']),'Invalid campaign path')
        run=root/'arms'/entry['run'];require(sha(run/'seal.json')==entry['sealSHA256'],'Campaign prediction seal changed')
        summary=laboratory.adapter_for_run(run).summary(run)
        require(summary['registration']['plan']==reg['plans'][i//reg['request']['simulationReplicates']],'Campaign arm plan changed')
    return seal


def reveal(root,runtime):
    root=Path(root);seal=verify_barrier(root);require(not (root/'comparison.json').exists(),'Campaign already compared')
    rows=[]
    for arm in seal['arms']:
        run=root/'arms'/arm['run'];adapter=laboratory.adapter_for_run(run)
        # Crash recovery resumes only already sealed/revealed arms, never refits.
        if not (run/'comparison.json').exists():adapter.reveal(run,runtime)
        result=adapter.summary(run)['comparison'];rows.append({'run':run.name,'verdict':result['verdict'],'comparisonSHA256':sha(run/'comparison.json')})
    result={'format':'numivivo-campaign-comparison/v1','predictionSealSHA256':sha(root/'prediction-seal.json'),'arms':rows,
            'criterionPassed':all(r['verdict']=='criterion-met-on-development-specimen' for r in rows),
            'biologicalQualification':'not-established'}
    write(root/'comparison.json',result);return result


def verify(root,runtime):
    root=Path(root);seal=verify_barrier(root)
    for arm in seal['arms']:
        run=root/'arms'/arm['run'];laboratory.adapter_for_run(run).verify(run,runtime)
    if (root/'comparison.json').exists():
        c=read(root/'comparison.json');require(c['predictionSealSHA256']==sha(root/'prediction-seal.json'),'Campaign comparison binding')
        require([r['run'] for r in c['arms']]==[a['run'] for a in seal['arms']],'Campaign comparison membership')
        for row in c['arms']:
            result=read(root/'arms'/row['run']/'comparison.json')
            require(row['comparisonSHA256']==sha(root/'arms'/row['run']/'comparison.json'),'Campaign arm comparison changed')
            require(row['verdict']==result['verdict'],'Campaign arm verdict changed')
        require(c['criterionPassed']==all(read(root/'arms'/r['run']/'comparison.json')['verdict']=='criterion-met-on-development-specimen' for r in c['arms']),'Campaign verdict changed')
    return {'status':'verified','campaign':root.name,'arms':len(seal['arms'])}


def propose_next(summary):
    """Uncalibrated model disagreement heuristic. Does not see held-out outcomes."""
    require(summary['family']=='molecular-perturbation','No proposal policy for this assay')
    ranked=[]
    for region in summary['regionalPredictions']:
        models={e['baseline']:e['predictedTreated'] for e in region['estimates']}
        disagreement=sum(abs(a-b) for a,b in zip(models['contextRidge'],models['meanResponse']))/len(models['contextRidge'])
        ranked.append({'regionID':region['regionID'],'score':disagreement})
    return {'policy':'unrevealed-model-disagreement-v1','evidenceState':'HYPOTHESIS','candidates':sorted(ranked,key=lambda x:(-x['score'],x['regionID'])),
            'limits':'Ranks where models disagree; not calibrated uncertainty or demonstrated information gain. No automatic real-world experiment execution.'}
