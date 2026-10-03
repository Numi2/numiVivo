"""Spatial assay orchestration over native Matter geometry and Vivo transport."""
from pathlib import Path
import json
import math
import shutil
import tempfile
import uuid
from adapters import WetLabExperimentAdapter
from wetlab import read, write, sha, require, native, inventory, timestamp

FORMAT='numivivo-wet-lab-spatial/v1'
LIMITS=('Synthetic deformable scaffold and hypothetical passive cell exchange. '
        'Discrete graph transport with one-way accepted geometry transfer. '
        'No biological calibration, phenotype prediction or two-way mechanical feedback.')


def load(config):
    config=Path(config).resolve();a=read(config)
    require(a['format']==FORMAT and a['family']=='spatial-tissue','Unsupported spatial assay')
    require(a['specimens'] and len(a['specimens'])<=32,'Specimen count')
    require(len({s['id'] for s in a['specimens']})==len(a['specimens']),'Duplicate specimen')
    def path(name): return (config.parent/name).resolve()
    for s in a['specimens']:
        for key in ('geometry','mesh','material','plan'):
            require(sha(path(s[key]))==s[key+'SHA256'],'Spatial specimen input changed: '+key)
        p=read(path(s['plan']));g=read(path(s['geometry']))
        require(p['specimenID']==s['id'] and g['status']=='accepted','Specimen or accepted state mismatch')
    return a,path


def identity(a):
    names=('spatialBinary','partitionShader','matterBinary','matterShader')
    result={n:sha(a['runtime'][n]) for n in names}
    result['spatialAdapter']=sha(__file__)
    result['adapterBoundary']=sha(Path(__file__).with_name('adapters.py'))
    result['recordUtilities']=sha(Path(__file__).with_name('wetlab.py'))
    return result


def check(run, replay=False):
    run=Path(run);seal=read(run/'seal.json')
    require(seal['format']==FORMAT,'Spatial seal format')
    for name,digest in seal['files'].items():
        p=run/name
        require(p.resolve().is_relative_to(run.resolve()) and not p.is_symlink() and sha(p)==digest,'Sealed spatial artifact changed: '+name)
    reg=read(run/'registration.json')
    require(reg['adapterID']=='spatial-tissue/v1','Adapter identity changed')
    if replay: require(identity(reg['assay'])==reg['runtimeIdentity'],'Exact spatial runtime or adapter changed')
    return reg


def compare(run, observations):
    run=Path(run);reg=check(run);pred=read(run/'prediction.json');o=observations
    require(o['format']=='wet-lab-spatial-observations/v1','Observation schema')
    require(o['evidenceClass'] in ('measured','numerical-reference'),'Observation evidence class')
    require(o['units']=='mol/m3' and o['specimenID']==pred['specimenID'],'Observation specimen or units')
    require(o['compartmentIDs']==[c['id'] for c in pred['compartments']],'Observation spatial identities/order')
    require(o['geometrySHA256']==pred['geometrySHA256'],'Observation geometry binding')
    require(isinstance(o['provenance'],str) and o['provenance'].strip(),'Observation provenance missing')
    require(len(o['arms'])==2 and [a['name'] for a in o['arms']]==[a['name'] for a in pred['arms']],'Observation arm mismatch')
    metrics=[]
    for pa,oa in zip(pred['arms'],o['arms']):
        require(len(pa['samples'])==len(oa['samples']),'Incomplete observation times')
        for ps,os in zip(pa['samples'],oa['samples']):
            require(os['timeSeconds']==ps['timeSeconds'],'Observation time differs')
            actual=os['concentrationsMolPerM3'];expected=ps['concentrationsMolPerM3']
            require(len(actual)==len(expected) and all(type(v) in (int,float) and math.isfinite(v) and v>=0 for v in actual),'Invalid observation field')
            errors=[a-b for a,b in zip(expected,actual)]
            metrics.append({'arm':pa['name'],'timeSeconds':ps['timeSeconds'],
                            'rmse':math.sqrt(math.fsum(v*v for v in errors)/len(errors)),
                            'maximumAbsoluteError':max(map(abs,errors))})
    return {'format':FORMAT,'predictionSealSHA256':sha(run/'seal.json'),'observationEvidenceClass':o['evidenceClass'],
            'metrics':metrics,'units':'mol/m3','status':'compared',
            'biologicalValidation':'not-established','limits':LIMITS,
            'interpretation':'Numerical-reference comparison only.' if o['evidenceClass']=='numerical-reference' else
                             'Measured-data comparison retained; no automatic biological qualification.'}


class SpatialTissueAdapter(WetLabExperimentAdapter):
    adapter_id='spatial-tissue/v1'
    family='spatial-tissue'

    def catalog(self, config):
        a,path=load(config)
        return {'id':a['id'],'title':a['title'],'adapterID':self.adapter_id,'family':self.family,
                'limits':LIMITS,'intervention':'Declared spatial pulse field','sourceCitation':'','provenance':a['provenance'],
                'specimens':[{'id':s['id'],'cells':len(read(path(s['plan']))['cells']),
                              'extracellularCompartments':len(read(path(s['geometry']))['tetrahedra']),
                              'sampleTimesSeconds':read(path(s['plan']))['sampleTimesSeconds'],
                              'observationEvidenceClass':s.get('observationEvidenceClass','unavailable')} for s in a['specimens']]}

    def predict(self, config, runtime, workspace, selection):
        a,path=load(config);require(set(selection)=={'specimen'},'Spatial selection only supports a registered specimen/protocol')
        s=next((s for s in a['specimens'] if s['id']==selection['specimen']),None);require(s is not None,'Unsupported spatial specimen')
        runtime_id=identity(a)
        workspace=Path(workspace).resolve();workspace.mkdir(parents=True,exist_ok=True);run=workspace/uuid.uuid4().hex;run.mkdir()
        for key in ('geometry','mesh','material','plan'): shutil.copy2(path(s[key]),run/({'material':'material.nmatter'}.get(key,key+'.json')))
        for name in ('adapters.py','spatial_adapter.py','wetlab.py'):shutil.copy2(Path(__file__).with_name(name),run/name)
        reg={'format':FORMAT,'adapterID':self.adapter_id,'family':self.family,'createdAt':timestamp(),'assay':a,
             'specimen':s['id'],'specimenDefinition':s,'runtimeIdentity':runtime_id,
             'observationsPath':str(path(s['observations'])) if s.get('observations') else None,
             'analysis':'RMSE and maximum absolute mol/m3 error by arm/time across every compartment; no implicit qualification',
             'randomness':'none; deterministic numerical replay, not biological replication','limits':LIMITS,
             'transitions':[{'from':'Matter accepted FEM geometry','to':'extracellular transport graph','status':'numerical-geometric-transfer'},
                            {'from':'extracellular tracer','to':'intracellular tracer','status':'hypothetical-passive-exchange'},
                            {'from':'intracellular tracer','to':'phenotype','status':'not-implemented'}]}
        write(run/'registration.json',reg)
        try:
            # Reconstruct the accepted mechanics; do not trust a status string alone.
            native(a['runtime']['matterBinary'],[run/'mesh.json',run/'material.nmatter',run/'replayed-geometry.json'],run/'mechanics-log.json')
            require(read(run/'replayed-geometry.json')==read(run/'geometry.json'),'Accepted mechanics does not replay')
            native(a['runtime']['spatialBinary'],[run/'plan.json',run/'geometry.json',run/'prediction.json'],run/'prediction-log.json')
            write(run/'seal.json',{'format':FORMAT,'sealedAt':timestamp(),'files':inventory(run)})
        except Exception as error:
            write(run/'failure.json',{'message':str(error),'at':timestamp()});raise
        return run

    def reveal(self,run,runtime):
        run=Path(run);reg=check(run,True)
        require(not (run/'comparison.json').exists(),'Observations already revealed')
        require(reg['observationsPath'] is not None,'No spatial observations registered. Supply an independent measured or numerical-reference artifact in a new assay record.')
        source=Path(reg['observationsPath'])
        require(sha(source)==reg['specimenDefinition']['observationsSHA256'],'Registered observations changed')
        obs=read(source)
        require(obs['evidenceClass']==reg['specimenDefinition']['observationEvidenceClass'],'Observation evidence class differs')
        result=compare(run,obs)
        # Atomic pair publication; interrupted staging never becomes a comparison.
        with tempfile.TemporaryDirectory(dir=run,prefix='.reveal-') as tmp:
            tmp=Path(tmp);shutil.copy2(source,tmp/'observations.json');write(tmp/'comparison.json',result)
            (tmp/'observations.json').replace(run/'observations.json');(tmp/'comparison.json').replace(run/'comparison.json')
        return result

    def verify(self,run,runtime):
        run=Path(run);reg=check(run,True);a=reg['assay']
        with tempfile.TemporaryDirectory(prefix='wet-lab-replay-') as temp:
            temp=Path(temp)
            native(a['runtime']['matterBinary'],[run/'mesh.json',run/'material.nmatter',temp/'geometry.json'],temp/'mechanics-log.json')
            require(read(temp/'geometry.json')==read(run/'geometry.json'),'Mechanics replay differs')
            native(a['runtime']['spatialBinary'],[run/'plan.json',temp/'geometry.json',temp/'prediction.json'],temp/'prediction-log.json')
            require(read(temp/'prediction.json')==read(run/'prediction.json'),'Transport replay differs')
        if (run/'comparison.json').exists():
            require(sha(run/'observations.json')==reg['specimenDefinition']['observationsSHA256'],'Revealed observations changed')
            require(compare(run,read(run/'observations.json'))==read(run/'comparison.json'),'Comparison does not reconstruct')
        write(run/('replay-'+uuid.uuid4().hex+'.json'),{'verifiedAt':timestamp(),'status':'exact-native-replay','runtimeIdentity':reg['runtimeIdentity']})
        return {'status':'verified','adapterID':self.adapter_id,'run':run.name}

    def summary(self,run):
        run=Path(run);reg=check(run);result={'id':run.name,'adapterID':self.adapter_id,'family':self.family,
            'registration':reg,'recordDirectory':str(run.resolve()),'spatial':read(run/'prediction.json'),
            'revealed':(run/'comparison.json').exists(),'observationsAvailable':reg['observationsPath'] is not None}
        if result['revealed']:
            require(sha(run/'observations.json')==reg['specimenDefinition']['observationsSHA256'],'Revealed observations changed')
            result['comparison']=read(run/'comparison.json')
            require(result['comparison']==compare(run,read(run/'observations.json')),'Comparison changed')
        return result
