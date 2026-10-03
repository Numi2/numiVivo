#!/usr/bin/env python3
"""Bounded native spatial acceptance, independently checked on the discrete graph.

Requires NumPy/SciPy only for verification, never for native prediction.
Every case and threshold is registered before execution; failures are retained.
"""
import argparse, copy, json, shutil, subprocess, traceback
from pathlib import Path
import numpy as np
from wetlab import read, write, sha, require
from spatial_reference import reference
from spatial_adapter import SpatialTissueAdapter, compare, check
from adapters import adapter_for_run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--assay',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rna-record',type=Path,required=True);p.add_argument('--rna-binary',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    acceptance={'purpose':'numerical instrument qualification, not a biological study',
        'cases':['baseline','half-step','no diffusion','no membrane','invalid geometry','exact replay','reference reveal','changed artifact rejection','RNA backward compatibility'],
        'stoppingRule':'one execution per declared case; retain failures',
        'reference':'independent FP64 matrix exponential on the admitted graph',
        'concentrationTolerance':2e-4,'massRelativeTolerance':2e-5}
    write(a.output/'acceptance-plan.json',acceptance)
    config=read(a.assay);specimen=config['specimens'][0];root=a.assay.parent
    native=config['runtime']['spatialBinary'];plan=read(root/specimen['plan']);geometry=root/specimen['geometry']
    results={};adapter=SpatialTissueAdapter();run=None
    def record(name,fn):
        try:results[name]={'passed':True,**(fn() or {})}
        except Exception as e:results[name]={'passed':False,'error':str(e),'traceback':traceback.format_exc()}
        write(a.output/(str(len(results)).zfill(2)+'-case.json'),{'case':name,**results[name]})
    def fail(fn):
        try:fn()
        except (ValueError,KeyError):return
        raise AssertionError('invalid record/input was accepted')
    def case(name,mutate=lambda p:None):
        q=copy.deepcopy(plan);mutate(q);folder=a.output/name;folder.mkdir();write(folder/'plan.json',q)
        proc=subprocess.run([native,str(folder/'plan.json'),str(geometry),str(folder/'result.json')],capture_output=True,text=True)
        write(folder/'native-log.json',{'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr})
        require(proc.returncode==0,'native case failed')
        result=read(folder/'result.json');oracle=reference(folder/'plan.json',geometry);write(folder/'reference.json',oracle)
        errors=[]
        for actual,expected in zip(result['arms'],oracle['arms']):
            require(actual['maximumRelativeMassError']<=acceptance['massRelativeTolerance'],'amount balance')
            for s,t in zip(actual['samples'],expected['samples']):
                error=float(np.max(np.abs(np.array(s['concentrationsMolPerM3'])-t['concentrationsMolPerM3'])));errors.append(error)
                require(error<=acceptance['concentrationTolerance'],'independent concentration tolerance')
        if name=='no-diffusion':
            for sample in result['arms'][1]['samples']:
                require(all(sample['concentrationsMolPerM3'][i]==0 for i in range(len(q['cells'])) if i not in (0,1)),'disabled diffusion leaked')
        if name=='no-membrane':
            require(all(all(v==0 for v in s['concentrationsMolPerM3'][len(q['extracellularInitialMolPerM3']):]) for s in result['arms'][1]['samples']),'disabled membrane leaked')
        return {'maximumAbsoluteErrorMolPerM3':max(errors),'maximumRelativeMassError':max(x['maximumRelativeMassError'] for x in result['arms']),
                'predictionSHA256':sha(folder/'result.json')}
    record('baseline',lambda:case('baseline'))
    record('half-step',lambda:case('half-step',lambda q:q.update(timeStepSeconds=q['timeStepSeconds']/2)))
    record('no diffusion',lambda:case('no-diffusion',lambda q:q.update(diffusionSquareMetresPerSecond=0)))
    def no_membrane(q):
        for c in q['cells']:c['membranePermeabilityMetresPerSecond']=0
    record('no membrane',lambda:case('no-membrane',no_membrane))
    def invalid():
        g=read(geometry);g['tetrahedra'][0][0],g['tetrahedra'][0][1]=g['tetrahedra'][0][1],g['tetrahedra'][0][0];write(a.output/'inverted.json',g)
        proc=subprocess.run([native,str(root/specimen['plan']),str(a.output/'inverted.json'),str(a.output/'invalid-result.json')],capture_output=True,text=True)
        write(a.output/'invalid-log.json',{'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr})
        require(proc.returncode!=0 and not (a.output/'invalid-result.json').exists(),'invalid geometry published')
        return {'rejection':proc.stderr.strip()}
    record('invalid geometry',invalid)
    def exact():
        nonlocal run
        run=adapter.predict(a.assay,{},a.output/'experiments',{'specimen':specimen['id']})
        require(not (run/'observations.json').exists(),'observations exposed before seal')
        return {'run':run.name,**adapter.verify(run,{})}
    record('exact replay',exact)
    def reveal():
        result=adapter.reveal(run,{})
        require(result['observationEvidenceClass']=='numerical-reference' and result['biologicalValidation']=='not-established','numerical evidence promoted')
        adapter.verify(run,{})
        return {'maximumAbsoluteErrorMolPerM3':max(x['maximumAbsoluteError'] for x in result['metrics']),'biologicalValidation':result['biologicalValidation']}
    record('reference reveal',reveal)
    def tamper():
        rejected=[]
        for file in ('prediction.json','geometry.json','plan.json','observations.json','comparison.json'):
            folder=a.output/('tamper-'+file.split('.')[0]);shutil.copytree(run,folder)
            with (folder/file).open('a') as f:f.write(' ')
            if file=='comparison.json':
                o=read(folder/file);o['biologicalValidation']='validated';(folder/file).write_text(json.dumps(o))
            fail(lambda:adapter.summary(folder));rejected.append(file)
        for key,value in [('units','mmol/L'),('evidenceClass','validated'),('geometrySHA256','wrong')]:
            o=read(run/'observations.json');o[key]=value;fail(lambda:compare(run,o));rejected.append(key)
        o=read(run/'observations.json');o['arms'][0]['samples'][0]['timeSeconds']=.01;fail(lambda:compare(run,o));rejected.append('time')
        return {'rejected':rejected}
    record('changed artifact rejection',tamper)
    record('RNA backward compatibility',lambda:adapter_for_run(a.rna_record).verify(a.rna_record,{'binary':a.rna_binary}))
    write(a.output/'checks.json',{'acceptance':acceptance,'results':results,'allPassed':len(results)==len(acceptance['cases']) and all(r['passed'] for r in results.values())})
    if not all(r['passed'] for r in results.values()):raise SystemExit('Some acceptance cases failed; inspect checks.json')
    print(a.output/'checks.json')
if __name__=='__main__':main()
