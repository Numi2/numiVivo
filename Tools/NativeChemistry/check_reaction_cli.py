#!/usr/bin/env python3
"""Run real nuclear/solvent/thermochemical CLI operations and cache rejection checks.
No numerical or storage substitutes are constructed. Results are finite-basis
local-saddle and ideal RRHO estimates, not biological rate measurements.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def run(binary, out):
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; existing artifacts are not removed')
    out.mkdir(parents=True, exist_ok=True)
    checks, counter = [], 0
    def check(condition, label):
        checks.append(dict(label=label, passed=bool(condition)))
        write(out/'checks.json', dict(checks=checks,passed=all(c['passed'] for c in checks)))
        if not condition: raise RuntimeError(label)
        print('PASS', label, flush=True)
    def execute(*args, failure=False):
        nonlocal counter
        counter += 1
        p = subprocess.run([str(binary), *map(str,args)], capture_output=True, text=True, timeout=1200)
        (out/f'command-{counter:02d}.log').write_text('args: '+json.dumps(list(map(str,args)))+'\n'+p.stdout+'\n'+p.stderr)
        if (p.returncode == 0) == failure:
            raise RuntimeError(f'Unexpected exit {p.returncode}: {args}\n{p.stderr}')
        return p
    store=out/'artifacts'
    def calculate(request,name):
        source=out/f'{name}.request.json'; destination=out/f'{name}.result.json'
        write(source,request)
        execute('reaction-run',source,'--store',store,'--output',destination)
        return json.loads(destination.read_text()),json.loads(Path(str(destination)+'.receipt.json').read_text())
    points, requests = {}, {}
    for template in ['h2-minimum','h3-saddle','h-atom']:
        path=out/f'{template}.template.json'
        execute('reaction-template',template,'--output',path)
        request=json.loads(path.read_text()); requests[template]=request
        result,receipt=calculate(request,template)
        point=result['qualified']['point'];points[template]=point
        check(not receipt['reused'] and point['electronicEvaluations']>0,f'{template}: native calculation executes through the production CLI and artifact store')
    h2,h3,h=points['h2-minimum'],points['h3-saddle'],points['h-atom']
    repeated,receipt=calculate(requests['h3-saddle'],'h3-repeated')
    check(receipt['reused'] and repeated['qualified']['point']==h3,'qualified saddle cache hit reconstructs acceptance and preserves output identity')
    check(len(h3['thermochemistry']['modes']['signedFrequenciesCM'])==4
          and sum(f<0 for f in h3['thermochemistry']['modes']['signedFrequenciesCM'])==1,
          'H3 output retains every stable and unstable vibrational mode')
    barrier_request=dict(schema='numivivo.org/reaction-calculation/v1',calculation={'harmonicBarrier':{'saddle':h3,'reactants':[h2,h]}})
    barrier,_=calculate(barrier_request,'barrier')
    b=barrier['harmonicBarrier']['result']
    check(b['reactantMolecularity']==2 and b['activationGibbsHartree']>0,
          'balanced separate-reactant harmonic Gibbs estimate executes without a rate export')
    check('connectivity' in b['model'] and 'transmission' in b['model'],
          'barrier output retains its unqualified connectivity and transmission assumptions')
    descent_request=dict(schema='numivivo.org/reaction-calculation/v1',calculation={'descent':{'saddle':h3,'configuration':{
        'initialDisplacementMassWeighted':0.03,'stepMassWeighted':0.08,'maximumStepsPerDirection':12,'endpointMaximumGradient':1e-4}}})
    branch,_=calculate(descent_request,'descent')
    d=branch['descent']['result']
    check(len(d['forward'])>5 and len(d['reverse'])>5 and 'no step-size-converged IRC' in d['interpretation'],
          'two native descent branches are exported with explicit endpoint and discretization limits')
    solvated_template=out/'solvated.template.json'
    execute('reaction-template','h2-solvated-path','--output',solvated_template)
    solvent_request=json.loads(solvated_template.read_text())
    solvated,solvent_receipt=calculate(solvent_request,'solvated')
    s=solvated['solvatedPath']['result']
    check(s['path']['converged'] and len(s['referenceSolventFields'])==5 and len(s['frozenFieldConstantsHartree'])==5,
          'smooth reference polarization and shared ECC execute through the common reaction workflow')
    check('not correlated self-consistent PCM' in s['convention'],
          'reference-frozen correlated solvent approximation is preserved in the result')
    repeated,second=calculate(solvent_request,'solvated-repeated')
    check(second['reused'] and second['result']==solvent_receipt['result'] and repeated==solvated,
          'solvated shared-path cache reuse rechecks field and electronic-path acceptance')
    def rejected(request,name):
        source=out/f'{name}.request.json';dest=out/f'{name}.rejected.json';write(source,request)
        execute('reaction-run',source,'--store',store,'--output',dest,failure=True)
        return not dest.exists() and not Path(str(dest)+'.receipt.json').exists()
    bad=copy.deepcopy(requests['h3-saddle']);bad['schema']='unsupported/v99'
    check(rejected(bad,'schema'),'unknown reaction schema cannot produce a fallback calculation')
    bad=copy.deepcopy(requests['h2-minimum']);bad['calculation']['qualify']['request']['differences']['maximumEnergyEvaluations']=1
    check(rejected(bad,'budget'),'exhausted electronic work cannot publish a qualified stationary result')
    bad=copy.deepcopy(requests['h2-minimum']);bad['calculation']['qualify']['request']['operation']='characterize'
    check(rejected(bad,'nonstationary'),'nonstationary nuclear coordinates cannot publish thermochemical qualification')
    bad=copy.deepcopy(barrier_request);bad['calculation']['harmonicBarrier']['saddle']['thermochemistry']['gibbsEnergyHartree']=99.0
    check(rejected(bad,'forged-gibbs'),'forged saddle thermochemistry is numerically rejected before barrier export')
    bad=copy.deepcopy(barrier_request);bad['calculation']['harmonicBarrier']['reactants']=[h2,h2]
    check(rejected(bad,'composition'),'unbalanced reaction composition cannot produce an activation estimate')
    bad=copy.deepcopy(requests['h2-minimum']);bad['calculation']={'unknown':{}}
    check(rejected(bad,'method'),'unknown nuclear method is not substituted with a different solver')
    source=out/'h3-saddle.request.json';original=source.read_bytes()
    execute('reaction-run',source,'--store',store,'--output',source,failure=True)
    check(source.read_bytes()==original,'output cannot overwrite its own nuclear request')
    alias=out/'hardlink.json';os.link(source,alias)
    execute('reaction-run',source,'--store',store,'--output',alias,failure=True)
    check(source.read_bytes()==original and os.path.samefile(source,alias),'hard-link output alias is rejected without overwriting input')
    inside=store/'result.json';execute('reaction-run',source,'--store',store,'--output',inside,failure=True)
    check(not inside.exists(),'output cannot overwrite an artifact-store file')
    digest=bytes(solvent_receipt['result']['bytes']).hex()
    object_path=store/'objects'/'sha256'/digest[:2]/digest[2:4]/digest
    data=object_path.read_bytes()
    if hashlib.sha256(data).hexdigest()!=digest:raise RuntimeError('Unrecognized artifact layout or digest')
    try:
        object_path.write_bytes(data+b'\n')
        check(rejected(solvent_request,'corrupt-cache'),'tampered cached solvent/ECC result fails digest verification')
    finally:object_path.write_bytes(data)
    write(out/'checks.json',dict(schema='numivivo.org/reaction-cli-checks/v1',passed=True,checks=checks,
          executableSHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),
          scope='real production nuclear qualification/solvent path classes and artifact store, scoped Apple build; no paper reaction or rates'))
    print('PASS',len(checks),'reaction CLI/store checks',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('binary',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();run(a.binary.resolve(),a.output.resolve())
