#!/usr/bin/env python3
"""Drive the real electronic CLI, artifact store and native numerical operations.

No API stand-ins. Requires the scoped production binary or the complete product.
Reference JSON is immutable output of the independent pinned PySCF exporter.
"""
import copy
import hashlib
import json
import pathlib
import subprocess
import sys


def run(binary, reference_directory, output_directory):
    binary=pathlib.Path(binary).resolve(strict=True)
    ref=pathlib.Path(reference_directory).resolve(strict=True)
    root=pathlib.Path(output_directory).resolve()
    root.mkdir(parents=True,exist_ok=True)
    store=root/'store'; checks=[]
    def check(flag,label):
        checks.append(dict(label=label,passed=bool(flag)))
        print(('PASS' if flag else 'FAIL'),label,flush=True)
        if not flag:
            raise AssertionError(label)
    def invoke(*args,expected=0):
        result=subprocess.run([str(binary),*map(str,args)],capture_output=True,text=True,timeout=600)
        with (root/'commands.log').open('a') as log:
            log.write(json.dumps(list(map(str,args)))+'\n'+result.stderr+'\n')
        if result.returncode!=expected:
            raise RuntimeError(f'Unexpected CLI exit {result.returncode}: {args}\n{result.stderr}\n{result.stdout}')
        return result.stdout
    def template(name):
        return json.loads(invoke('chemistry-template',name))
    def write(name,value):
        path=root/name
        path.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')
        return path
    def execute(label,request):
        source=write(label+'.request.json',request);output=root/(label+'.result.json')
        invoke('chemistry-run',source,'--store',store,'--output',output)
        return json.loads(output.read_text()),json.loads(pathlib.Path(str(output)+'.receipt.json').read_text())
    h2=json.loads((ref/'h2.json').read_text())
    budget=dict(maximumBytes=512*1024*1024,maximumBasisFunctions=256,maximumDeterminants=100000,maximumOperatorApplications=1_000_000_000)
    request=template('h2-tensor-ccsd')
    request.update(system=h2['system'],basis=h2['basis'],budget=budget)
    tensor,receipt=execute('tensor',request)
    check(abs(tensor['tensorCCSD']['result']['energyHartree']-h2['reference']['ccsd'])<1e-8,'advanced tensor CCSD executes through chemistry-run')
    _,repeat=execute('tensor-repeat',request)
    check(all(n['reused'] for n in repeat['nodes']),'unchanged calculation reuses every verified stage')
    direct=copy.deepcopy(request)
    reference=direct['calculation']['correlated']['reference']
    direct['calculation']['correlated']['solver']=template('solver-direct-fci')
    result,reused=execute('direct',direct)
    check(abs(result['directCI']['result']['roots'][0]['energyHartree']-h2['reference']['fciRoots'][0])<1e-8,'matrix-free CI uses common molecular workflow')
    check(all(n['reused'] for n in reused['nodes'] if n['identifier']!='result') and not next(n for n in reused['nodes'] if n['identifier']=='result')['reused'],
          'changing correlation solver preserves integral/RHF/Hamiltonian reuse')
    legacy=template('h2-ccsd');legacy.update(system=h2['system'],basis=h2['basis'],budget=budget)
    old,old_receipt=execute('legacy',legacy)
    check(abs(old['ccsd']['result']['energyHartree']-h2['reference']['ccsd'])<1e-8,'legacy electronic schema remains executable')
    check(all(n['reused'] for n in old_receipt['nodes'] if n['identifier']!='result'),'legacy and advanced routes share prepared artifacts')
    ri=copy.deepcopy(request)
    ri['calculation']={'densityFittedMP2':{'configuration':{'auxiliaryBasis':h2['auxiliaryBasis'],
        'fitting':{'metricEigenvalueThreshold':1e-10,'maximumAuxiliaryFunctions':4096,'pairBlockSize':256},'scf':reference}}}
    result,_=execute('ri-mp2',ri)
    check(abs(result['referenceEnergyHartree']+result['correlationEnergyHartree']-h2['reference']['riMP2'])<1e-8,
          'RI integrals, factorized RHF and RI-MP2 execute without dense AO-ERI workflow')
    smooth=template('h2-smooth-cpcm');smooth.update(system=h2['system'],basis=h2['basis'],budget=budget)
    result,_=execute('smooth-solvent',smooth)
    check(abs(result['scf']['energyHartree']-h2['reference']['smoothCPCM'])<1e-8,'smooth solvent workflow matches independent reference')
    ecc=template('h2-ecc-dmet');ecc.update(system=h2['system'],basis=h2['basis'],budget=budget)
    result,_=execute('ecc',ecc)
    check(result['converged'] and result['configuration']['mode']=='singleFragment' and result.get('matching') is None,
          'ECC molecular workflow preserves single-fragment no-matching semantics')
    check(abs(result['frame']['energyHartree']-h2['reference']['fciRoots'][0])<1e-8,'ECC physical energy survives workflow validation')
    _,ecc_repeat=execute('ecc-repeat',ecc)
    check(all(n['reused'] for n in ecc_repeat['nodes']),'ECC cache hit rebuilds its numerical acceptance checks')
    lih=json.loads((ref/'lih.json').read_text())
    sa=copy.deepcopy(request);sa.update(system=lih['system'],basis=lih['basis'])
    sa['calculation']['correlated']['solver']=template('solver-lih-sa-casscf')
    result,_=execute('state-average',sa)
    check(abs(result['multistateCASSCF']['result']['weightedEnergyHartree']-lih['reference']['saEnergy'])<5e-7,
          'state-averaged root-followed CASSCF executes through the shared molecular workflow')
    source=write('embedded.json',h2['embedded'])
    for solver,label in [('solver-tensor-ccsd','direct-tensor'),('solver-direct-fci','direct-spectrum'),('solver-ecc-dmet','direct-ecc')]:
        solver_path=write(label+'.solver.json',template(solver));out=root/(label+'.json')
        invoke('chemistry-solve',source,'--solver',solver_path,'--store',store,'--output',out)
        check(out.is_file() and pathlib.Path(str(out)+'.receipt.json').is_file(),label+' schema dispatched by chemistry-solve')
    bad=copy.deepcopy(request);bad['schema']='numivivo.org/unknown/v99'
    bad_path=write('invalid-schema.json',bad);bad_output=root/'invalid.result.json'
    invoke('chemistry-run',bad_path,'--store',store,'--output',bad_output,expected=1)
    check(not bad_output.exists(),'unknown request schema is rejected without a fallback result')
    bad_solver=write('invalid-solver.json',{'unknownMethod':{}})
    invoke('chemistry-solve',source,'--solver',bad_solver,'--store',store,'--output',bad_output,expected=1)
    check(not bad_output.exists(),'unknown solver is rejected without method substitution')
    source_path=root/'tensor.request.json';before=source_path.read_bytes()
    invoke('chemistry-run',source_path,'--store',store,'--output',source_path,expected=1)
    check(source_path.read_bytes()==before,'input/output aliasing is rejected before an overwrite')
    # Verify byte integrity at the real cache boundary, not in a mock store.
    digest=bytes(receipt['result']['bytes']).hex()
    object_path=store/'objects'/'sha256'/digest[:2]/digest[2:4]/digest
    original=object_path.read_bytes()
    object_path.write_bytes(original+b' ')
    try:
        invoke('chemistry-run',source_path,'--store',store,'--output',root/'corrupt.json',expected=1)
        check(not (root/'corrupt.json').exists(),'tampered cached object fails SHA-256 verification')
    finally:
        object_path.write_bytes(original)
    report=dict(schema='numivivo.org/electronic-cli-integration/v1',passed=all(c['passed'] for c in checks),checks=checks,
        binarySHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        scope='Real production electronic CLI, artifact store, planners, validators and numerical sources; scoped module unless supplied binary is the complete product')
    (root/'results.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    if len(sys.argv)!=4:
        raise SystemExit('usage: check_cli_integration.py binary reference-directory output-directory')
    run(*sys.argv[1:])
