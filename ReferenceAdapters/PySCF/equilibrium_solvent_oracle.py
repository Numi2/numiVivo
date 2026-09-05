#!/usr/bin/env python3
"""Pinned, independent equilibrium FCI/CASCI PCM references and real CLI checks.

Export uses only declared model inputs and PySCF. Check mode uses only the Python
standard library and a compiled native executable. No production Python runtime.
"""
import argparse
import hashlib
import json
from pathlib import Path

CASES = ('h2-full-eps4', 'lih-full-eps78', 'lih-cas22-eps78')


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')


def export(directory):
    import numpy as np
    import pyscf
    from pyscf import gto, scf, ao2mo, fci
    from pyscf.solvent import pcm
    from export_native_reference import basis, matrix
    if pyscf.__version__ != '2.8.0':
        raise RuntimeError('Expected PySCF 2.8.0')
    pyscf.lib.num_threads(1)
    records = []
    for name in CASES:
        is_h2 = name.startswith('h2')
        mol = gto.M(atom='H 0 0 -0.7; H 0 0 0.7' if is_h2 else 'Li 0 0 0; H 0 0 3.0',
                    unit='Bohr', basis='sto-3g', cart=True, verbose=0)
        mf = scf.RHF(mol).run(conv_tol=1e-13)
        if not mf.converged:
            raise RuntimeError('Unconverged orbital-frame reference')
        c = mf.mo_coeff
        n = mol.nao
        g = ao2mo.kernel(mol, c, compact=False).reshape(n,n,n,n)
        gas_h = c.T @ mf.get_hcore() @ c
        core = [0] if 'cas22' in name else []
        active = [1,2] if core else list(range(n))
        nelec = tuple(e-len(core) for e in mol.nelec)
        ga = g[np.ix_(active, active, active, active)]
        def solve(potential):
            h = gas_h + c.T @ potential @ c
            ha = h[np.ix_(active, active)].copy()
            for i,p in enumerate(active):
                for j,q in enumerate(active):
                    for a in core:
                        ha[i,j] += 2*g[p,q,a,a]-g[p,a,a,q]
            constant = mol.energy_nuc() + sum(2*h[a,a] for a in core)
            constant += sum(2*g[a,a,b,b]-g[a,b,b,a] for a in core for b in core)
            engine = fci.direct_spin1.FCI()
            engine.conv_tol = 1e-13
            energy, state = engine.kernel(ha, ga, len(active), nelec, ecore=constant)
            if not engine.converged:
                raise RuntimeError('Unconverged impurity oracle')
            p = np.zeros((n,n))
            for a in core:
                p[a,a] = 2
            p[np.ix_(active, active)] = engine.make_rdm1(state, len(active), nelec)
            return float(energy), p, c @ p @ c.T
        dielectric = 4.0 if is_h2 else 78.3
        radii = {1:1.2, 3:1.82}
        pcm_config = dict(dielectricConstant=dielectric, angularPoints=50, radiusScale=1.2,
                          radiiAngstrom={str(k):v for k,v in radii.items()}, extraCavitySpheres=[],
                          maximumTesserae=8192, tolerance=1e-11, maximumIterations=2000)
        solvent = pcm.PCM(mol)
        solvent.method = 'C-PCM'
        solvent.eps = dielectric
        solvent.lebedev_order = 11
        table = np.ones(119)
        for element,radius in radii.items():
            table[element] = radius
        solvent.radii_table = table*1.2/0.529177210544
        _, _, density = solve(np.zeros((n,n)))
        original = density.copy()
        previous = None
        for iteration in range(1,257):
            _, potential = solvent._get_vind(density)
            energy, p, actual = solve(potential)
            pol, actual_potential = solvent._get_vind(actual)
            gas = energy - np.einsum('ij,ji',actual,potential)
            total = gas + pol
            delta = np.linalg.norm(c.T @ mol.intor('int1e_ovlp') @ (actual-density) @ mol.intor('int1e_ovlp') @ c)
            if previous is not None and abs(total-previous)<1e-12 and delta<1e-11:
                break
            density = actual
            previous = total
        else:
            raise RuntimeError('Correlated PCM oracle did not converge')
        system = dict(nuclei=[dict(atomicNumber=int(z), positionBohr=r.tolist()) for z,r in zip(mol.atom_charges(), mol.atom_coords())],
                      pointCharges=[], alphaElectrons=mol.nelec[0], betaElectrons=mol.nelec[1])
        config = dict(maximumIterations=128, densityTolerance=1e-8, potentialToleranceHartree=1e-8,
                      energyToleranceHartree=1e-10, ciResidualTolerance=1e-11, damping=0.0)
        if core:
            config['partition'] = dict(doublyOccupiedCore=core, active=active, frozenOrbitals=[])
        request = dict(system=system, basis=basis(mol,'sto-3g'), coefficients=matrix(c), solvent=pcm_config,
                       configuration=config, budget=dict(maximumBasisFunctions=64, maximumBytes=268435456,
                       maximumDeterminants=512, maximumOperatorApplications=100000000))
        record = dict(identifier=name, request=request, energyHartree=float(total), gasEnergyHartree=float(gas),
                      polarizationEnergyHartree=float(pol), densityAO=matrix(actual), iterations=iteration,
                      densityChangeFromGas=float(np.linalg.norm(actual-original)))
        path = directory/(name+'.json')
        write(path,record)
        records.append(dict(file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    write(directory/'manifest.json',dict(schema='numivivo.org/equilibrium-solvent-oracle/v1',
          oracle='PySCF',version=pyscf.__version__,fixtures=records))


def check(binary, directory, out):
    import copy
    import math
    import subprocess
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; no existing results will be deleted')
    out.mkdir(parents=True,exist_ok=True)
    manifest = json.loads((directory/'manifest.json').read_text())
    if manifest['schema'] != 'numivivo.org/equilibrium-solvent-oracle/v1' or manifest['version'] != '2.8.0':
        raise RuntimeError('Unexpected oracle contract')
    if sorted(f['file'] for f in manifest['fixtures']) != sorted(n+'.json' for n in CASES):
        raise RuntimeError('Missing equilibrium fixture')
    checks, records = [], []
    def require(condition,label,value=None):
        item=dict(label=label,passed=bool(condition))
        if value is not None:item['value']=value
        checks.append(item)
        write(out/'checks.json',dict(checks=checks,passed=all(x['passed'] for x in checks)))
        if not condition:raise RuntimeError(label)
    def execute(request,name,fail=False):
        source=out/(name+'.request.json');result=out/(name+'.result.json')
        write(source,dict(schema='numivivo.org/reaction-calculation/v1',calculation=dict(correlatedSolvent=dict(request=request))))
        process=subprocess.run([str(binary),'reaction-run',str(source),'--store',str(out/'store'),'--output',str(result)],
                               capture_output=True,text=True,timeout=600)
        (out/(name+'.log')).write_text(process.stdout+'\n'+process.stderr)
        if fail:
            require(process.returncode!=0 and not result.exists() and not Path(str(result)+'.receipt.json').exists(),name)
            return None
        if process.returncode:raise RuntimeError(process.stderr)
        return json.loads(result.read_text())['correlatedSolvent']['result'],json.loads(Path(str(result)+'.receipt.json').read_text())
    for item in manifest['fixtures']:
        path=directory/item['file'];data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()!=item['sha256']:raise RuntimeError('Corrupted oracle')
        record=json.loads(data);name=record['identifier'];records.append(record)
        actual, receipt=execute(record['request'],name)
        for field in ('energyHartree','gasEnergyHartree'):
            error=abs(actual[field]-record[field])
            require(math.isfinite(error) and error<2e-8,name+' '+field,error)
        error=abs(actual['equilibriumField']['polarizationEnergyHartree']-record['polarizationEnergyHartree'])
        require(math.isfinite(error) and error<2e-8,name+' polarization',error)
        a,b=actual['totalDensityAO'],record['densityAO']
        require(a['rows']==b['rows'] and a['columns']==b['columns'] and len(a['values'])==len(b['values']),name+' density shape')
        error=math.sqrt(sum((x-y)**2 for x,y in zip(a['values'],b['values'])))
        require(math.isfinite(error) and error<2e-7,name+' density agreement',error)
        require(actual['densityResidual']<=record['request']['configuration']['densityTolerance'],name+' density closure')
        require(not receipt['reused'],name+' fresh native calculation')
        again,second=execute(record['request'],name+'-cached')
        require(second['reused'] and again==actual and second['result']==receipt['result'],name+' verified cache hit')
        if not name.startswith('h2'):
            require(record['densityChangeFromGas']>1e-4,name+' nontrivial correlated density polarization',record['densityChangeFromGas'])
    request=copy.deepcopy(records[-1]['request']);request['configuration']['maximumIterations']=2
    execute(request,'insufficient-equilibration',fail=True)
    request=copy.deepcopy(records[0]['request']);request['coefficients']['values'][0]+=0.1
    execute(request,'nonorthogonal-frame',fail=True)
    request=copy.deepcopy(records[0]['request']);request['configuration']['damping']=1
    execute(request,'invalid-density-damping',fail=True)
    request=copy.deepcopy(records[-1]['request']);request['configuration']['partition']['active']=[0,1]
    execute(request,'overlapping-core-active',fail=True)
    # Actual nuclear displacements use the equilibrium method, not the legacy
    # frozen-reference approximation. Acceptance reconstructs gradients/modes.
    source=out/'equilibrium-minimum.request.json';result=out/'equilibrium-minimum.result.json'
    p=subprocess.run([str(binary),'reaction-template','h2-equilibrium-minimum','--output',str(source)],capture_output=True,text=True,timeout=60)
    if p.returncode:raise RuntimeError(p.stderr)
    for attempt in range(2):
        dest=result if attempt==0 else out/'equilibrium-minimum-cached.result.json'
        p=subprocess.run([str(binary),'reaction-run',str(source),'--store',str(out/'store'),'--output',str(dest)],capture_output=True,text=True,timeout=600)
        (out/f'equilibrium-minimum-{attempt}.log').write_text(p.stdout+'\n'+p.stderr)
        if p.returncode:raise RuntimeError(p.stderr)
        point=json.loads(dest.read_text())['qualified']['point']
        receipt=json.loads(Path(str(dest)+'.receipt.json').read_text())
        require(receipt['reused']==bool(attempt),'equilibrium nuclear cache state '+str(attempt))
        require(point['request']['model']['solver']=='equilibriumFullCI','explicit correlated nuclear solver '+str(attempt))
        modes=point['thermochemistry']['modes']
        require(len(modes['signedFrequenciesCM'])==1 and modes['signedFrequenciesCM'][0]>0
                and modes['maximumGradient']<=point['request']['thermochemistry']['maximumGradientTolerance'],
                'equilibrium nuclear minimum and full Hessian acceptance '+str(attempt))
        if attempt==0:first=point
        else:require(point==first,'identical reconstructed equilibrium nuclear result on reuse')
    write(out/'checks.json',dict(schema='numivivo.org/equilibrium-solvent-cli-checks/v1',checks=checks,passed=True,
          executableSHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),
          scope='fixed-orbital FCI/CASCI equilibrium PCM, full-CI equilibrium nuclear minimum and actual CLI/cache; no overlapping-fragment density, paper reproduction or rate certificate'))
    print('PASS',len(checks),'equilibrium solvent comparisons and CLI checks')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('export');a.add_argument('directory',type=Path)
    a=sub.add_parser('check');a.add_argument('binary',type=Path);a.add_argument('directory',type=Path);a.add_argument('out',type=Path)
    args=p.parse_args()
    if args.command=='export':export(args.directory)
    else:check(args.binary.resolve(),args.directory,args.out)
