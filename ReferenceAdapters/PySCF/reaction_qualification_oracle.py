#!/usr/bin/env python3
"""Independent PySCF references for local saddle/RRHO and frozen-field ECC.
No native result is read during export. Primitive constants are PySCF 2.8.0's;
small native/PySCF basis rounding and CODATA differences are not fitted away.
"""
import argparse
import hashlib
import json
from pathlib import Path


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, sort_keys=True, indent=2, allow_nan=False) + '\n')


def export(out):
    import numpy as np
    import pyscf
    from scipy.linalg import eigh
    from scipy.optimize import minimize, minimize_scalar
    from pyscf import ao2mo, fci, gto, scf
    from pyscf.hessian import thermo
    from pyscf.solvent import pcm
    from types import SimpleNamespace
    if pyscf.__version__ != '2.8.0':
        raise RuntimeError('Expected PySCF 2.8.0')
    pyscf.lib.num_threads(1)
    def molecule(positions, spin=0, elements=None):
        return gto.M(atom=list(zip(elements or ['H']*len(positions), positions)), spin=spin,
                     basis='sto-3g', unit='Bohr', verbose=0)
    def ci(m, core=None, constant=None):
        d, c = eigh(m.intor('int1e_ovlp'))
        x = c / np.sqrt(d)
        h = m.intor('int1e_kin') + m.intor('int1e_nuc') if core is None else core
        g = ao2mo.restore(1, ao2mo.kernel(m, x), m.nao)
        solver = fci.direct_spin1.FCI(); solver.conv_tol = 1e-13
        e, _ = solver.kernel(x.T @ h @ x, g, m.nao, m.nelec)
        if not solver.converged: raise RuntimeError('FCI oracle did not converge')
        return e + (m.energy_nuc() if constant is None else constant)
    def numerical_hessian(positions, spin):
        x = np.asarray(positions).reshape(-1)
        def energy(x): return ci(molecule(x.reshape(-1, 3), spin))
        center = energy(x)
        def hessian(step):
            h = np.zeros((len(x), len(x))); unit = np.eye(len(x))*step
            for i in range(len(x)):
                h[i, i] = (energy(x+unit[i])-2*center+energy(x-unit[i]))/step**2
                for j in range(i):
                    h[i, j] = h[j, i] = (energy(x+unit[i]+unit[j])-energy(x+unit[i]-unit[j])
                                        -energy(x-unit[i]+unit[j])+energy(x-unit[i]-unit[j]))/(4*step**2)
            return h
        return center, (4*hessian(0.0005)-hessian(0.001))/3
    def summary(m, energy, hessian):
        n = m.natm
        vibrations = thermo.harmonic_analysis(m, hessian.reshape(n,3,n,3).transpose(0,2,1,3))
        frequencies = [float(-abs(z.imag) if abs(z.imag)>1e-8 else z.real) for z in vibrations['freq_wavenumber']]
        thermal = thermo.thermo(SimpleNamespace(mol=m, e_tot=energy), vibrations['freq_au'], temperature=298.15, pressure=101325)
        return dict(positionsBohr=m.atom_coords().tolist(), massesDa=m.atom_mass_list(isotope_avg=True).tolist(),
                    electronicEnergy=float(energy), frequenciesCM=frequencies,
                    zeroPoint=float(thermal['ZPE'][0]), enthalpyCorrection=float(thermal['H_tot'][0]-energy),
                    entropy=float(thermal['S_tot'][0]), gibbs=float(thermal['G_tot'][0]),
                    hessian=dict(rows=3*n,columns=3*n,values=hessian.reshape(-1).tolist()))
    records={}
    for name in ['h2','h3']:
        spin = 0 if name=='h2' else 1
        def positions(distance):
            return [[0,0,-distance/2],[0,0,distance/2]] if name=='h2' else [[0,0,-distance],[0,0,0],[0,0,distance]]
        solved=minimize_scalar(lambda r:ci(molecule(positions(r),spin)),bounds=(1.1,2.3),method='bounded',options={'xatol':1e-10})
        if not solved.success:raise RuntimeError('Oracle geometry optimization failed')
        x=positions(solved.x);energy,hessian=numerical_hessian(x,spin)
        records[name]=summary(molecule(x,spin),energy,hessian)
        records[name]['distanceBohr']=float(solved.x)
    atom=molecule([[0,0,0]],1)
    records['h']=summary(atom,ci(atom),np.zeros((3,3)))
    # Independent analytic HF water Hessian validates the nonlinear rotor/mode
    # treatment, not a claim that this fixture executes native water FCI.
    def water(parameters):
        r,angle=parameters
        return molecule([[0,0,0],[r*np.sin(angle/2),0,r*np.cos(angle/2)],[-r*np.sin(angle/2),0,r*np.cos(angle/2)]],elements=['O','H','H'])
    def objective(parameters):
        r,t=parameters;m=water(parameters);mf=scf.RHF(m).run(conv_tol=1e-13)
        if not mf.converged:raise RuntimeError('Water HF oracle did not converge')
        g=mf.nuc_grad_method().kernel()
        dr=np.array([[np.sin(t/2),0,np.cos(t/2)],[-np.sin(t/2),0,np.cos(t/2)]])
        dt=np.array([[r/2*np.cos(t/2),0,-r/2*np.sin(t/2)],[-r/2*np.cos(t/2),0,-r/2*np.sin(t/2)]])
        return mf.e_tot,np.array([(g[1:]*dr).sum(),(g[1:]*dt).sum()])
    solved=minimize(objective,[1.8,1.9],jac=True,method='BFGS',options={'gtol':1e-8})
    if np.max(np.abs(objective(solved.x)[1]))>1e-7:raise RuntimeError('Water oracle is not stationary')
    m=water(solved.x);mf=scf.RHF(m).run(conv_tol=1e-13)
    h=mf.Hessian().kernel().transpose(0,2,1,3).reshape(9,9)
    records['water']=summary(m,mf.e_tot,h)
    records['water']['gradient']=mf.nuc_grad_method().kernel().tolist()
    # Explicitly matched frozen RHF-reference Gaussian surface fields. The
    # correlated FCI density is not used to re-equilibrate these fields.
    path=[]
    for distance in [1.2,1.4,1.6,1.8,2.0]:
        m=molecule([[0,0,-distance/2],[0,0,distance/2]])
        solvent=pcm.PCM(m);solvent.method='C-PCM';solvent.eps=4;solvent.lebedev_order=11
        radii=np.ones(119);radii[1]=1.2;solvent.radii_table=radii*1.2/0.529177210544
        mf=scf.RHF(m).PCM(solvent).run(conv_tol=1e-12)
        if not mf.converged:raise RuntimeError('Smooth reference SCF failed')
        density=mf.make_rdm1();polarization,potential=solvent._get_vind(density)
        correction=float(polarization-np.einsum('ij,ji',density,potential))
        core=m.intor('int1e_kin')+m.intor('int1e_nuc')+potential
        correlated=ci(m,core=core,constant=m.energy_nuc()+correction)
        path.append(dict(distanceBohr=distance,energy=float(correlated),rhf=float(mf.e_tot),
                         polarization=float(polarization),frozenFieldConstant=correction))
    records['solvatedPath']=path
    write(out/'reference.json',dict(schema='numivivo.org/reaction-oracle/v1',oracle='PySCF',version=pyscf.__version__,records=records))
    write(out/'manifest.json',dict(sha256=hashlib.sha256((out/'reference.json').read_bytes()).hexdigest(),required=['h','h2','h3','water','solvatedPath']))
    print('Generated independent nuclear, harmonic and reference-polarized path fixtures')


def compare(native, reference, out):
    import math
    data=(reference/'reference.json').read_bytes();manifest=json.loads((reference/'manifest.json').read_text())
    if hashlib.sha256(data).hexdigest()!=manifest['sha256']:raise RuntimeError('Reference hash mismatch')
    oracle=json.loads(data)
    if oracle['version']!='2.8.0' or oracle['schema']!='numivivo.org/reaction-oracle/v1':raise RuntimeError('Unknown reference protocol')
    expected=oracle['records'];checks=[]
    if set(manifest['required'])!=set(expected):raise RuntimeError('Incomplete reaction references')
    def check(label,a,b,tolerance):
        error=abs(a-b)
        checks.append(dict(label=label,actual=a,expected=b,error=error,tolerance=tolerance,passed=math.isfinite(a) and math.isfinite(b) and error<=tolerance))
    for name in ['h','h2','h3']:
        actual=json.loads((native/(name+'.json')).read_text())['thermochemistry'];r=expected[name]
        check(name+' electronic',actual['electronicEnergyHartree'],r['electronicEnergy'],1e-8)
        for field,key,tol in [('zeroPointEnergyHartree','zeroPoint',1e-7),('enthalpyCorrectionHartree','enthalpyCorrection',1e-7),
                              ('gasEntropyHartreePerK','entropy',1e-10),('gibbsEnergyHartree','gibbs',1e-7)]:
            check(name+' '+key,actual[field],r[key],tol)
        frequencies=actual['modes']['signedFrequenciesCM']
        if len(frequencies)!=len(r['frequenciesCM']):raise RuntimeError('Missing modes')
        for i,(a,b) in enumerate(zip(frequencies,r['frequenciesCM'])):check(name+f' frequency {i}',a,b,0.05)
    water=json.loads((native/'water.json').read_text())
    for field,key,tol in [('zeroPointEnergyHartree','zeroPoint',1e-7),('enthalpyCorrectionHartree','enthalpyCorrection',1e-7),('gasEntropyHartreePerK','entropy',1e-10),('gibbsEnergyHartree','gibbs',1e-7)]:
        check('water '+key,water[field],expected['water'][key],tol)
    barrier=json.loads((native/'barrier.json').read_text())
    check('electronic separated-reactant barrier',barrier['electronicBarrierHartree'],expected['h3']['electronicEnergy']-expected['h2']['electronicEnergy']-expected['h']['electronicEnergy'],1e-8)
    check('RRHO separated-reactant Gibbs estimate',barrier['activationGibbsHartree'],expected['h3']['gibbs']-expected['h2']['gibbs']-expected['h']['gibbs'],3e-7)
    path=json.loads((native/'solvated-path.json').read_text())
    if len(path['path']['pointResults'])!=5:raise RuntimeError('Missing solvated path geometries')
    for i,ref in enumerate(expected['solvatedPath']):
        check(f'solvated ECC energy {i}',path['path']['pointResults'][i]['frame']['energyHartree'],ref['energy'],1e-8)
        check(f'field constant {i}',path['frozenFieldConstantsHartree'][i],ref['frozenFieldConstant'],1e-8)
    passed=all(c['passed'] for c in checks)
    write(out,dict(schema='numivivo.org/reaction-comparison/v1',checks=checks,passed=passed,referenceSHA256=manifest['sha256'],
                   scope='finite STO-3G local molecular saddle/RRHO and HF-reference-polarized full-bath ECC; no paper or rate reproduction'))
    if not passed:raise RuntimeError('Reaction conformance failed; see comparison report')
    print('PASS',len(checks),'independent reaction-qualification comparisons')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    a=sub.add_parser('export');a.add_argument('out',type=Path)
    a=sub.add_parser('compare');a.add_argument('native',type=Path);a.add_argument('reference',type=Path);a.add_argument('out',type=Path)
    args=parser.parse_args()
    if args.command=='export':export(args.out)
    else:compare(args.native,args.reference,args.out)
