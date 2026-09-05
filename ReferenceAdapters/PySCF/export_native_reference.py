#!/usr/bin/env python3
"""Independent pinned PySCF data. Never imported by the native runtime.

Every reference solver must converge before its output is published. Cartesian
AO functions are individually normalized to match NumiVivo's explicit basis
contract. RI uses native PySCF two/three-centre integrals and a declared metric
cutoff; the reference is not calculated from any NumiVivo result.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import argparse
import hashlib
import json
import pathlib
import platform
import sys
import numpy as np
import pyscf
from pyscf import gto, scf, cc, mp, fci, df, ao2mo, mcscf
from pyscf.solvent import pcm


def matrix(a):
    a = np.asarray(a)
    return dict(rows=a.shape[0], columns=a.shape[1], values=a.ravel().tolist())


def basis(molecule, name):
    shells = []
    for atom in range(molecule.natm):
        for shell in molecule._basis[molecule.atom_symbol(atom)]:
            if not isinstance(shell[1], (list, tuple)):
                raise ValueError('Spinor/kappa bases are outside Cartesian conformance')
            for k in range(len(shell[1]) - 1):
                shells.append(dict(nucleusIndex=atom, angularMomentum=shell[0],
                    primitives=[dict(exponent=row[0], coefficient=row[k+1]) for row in shell[1:]]))
    return dict(identifier=name, representation='normalized-cartesian', shells=shells,
                source='PySCF 2.8.0 explicit basis coefficients; reference export')


def generate(name, atom):
    molecule = gto.M(atom=atom, basis='sto-3g', unit='Bohr', cart=True, verbose=0)
    overlap = molecule.intor('int1e_ovlp')
    core = molecule.intor('int1e_kin') + molecule.intor('int1e_nuc')
    eri = molecule.intor('int2e')
    scale = 1/np.sqrt(np.diag(overlap))
    auxiliary = df.addons.make_auxmol(molecule, 'weigend')
    auxiliary.cart = True
    auxiliary_scale = 1/np.sqrt(np.diag(auxiliary.intor('int1e_ovlp')))
    metric = auxiliary.intor('int2c2e')*auxiliary_scale[:, None]*auxiliary_scale[None, :]
    three = np.einsum('pqL,p,q,L->pqL', df.incore.aux_e2(molecule, auxiliary, aosym='s1'),
                      scale, scale, auxiliary_scale)
    eigenvalues, vectors = np.linalg.eigh(metric)
    retained = eigenvalues > 1e-10
    factors = np.einsum('pqL,LQ->pqQ', three, vectors[:, retained]/np.sqrt(eigenvalues[retained]))
    fitted = np.einsum('pqL,rsL->pqrs', factors, factors)
    hf = scf.RHF(molecule).run(conv_tol=1e-13)
    if not hf.converged:
        raise RuntimeError(name + ': reference HF did not converge')
    coupled = cc.CCSD(hf)
    coupled.conv_tol = 1e-12
    coupled.conv_tol_normt = 1e-10
    coupled.kernel()
    if not coupled.converged:
        raise RuntimeError(name + ': reference CCSD did not converge')
    rif = scf.RHF(molecule)
    physical_fitted = np.einsum('pqrs,p,q,r,s->pqrs', fitted, 1/scale, 1/scale, 1/scale, 1/scale)
    rif._eri = ao2mo.restore(8, physical_fitted, molecule.nao_nr())
    rif.conv_tol = 1e-13
    rif.kernel()
    if not rif.converged:
        raise RuntimeError(name + ': reference RI-HF did not converge')
    fci_solver = fci.FCI(hf)
    root_energies = np.asarray(fci_solver.kernel(nroots=3)[0])
    if not np.all(fci_solver.converged):
        raise RuntimeError(name + ': reference FCI roots did not converge')
    mo_eri = ao2mo.incore.full(eri, hf.mo_coeff, compact=False).reshape([molecule.nao_nr()]*4)
    embedded = dict(schema='numivivo.org/embedded-hamiltonian/v1',
        orbitalIdentifiers=[name+'-'+str(i) for i in range(molecule.nao_nr())],
        alphaElectrons=molecule.nelec[0], betaElectrons=molecule.nelec[1],
        oneElectron=matrix(hf.mo_coeff.T@core@hf.mo_coeff), twoElectron=mo_eri.ravel().tolist(),
        constantEnergyHartree=molecule.energy_nuc(), energyReference='Born-Oppenheimer electronic; nuclear scalar included',
        provenance=dict(oracle='PySCF 2.8.0'))
    reference = dict(hf=hf.e_tot, mp2=mp.MP2(hf).run().e_tot, ccsd=coupled.e_tot,
        fciRoots=root_energies.tolist(), riHF=rif.e_tot, riMP2=mp.MP2(rif).run().e_tot, riRank=int(retained.sum()))
    if name == 'lih':
        cas = mcscf.CASSCF(hf, 2, 2).state_average_([0.5, 0.5])
        cas.conv_tol = 1e-11
        cas.conv_tol_grad = 1e-6
        cas.max_cycle_macro = 100
        cas.kernel()
        if not cas.converged:
            raise RuntimeError('Independent state-average did not converge')
        reference.update(saConverged=True, saEnergy=cas.e_tot, saStates=np.asarray(cas.e_states).tolist())
    if name in ('h2', 'water'):
        solvent = pcm.PCM(molecule)
        solvent.method = 'C-PCM'
        solvent.eps = 78.3
        solvent.lebedev_order = 29
        radii = np.ones(119)
        radii[1], radii[8] = 1.2, 1.52
        solvent.radii_table = radii*1.2/0.529177210544
        solvated = scf.RHF(molecule).PCM(solvent)
        solvated.conv_tol = 1e-12
        solvated.kernel()
        if not solvated.converged:
            raise RuntimeError('Independent solvent SCF did not converge')
        reference.update(smoothCPCM=solvated.e_tot,
                         smoothPolarization=solvent._get_vind(solvated.make_rdm1())[0])
    system = dict(nuclei=[dict(atomicNumber=int(z), positionBohr=r.tolist())
        for z, r in zip(molecule.atom_charges(), molecule.atom_coords())], pointCharges=[],
        alphaElectrons=molecule.nelec[0], betaElectrons=molecule.nelec[1])
    return dict(identifier=name, system=system, basis=basis(molecule, 'sto-3g'),
        auxiliaryBasis=basis(auxiliary, 'weigend'), overlap=matrix(overlap*scale[:, None]*scale[None, :]),
        coreHamiltonian=matrix(core*scale[:, None]*scale[None, :]), fittedERI=fitted.ravel().tolist(),
        embedded=embedded, reference=reference)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    if pyscf.__version__ != '2.8.0':
        raise RuntimeError('Expected pinned PySCF 2.8.0')
    pyscf.lib.num_threads(1)
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    systems = [('h2', 'H 0 0 -0.7; H 0 0 0.7'), ('lih', 'Li 0 0 0; H 0 0 3.0'),
        ('water', 'O 0 0 0; H 0 -1.43233673 1.10715266; H 0 1.43233673 1.10715266')]
    for name, atom in systems:
        data = generate(name, atom)
        path = args.output/(name+'.json')
        path.write_text(json.dumps(data, sort_keys=True, allow_nan=False)+'\n')
        records.append(dict(file=path.name, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        print(name, data['reference'], flush=True)
    manifest = dict(schema='numivivo.org/external-conformance/v1', oracle='PySCF', oracleVersion=pyscf.__version__,
        numpy=np.__version__, python=sys.version, platform=platform.platform(), fixtures=records)
    (args.output/'manifest.json').write_text(json.dumps(manifest, sort_keys=True, indent=2)+'\n')


if __name__ == '__main__':
    main()
