#!/usr/bin/env python3
"""Independent fixed-frame reconstruction of native ECC-DMET observations.

PySCF contracts the global Hamiltonian; a projected CISD reference is obtained
by restricting PySCF's independently generated FCI action. NumPy performs its
own bath SVD and tensor projections. Native state coefficients, Hamiltonians,
RDMs and energies are never inputs to the expected-energy calculation. Returned
orbital rotations and fitted potentials specify the fixed point to verify.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import hashlib
import json
import pathlib
import sys
import numpy as np
import pyscf
from pyscf import fci


def matrix(value):
    return np.array(value['values'], dtype=float).reshape(value['rows'], value['columns'])


def reference(one, eri, nelec, method):
    n = len(one)
    alpha = fci.cistring.make_strings(range(n), nelec[0])
    beta = fci.cistring.make_strings(range(n), nelec[1])
    dimension = len(alpha) * len(beta)
    if dimension > 1024:
        raise ValueError('Independent dense ECC oracle is deliberately bounded')
    absorbed = fci.direct_spin1.absorb_h1e(one, eri, n, nelec, .5)
    h = np.zeros((dimension, dimension))
    for column in range(dimension):
        vector = np.zeros((len(alpha), len(beta)))
        vector.ravel()[column] = 1
        h[:, column] = fci.direct_spin1.contract_2e(absorbed, vector, n, nelec).ravel()
    if np.max(np.abs(h-h.T)) > 1e-10:
        raise ValueError('Independent Hamiltonian is not Hermitian')
    if method == 'cisd':
        a0, b0 = (1 << nelec[0])-1, (1 << nelec[1])-1
        selected = [i*len(beta)+j for i, a in enumerate(alpha) for j, b in enumerate(beta)
                    if ((int(a)^a0).bit_count()+(int(b)^b0).bit_count())//2 <= 2]
    elif method == 'fci':
        selected = list(range(dimension))
    else:
        raise ValueError('Unknown reference method')
    energies, vectors = np.linalg.eigh(h[np.ix_(selected, selected)])
    state = np.zeros(dimension)
    state[selected] = vectors[:, 0]
    dm1, dm2 = fci.direct_spin1.make_rdm12(state.reshape(len(alpha),len(beta)), n, nelec)
    evaluated = np.einsum('pq,pq', one, dm1)+.5*np.einsum('pqrs,pqrs', eri, dm2)
    if abs(evaluated-energies[0]) > 1e-9:
        raise ValueError('Independent RDM convention does not reconstruct its energy')
    return energies[0], dm1, dm2


def transform(g, c):
    return np.einsum('pqrs,pi,qj,rk,sl->ijkl', g, c, c, c, c, optimize=True)


def compare(directory):
    if pyscf.__version__ != '2.8.0':
        raise ValueError('PySCF 2.8.0 is required')
    root = pathlib.Path(directory)
    checks = []
    def check(label, actual, expected, tolerance=1e-8):
        error = abs(float(actual)-float(expected))
        passed = bool(np.isfinite(error) and error <= tolerance)
        checks.append(dict(observable=label, actual=float(actual), expected=float(expected), error=error,
                           tolerance=tolerance, passed=passed))
        print(('PASS' if passed else 'FAIL'), label, error)
    hashes = {}
    for name, hname in [('partition', 'hamiltonian'), ('environment', 'environment-hamiltonian')]:
        source = root / (hname+'.json')
        output = root / (name+'-result.json')
        for path in [source,output]:
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        h, observed = json.loads(source.read_text()), json.loads(output.read_text())
        cfg = observed['configuration']
        n = len(h['orbitalIdentifiers'])
        u = matrix(observed['orbitalRotation'])
        physical_one = u.T@matrix(h['oneElectron'])@u
        physical_two = transform(np.array(h['twoElectron']).reshape([n]*4), u)
        ref_one, ref_two = physical_one.copy(), physical_two.copy()
        for amplitude, operator in zip(observed['correlationPotentialHartree'], cfg['correlationOperators']):
            ref_one += amplitude*matrix(operator['one'])
            ref_two += amplitude*np.array(operator['two']).reshape([n]*4)
        nelec = (h['alphaElectrons'], h['betaElectrons'])
        _, dm1, dm2 = reference(ref_one,ref_two,nelec,cfg['matching']['referenceMethod'])
        base = h['constantEnergyHartree']+np.einsum('pq,pq',physical_one,dm1)+.5*np.einsum('pqrs,pqrs',physical_two,dm2)
        check(name+' reference physical energy',observed['frame']['referencePhysicalEnergyHartree'],base)
        correction, populations = 0.0, 0.0
        fragment_moments = np.zeros(len(cfg['correlationOperators']))
        for index, fragment in enumerate(cfg['fragments']):
            o = observed['frame']['clusters'][index]
            c = matrix(o['coefficients'])
            f = fragment['orbitals']; environment = [i for i in range(n) if i not in f]
            # Left SVD vectors, not the right vectors, span the environment bath.
            left, singular, _ = np.linalg.svd(dm1[np.ix_(environment,f)],full_matrices=False)
            remaining = float(np.sum(singular**2)); rank = 0
            while remaining > cfg['matching']['bathDiscardedWeight'] and rank < min(fragment['maximumBathOrbitals'],len(singular)):
                if singular[rank]**2 <= 1e-13*max(1,float(np.sum(singular**2))):
                    break
                remaining = max(0,remaining-singular[rank]**2); rank += 1
            check(name+' spectral bath rank '+str(index),o['spectralBathRank'],rank,0)
            if rank:
                embedded_bath = np.zeros((n,rank)); embedded_bath[environment,:] = left[:,:rank]
                native_bath = c[:,len(f):len(f)+rank]
                defect = np.linalg.norm(native_bath@native_bath.T-embedded_bath@embedded_bath.T)
                check(name+' spectral bath projector '+str(index),defect,0,1e-7)
            q = np.eye(n)-c@c.T; env_density=q@dm1@q
            field = np.einsum('rs,pqrs->pq',env_density,physical_two)-.5*np.einsum('rs,prqs->pq',env_density,physical_two)
            bare=c.T@physical_one@c; local_field=c.T@field@c; effective=bare+local_field
            eri=transform(physical_two,c); k=c.shape[1]
            check(name+' environment field '+str(index),np.max(np.abs(local_field-matrix(o['environmentPotential']))),0)
            mu=observed['frame'].get('numberMatching',{}).get('chemicalPotentialHartree',0)
            projector=np.zeros((k,k));projector[range(len(f)),range(len(f))]=1
            _, high_one, high_two = reference(effective-mu*projector,eri,(fragment['clusterAlphaElectrons'],fragment['clusterBetaElectrons']),'fci')
            populations += np.einsum('pq,pq',projector,high_one)
            physical_cluster = np.einsum('pq,pq',effective,high_one)+.5*np.einsum('pqrs,pqrs',eri,high_two)
            check(name+' physical impurity '+str(index),observed['frame']['states'][index]['physicalClusterEnergyHartree'],physical_cluster)
            low_one=c.T@dm1@c; low_two=transform(dm2,c)
            if cfg['mode']=='selfConsistentPartition':
                w=np.zeros(k);w[:len(f)]=1
                energy_one=.5*(w[:,None]+w[None,:])*(bare+.5*local_field)
                energy_two=.25*(w[:,None,None,None]+w[None,:,None,None]+w[None,None,:,None]+w[None,None,None,:])*eri
            else:
                energy_one,energy_two=effective,eri
            low=np.einsum('pq,pq',energy_one,low_one)+.5*np.einsum('pqrs,pqrs',energy_two,low_two)
            high=np.einsum('pq,pq',energy_one,high_one)+.5*np.einsum('pqrs,pqrs',energy_two,high_two)
            reported=observed['frame']['energyContributions'][index]
            check(name+' reference contribution '+str(index),reported['referenceContributionHartree'],low)
            check(name+' high-level contribution '+str(index),reported['impurityContributionHartree'],high)
            correction += high-low
            for j, operator in enumerate(cfg['correlationOperators']):
                op1=matrix(operator['one'])[np.ix_(f,f)]
                op2=np.array(operator['two']).reshape([n]*4)[np.ix_(f,f,f,f)]
                fragment_moments[j] += np.einsum('pq,pq',op1,high_one[:len(f),:len(f)])+.5*np.einsum('pqrs,pqrs',op2,high_two[:len(f),:len(f),:len(f),:len(f)])
        check(name+' restored total energy',observed['frame']['energyHartree'],base+correction)
        if cfg['mode']=='selfConsistentPartition':
            check(name+' physical particle equation',populations,sum(nelec),cfg['matching']['numberMatching']['populationTolerance']*1.1)
            for i,operator in enumerate(cfg['correlationOperators']):
                value=np.einsum('pq,pq',matrix(operator['one']),dm1)+.5*np.einsum('pqrs,pqrs',np.array(operator['two']).reshape([n]*4),dm2)
                check(name+' matched moment '+str(i),value,fragment_moments[i],cfg['matching']['momentTolerance']*1.1)
    report=dict(schema='numivivo.org/ecc-external-conformance/v1',oracle='PySCF',oracleVersion=pyscf.__version__,
                passed=bool(checks and all(c['passed'] for c in checks)), checks=checks,inputSHA256=hashes,
                scope='Independent fixed-point and energy reconstruction; no assertion of globally N-representable patched RDMs or a validated reaction barrier')
    (root/'independent-results.json').write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n')
    if not report['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    if len(sys.argv)!=2:
        raise SystemExit('usage: compare_ecc_result.py native-ecc-results-directory')
    compare(sys.argv[1])
