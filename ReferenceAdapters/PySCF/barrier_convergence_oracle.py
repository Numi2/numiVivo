#!/usr/bin/env python3
"""Independent PySCF reconstruction from declared molecular inputs, never native energies.

export: pinned PySCF 2.8.0, normalized Cartesian integrals, Gaussian cross overlaps,
parallel transport and direct FCI/CASCI at every specified geometry and level.
check: standard-library-only production CLI/verified artifact cache comparisons.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')


def load_request(path):
    doc = json.loads(path.read_text())
    if doc['schema'] != 'numivivo.org/reaction-calculation/v1' or set(doc['calculation']) != {'barrierConvergence'}:
        raise ValueError('Expected the explicit barrier-convergence reaction request')
    request = doc['calculation']['barrierConvergence']['request']
    if request['schema'] != 'numivivo.org/barrier-convergence/v1':
        raise ValueError('Unsupported campaign schema')
    return doc, request


def export(request_file, out):
    import numpy as np
    import pyscf
    from scipy.linalg import eigh
    from pyscf import gto, ao2mo, fci
    if pyscf.__version__ != '2.8.0':
        raise RuntimeError('This conformance protocol pins PySCF 2.8.0')
    pyscf.lib.num_threads(1)
    _, r = load_request(request_file)
    # General Cartesian scaling supported, but fixed-spherical, point-charge and
    # ECC oracles are separate methods and cannot be silently substituted here.
    if r['basis']['representation'] != 'normalized-cartesian' or any(p['system']['pointCharges'] for p in r['snapshots']):
        raise ValueError('Oracle requires explicit normalized Cartesian gas-phase inputs')
    if any(set(level['method']) != {'casci'} for level in r['levels']):
        raise ValueError('This independent adapter checks CAS ladders, not ECC through a CAS substitute')
    molecules, scales, overlaps, cores, coefficients = [], [], [], [], []
    for point in r['snapshots']:
        system = point['system']
        atom_list, basis = [], {}
        for i, nucleus in enumerate(system['nuclei']):
            symbol = gto.mole._atom_symbol(nucleus['atomicNumber'])
            label = symbol + str(i)
            atom_list.append((label, nucleus['positionBohr']))
            basis[label] = [[s['angularMomentum'], *[[p['exponent'], p['coefficient']] for p in s['primitives']]]
                            for s in r['basis']['shells'] if s['nucleusIndex'] == i]
        electrons = system['alphaElectrons'] + system['betaElectrons']
        m = gto.M(atom=atom_list, basis=basis, unit='Bohr', cart=True, verbose=0,
                  charge=sum(n['atomicNumber'] for n in system['nuclei'])-electrons,
                  spin=system['alphaElectrons']-system['betaElectrons'])
        raw_s = m.intor('int1e_ovlp')
        scale = 1/np.sqrt(np.diag(raw_s))
        s = raw_s * scale[:,None] * scale[None,:]
        h = (m.intor('int1e_kin')+m.intor('int1e_nuc')) * scale[:,None] * scale[None,:]
        _, c = eigh(h, s)
        molecules.append(m); scales.append(scale); overlaps.append(s); cores.append(h); coefficients.append(c)
    anchor = next(i for i,p in enumerate(r['snapshots']) if p['identifier'] == r['anchorPointIdentifier'])
    candidates = [c.copy() for c in coefficients]
    if 'anchorCoefficients' in r and r['anchorCoefficients'] is not None:
        item = r['anchorCoefficients']; coefficients[anchor] = np.array(item['values']).reshape(item['rows'], item['columns'])
    if np.linalg.norm(coefficients[anchor].T @ overlaps[anchor] @ coefficients[anchor]-np.eye(molecules[anchor].nao)) > 1e-8:
        raise ValueError('Oracle anchor is not orthonormal')
    cross = [gto.intor_cross('int1e_ovlp', a, b)*sa[:,None]*sb[None,:]
             for a,b,sa,sb in zip(molecules,molecules[1:],scales,scales[1:])]
    minima = [1.0]*len(molecules)
    def align(source, dest, s):
        matrix = coefficients[source].T @ s @ candidates[dest]
        rotation = np.eye(matrix.shape[0]); minimum = 1.0
        for group in r['transportGroups']:
            u, singular, vh = np.linalg.svd(matrix[np.ix_(group, group)])
            minimum = min(minimum, float(singular[-1]))
            if singular[-1] < r['minimumTransportSingularValue']:
                raise ValueError('Independent oracle also lost the declared transport subspace')
            rotation[np.ix_(group, group)] = vh.T @ u.T
        coefficients[dest] = candidates[dest] @ rotation; minima[dest] = minimum
    for i in range(anchor-1, -1, -1): align(i+1, i, cross[i].T)
    for i in range(anchor+1, len(molecules)): align(i-1, i, cross[i-1])
    packed = []
    full, levels, states = [], {x['identifier']: [] for x in r['levels']}, []
    for i, m in enumerate(molecules):
        c = coefficients[i]
        h = c.T @ cores[i] @ c
        raw_c = scales[i][:,None] * c
        eri = ao2mo.restore(1, ao2mo.kernel(m, raw_c), m.nao)
        def solve(active, core):
            ha = h[np.ix_(active,active)].copy()
            for p_idx,p in enumerate(active):
                for q_idx,q in enumerate(active):
                    ha[p_idx,q_idx] += sum(2*eri[p,q,a,a]-eri[p,a,a,q] for a in core)
            constant = m.energy_nuc()+sum(2*h[a,a] for a in core)
            constant += sum(2*eri[a,a,b,b]-eri[a,b,b,a] for a in core for b in core)
            engine = fci.direct_spin1.FCI(); engine.conv_tol = 1e-13
            ne = tuple(n-len(core) for n in m.nelec)
            energy, vector = engine.kernel(ha, eri[np.ix_(active,active,active,active)], len(active), ne, ecore=constant)
            if not engine.converged or not np.isfinite(energy): raise RuntimeError('Oracle CI did not converge')
            return float(energy), vector
        energy, vector = solve(list(range(m.nao)), [])
        full.append(energy); states.append(vector); packed.append((m,h,eri))
    occupations = None
    rotation = np.eye(molecules[0].nao)
    if r.get('ensembleOrbitals') is not None:
        policy = r['ensembleOrbitals']
        density = sum(weight*fci.direct_spin1.make_rdm1(vector,m.nao,m.nelec)
                      for weight,vector,m in zip(policy['pointWeights'],states,molecules))
        occupations, rotation = eigh(density)
        occupations = occupations[::-1].tolist(); rotation = rotation[:,::-1]
    for m, h0, eri0 in packed:
        h = rotation.T @ h0 @ rotation
        eri = ao2mo.incore.full(eri0,rotation,compact=False).reshape([m.nao]*4)
        for level in r['levels']:
            p = level['method']['casci']['partition']
            energy, _ = solve(p['active'], p['doublyOccupiedCore'])
            levels[level['identifier']].append(energy)
    ref_overlaps = []
    for i,s in enumerate(cross):
        orbital_overlap = coefficients[i].T @ s @ coefficients[i+1]
        value = fci.addons.overlap(states[i], states[i+1], molecules[i].nao, molecules[i].nelec, s=orbital_overlap)
        ref_overlaps.append(float(abs(value)**2))
    b = next(i for i,p in enumerate(r['snapshots']) if p['identifier']==r['barrierPointIdentifier'])
    def diff(e): return dict(forwardHartree=e[b]-e[0], reverseHartree=e[b]-e[-1], reactionHartree=e[-1]-e[0], relativeProfileHartree=[v-e[0] for v in e])
    record = dict(schema='numivivo.org/barrier-oracle/v1', oracle='PySCF', version=pyscf.__version__,
        pointIdentifiers=[p['identifier'] for p in r['snapshots']], referenceEnergiesHartree=full,
        ensembleOccupations=occupations,
        referenceDifferences=diff(full), transportMinimumSingularValues=minima, adjacentReferenceOverlapsSquared=ref_overlaps,
        levels=[dict(identifier=k,energiesHartree=e,differences=diff(e)) for k,e in levels.items()],
        sourceRequestSHA256=hashlib.sha256(request_file.read_bytes()).hexdigest())
    write(out/'reference.json',record)
    write(out/'manifest.json',dict(referenceSHA256=hashlib.sha256((out/'reference.json').read_bytes()).hexdigest(),
                                   sourceRequestSHA256=record['sourceRequestSHA256'], pointCount=len(full),levelCount=len(levels)))
    print('Generated independent barrier/reference/transport values for',len(full),'geometries')


def check(binary, request_file, oracle_dir, out):
    import copy
    import subprocess
    if out.exists() and any(out.iterdir()): raise RuntimeError('Use an empty output directory; artifacts are not deleted')
    out.mkdir(parents=True,exist_ok=True)
    document, request = load_request(request_file)
    raw = (oracle_dir/'reference.json').read_bytes(); manifest=json.loads((oracle_dir/'manifest.json').read_text()); ref=json.loads(raw)
    if hashlib.sha256(raw).hexdigest()!=manifest['referenceSHA256'] or hashlib.sha256(request_file.read_bytes()).hexdigest()!=manifest['sourceRequestSHA256']:
        raise RuntimeError('Input or independent reference hash mismatch')
    if ref['version']!='2.8.0' or ref['schema']!='numivivo.org/barrier-oracle/v1': raise RuntimeError('Oracle method/schema mismatch')
    if manifest['pointCount']!=len(request['snapshots']) or manifest['levelCount']!=len(request['levels']): raise RuntimeError('Missing oracle entries')
    if ref['pointIdentifiers']!=[p['identifier'] for p in request['snapshots']] or [p['identifier'] for p in ref['levels']]!=[p['identifier'] for p in request['levels']]:
        raise RuntimeError('Missing or reordered point/level oracle identities')
    checks=[]
    def require(value,label,error=None):
        item=dict(label=label,passed=bool(value))
        if error is not None: item['absoluteError']=error
        checks.append(item)
        write(out/'checks.json',dict(checks=checks,passed=all(c['passed'] for c in checks)))
        if not value: raise RuntimeError(label)
    def compare(a,b,label,tol=1e-8):
        error=abs(a-b); require(math.isfinite(error) and error<=tol,label,error)
    def invoke(doc,name,fail=False):
        path=out/(name+'.request.json'); dest=out/(name+'.result.json');write(path,doc)
        p=subprocess.run([str(binary),'reaction-run',str(path),'--store',str(out/'store'),'--output',str(dest)],capture_output=True,text=True,timeout=600)
        (out/(name+'.log')).write_text(p.stdout+'\n'+p.stderr)
        if fail:
            require(p.returncode!=0 and not dest.exists() and not Path(str(dest)+'.receipt.json').exists(),name);return None
        if p.returncode: raise RuntimeError(p.stderr)
        return json.loads(dest.read_text())['barrierConvergence']['result'],json.loads(Path(str(dest)+'.receipt.json').read_text())
    template = 'h3-barrier-convergence-ensemble' if request.get('ensembleOrbitals') is not None else 'h3-barrier-convergence'
    generated = out/'generated-template.json'
    process = subprocess.run([str(binary),'reaction-template',template,'--output',str(generated)],capture_output=True,text=True,timeout=60)
    (out/'template.log').write_text(process.stdout+'\n'+process.stderr)
    require(process.returncode == 0 and generated.exists(), 'production template command succeeds')
    require(json.loads(generated.read_text()) == document, 'committed inputs match the native template exactly')
    result, receipt=invoke(document,'first')
    require(not receipt['reused'],'fresh bounded campaign executes on the production CLI')
    require(result['request']==request,'campaign result binds the exact declared request')
    require(len(result['referenceEnergiesHartree'])==len(ref['referenceEnergiesHartree']) and len(result['levels'])==len(ref['levels'])
            and len(result['transportMinimumSingularValues'])==manifest['pointCount']
            and len(result['adjacentReferenceOverlapsSquared'])==manifest['pointCount']-1,'native result contains every point, level and continuity record')
    if ref['ensembleOccupations'] is not None:
        require(len(result.get('ensembleOccupations',[]))==len(ref['ensembleOccupations']),'all ensemble occupations retained')
        for i,(a,b) in enumerate(zip(result['ensembleOccupations'],ref['ensembleOccupations'])): compare(a,b,f'ensemble occupation {i}')
    else: require(result.get('ensembleOccupations') is None,'plain transport does not claim ensemble-derived orbitals')
    for i,(a,b) in enumerate(zip(result['referenceEnergiesHartree'],ref['referenceEnergiesHartree'])): compare(a,b,f'full FCI energy {i}')
    for key in ('forwardHartree','reverseHartree','reactionHartree'):
        compare(result['referenceDifferences'][key],ref['referenceDifferences'][key],'reference '+key)
    for i,(a,b) in enumerate(zip(result['transportMinimumSingularValues'],ref['transportMinimumSingularValues'])): compare(a,b,f'physical orbital transport {i}')
    for i,(a,b) in enumerate(zip(result['adjacentReferenceOverlapsSquared'],ref['adjacentReferenceOverlapsSquared'])): compare(a,b,f'physical CI state overlap {i}')
    for level, reference in zip(result['levels'],ref['levels']):
        require(level['identifier']==reference['identifier'] and 'evaluated' in level['outcome'],level['identifier']+' independently completed')
        value=level['outcome']['evaluated']['result']
        require(len(value['energiesHartree'])==len(ref['referenceEnergiesHartree']),level['identifier']+' all energies present')
        for i,(a,b) in enumerate(zip(value['energiesHartree'],reference['energiesHartree'])): compare(a,b,f"{level['identifier']} energy {i}")
        for key in ('forwardHartree','reverseHartree','reactionHartree'):
            compare(value['differences'][key],reference['differences'][key],level['identifier']+' '+key)
    # The input fixes the accuracy threshold before execution. All deliberately
    # naive truncated spaces fail for the plain frame. With ensemble orbitals,
    # one reduced level passes but adjacent-level stability still fails. The
    # full-space endpoint cannot substitute for that missing evidence.
    require(result['assessment']=='reducedAccuracyNotEstablished' and result.get('acceptedReducedLevelIdentifier') is None,
            'full-space success does not certify inaccurate reduced models')
    if request.get('ensembleOrbitals') is None:
        require(all(not l['outcome']['evaluated']['result']['meetsReferenceAccuracy'] for l in result['levels'][:-1]),'all inaccurate reduced levels remain explicitly failed accuracy assessments')
    else:
        require(result['levels'][-2]['outcome']['evaluated']['result']['meetsReferenceAccuracy']
                and not result['levels'][-3]['outcome']['evaluated']['result']['meetsReferenceAccuracy'],
                'improved five-orbital accuracy does not bypass the adjacent-level stability requirement')
    require(result['levels'][-1]['outcome']['evaluated']['result']['meetsReferenceAccuracy'] and not result['levels'][-1]['outcome']['evaluated']['result']['isGenuinelyReduced'],
            'full-space endpoint is identified separately')
    require('not a saddle' in result['meaning'] and 'rate' in result['meaning'],'assessment retains nuclear/kinetic qualification boundary')
    repeated, second=invoke(document,'cached')
    require(second['reused'] and second['result']==receipt['result'] and repeated==result,'verified cache reconstructs identical negative scientific assessment')
    changed=copy.deepcopy(document); changed['calculation']['barrierConvergence']['request']['levels'][1]['method']['casci']['partition']['active']=[0,1,2]
    invoke(changed,'duplicate-active-space-rejected',True)
    changed=copy.deepcopy(document); changed['calculation']['barrierConvergence']['request']['maximumPointEvaluations']=1
    invoke(changed,'unfunded-reference-and-ladder-rejected',True)
    changed=copy.deepcopy(document); changed['calculation']['barrierConvergence']['request']['snapshots'][0]['system']['betaElectrons']=0
    invoke(changed,'electron-sector-change-rejected',True)
    digest=bytes(receipt['result']['bytes']).hex();object_path=out/'store'/'objects'/'sha256'/digest[:2]/digest[2:4]/digest;stored=object_path.read_bytes()
    if hashlib.sha256(stored).hexdigest()!=digest: raise RuntimeError('Unexpected artifact layout')
    try:
        object_path.write_bytes(stored+b'\n');invoke(document,'corrupt-campaign-artifact-rejected',True)
    finally: object_path.write_bytes(stored)
    write(out/'checks.json',dict(schema='numivivo.org/barrier-conformance/v1',checks=checks,passed=True,
        executableSHA256=hashlib.sha256(binary.read_bytes()).hexdigest(), referenceSHA256=manifest['referenceSHA256'],
        scientificAssessment=result['assessment'], scope='mapped H3/6-31G fixed-geometry energy profile and deliberately inaccurate nested CAS spaces; not a stationary barrier, paper reproduction or rate'))
    print('PASS',len(checks),'independent barrier comparisons and production CLI checks')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('export');p.add_argument('request',type=Path);p.add_argument('out',type=Path)
    p=sub.add_parser('check');p.add_argument('binary',type=Path);p.add_argument('request',type=Path);p.add_argument('reference',type=Path);p.add_argument('out',type=Path)
    a=parser.parse_args()
    if a.command=='export':export(a.request,a.out)
    else:check(a.binary.resolve(),a.request,a.reference,a.out)
