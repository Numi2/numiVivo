#!/usr/bin/env python3
"""Independent bounded variational-space references, not production dependencies.

Export reads declared inputs, reconstructs Gaussian integrals with PySCF 2.8.0,
then independently implements the shared projected eigenproblem and correlated
solvent functional with NumPy/SciPy. No native energies or states enter export.
Check uses only the standard library and the real production executable.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')


def load(path, kind):
    doc = json.loads(path.read_text())
    if doc['schema'] != 'numivivo.org/reaction-calculation/v1' or set(doc['calculation']) != {kind}:
        raise ValueError('Unexpected calculation input')
    return doc['calculation'][kind]['request']


def full_matrix(h, eri, n, ne, constant):
    import numpy as np
    from pyscf import fci
    a = fci.cistring.make_strings(range(n), ne[0])
    b = fci.cistring.make_strings(range(n), ne[1])
    shape = (len(a), len(b))
    effective = fci.direct_spin1.absorb_h1e(h, eri, n, ne, 0.5)
    columns = []
    for k in range(len(a)*len(b)):
        unit = np.zeros(shape); unit.flat[k] = 1
        columns.append(fci.direct_spin1.contract_2e(effective, unit, n, ne).reshape(-1))
    matrix = np.column_stack(columns)+np.eye(len(columns))*constant
    if np.linalg.norm(matrix-matrix.T)>1e-9:
        raise RuntimeError('Independent Hamiltonian is not Hermitian')
    return matrix, a, b


def seed_indices(partition, alpha, beta, n):
    core, active = partition['doublyOccupiedCore'], partition['active']
    if partition['frozenOrbitals'] or set(core) & set(active):
        raise ValueError('Invalid fixed subspace partition')
    fixed = sum(1 << i for i in core)
    allowed = fixed | sum(1 << i for i in active)
    return [i*len(beta)+j for i,a in enumerate(alpha) for j,b in enumerate(beta)
            if int(a)&fixed == fixed and int(b)&fixed == fixed
            and int(a)&~allowed == 0 and int(b)&~allowed == 0]



def barrier_workspace(r):
    """Rebuild the common frame and full reference from declared physical inputs."""
    import numpy as np
    from scipy.linalg import eigh
    from pyscf import gto, ao2mo, fci
    molecules, overlaps, cores, candidates = [], [], [], []
    # These benchmark shells are normalized s functions; reject another
    # representation rather than approximating its convention.
    if r['basis']['representation']!='normalized-cartesian' or any(s['angularMomentum']!=0 for s in r['basis']['shells']):
        raise ValueError('This residual fixture oracle requires explicit Cartesian s shells')
    for point in r['snapshots']:
        system=point['system'];atoms=[];basis={}
        if system['pointCharges']: raise ValueError('External-charge oracle not enabled')
        for i,a in enumerate(system['nuclei']):
            label=gto.mole._atom_symbol(a['atomicNumber'])+str(i);atoms.append((label,a['positionBohr']))
            basis[label]=[[s['angularMomentum'],*[[v['exponent'],v['coefficient']] for v in s['primitives']]]
                          for s in r['basis']['shells'] if s['nucleusIndex']==i]
        m=gto.M(atom=atoms,basis=basis,unit='Bohr',cart=True,verbose=0,
                charge=sum(a['atomicNumber'] for a in system['nuclei'])-system['alphaElectrons']-system['betaElectrons'],
                spin=system['alphaElectrons']-system['betaElectrons'])
        overlap=m.intor('int1e_ovlp');h=m.intor('int1e_kin')+m.intor('int1e_nuc')
        _,c=eigh(h,overlap);molecules.append(m);overlaps.append(overlap);cores.append(h);candidates.append(c)
    frames=[c.copy() for c in candidates]
    anchor=next(i for i,p in enumerate(r['snapshots']) if p['identifier']==r['anchorPointIdentifier'])
    if r.get('anchorCoefficients') is not None:
        c=r['anchorCoefficients'];frames[anchor]=np.array(c['values']).reshape(c['rows'],c['columns'])
    cross=[gto.intor_cross('int1e_ovlp',a,b) for a,b in zip(molecules,molecules[1:])]
    def align(a,b,s):
        overlap=frames[a].T@s@candidates[b];rotation=np.eye(molecules[a].nao)
        for group in r['transportGroups']:
            u,singular,vh=np.linalg.svd(overlap[np.ix_(group,group)])
            if singular[-1]<r['minimumTransportSingularValue']: raise ValueError('Physical transport rank lost')
            rotation[np.ix_(group,group)]=vh.T@u.T
        frames[b]=candidates[b]@rotation
    for i in range(anchor-1,-1,-1):align(i+1,i,cross[i].T)
    for i in range(anchor+1,len(frames)):align(i-1,i,cross[i-1])
    packed=[];energies=[];density=np.zeros((molecules[0].nao,)*2)
    for i,(m,c) in enumerate(zip(molecules,frames)):
        h=c.T@cores[i]@c;g=ao2mo.restore(1,ao2mo.kernel(m,c),m.nao)
        engine=fci.direct_spin1.FCI();engine.conv_tol=1e-13
        e,v=engine.kernel(h,g,m.nao,m.nelec,ecore=m.energy_nuc())
        if not engine.converged:raise RuntimeError('Unconverged independent reference')
        packed.append((m,h,g));energies.append(float(e))
        if r.get('ensembleOrbitals') is not None:
            density+=r['ensembleOrbitals']['pointWeights'][i]*engine.make_rdm1(v,m.nao,m.nelec)
    rotation=np.eye(molecules[0].nao)
    if r.get('ensembleOrbitals') is not None:
        _,rotation=eigh(density);rotation=rotation[:,::-1]
    return dict(packed=packed,rotation=rotation,reference=dict(referenceEnergiesHartree=energies))


def residual_reference(path, out):
    import numpy as np
    from scipy.linalg import eigh
    from pyscf import ao2mo
    r = load(path, 'residualBarrier')
    workspace = barrier_workspace(r['baseline'])
    matrices, alpha, beta = [], None, None
    rotation = workspace['rotation']
    for m, h0, eri0 in workspace['packed']:
        h = rotation.T @ h0 @ rotation
        eri = ao2mo.incore.full(eri0, rotation, compact=False).reshape([m.nao]*4)
        matrix, alpha, beta = full_matrix(h, eri, m.nao, m.nelec, m.energy_nuc())
        matrices.append(matrix)
    dimension = len(alpha)*len(beta)
    indices = seed_indices(r['seed'], alpha, beta, matrices[0].shape[0])
    columns = [np.eye(dimension)[:,i] for i in indices]
    rank_limit = r['space']['maximumDimension']
    rank_tolerance = r['space']['linearDependenceTolerance']
    residual_tolerance = r['space']['projectedResidualToleranceHartree']
    truth = np.array(workspace['reference']['referenceEnergiesHartree'])
    points = r['baseline']['snapshots']
    barrier = next(i for i,p in enumerate(points) if p['identifier']==r['baseline']['barrierPointIdentifier'])
    tprofile = truth-truth[0]
    levels = []
    added = len(columns)
    def outside(v, space):
        v = v.copy()
        for _ in range(2):
            for q in space: v -= np.dot(q,v)*q
        return v
    for step in range(r['maximumRefinementRounds']+1):
        b = np.column_stack(columns)
        if np.linalg.norm(b.T@b-np.eye(len(columns)))>1e-9:
            raise RuntimeError('Independent shared space lost rank')
        energies, residuals, projected, external = [], [], [], []
        for h in matrices:
            e, v = eigh(b.T@h@b)
            vector = b@v[:,0]
            energy = float(vector@h@vector)
            residual = h@vector-energy*vector
            energies.append(energy)
            projected.append(float(np.linalg.norm(b.T@residual)))
            ext = outside(residual, columns)
            external.append(float(np.linalg.norm(ext)))
            residuals.append(ext)
        energies = np.array(energies)
        profile = energies-energies[0]
        barrier_error = max(abs(profile[barrier]-tprofile[barrier]),
                            abs(energies[barrier]-energies[-1]-truth[barrier]+truth[-1]),
                            abs(profile[-1]-tprofile[-1]))
        successive = None
        if levels:
            old = np.array(levels[-1]['energiesHartree'])
            successive = float(max(np.max(abs(profile-old+old[0])),
                                   abs(energies[barrier]-energies[-1]-old[barrier]+old[-1])))
        levels.append(dict(round=step, variationalDimension=len(columns), fullSectorDimension=dimension,
            newDirections=added, energiesHartree=energies.tolist(),
            maximumBarrierErrorHartree=float(barrier_error), maximumProfileErrorHartree=float(max(abs(profile-tprofile))),
            maximumAbsoluteEnergyErrorHartree=float(max(abs(energies-truth))), maximumSuccessiveChangeHartree=successive,
            projectedResidualHartree=projected, externalResidualHartree=external))
        if step==r['maximumRefinementRounds'] or len(columns)==dimension: break
        before = len(columns)
        # All residuals were formed in the same previous projector.
        for residual in residuals:
            length = np.linalg.norm(residual)
            if length<=residual_tolerance: continue
            q = outside(residual/length, columns); norm = np.linalg.norm(q)
            if norm<=rank_tolerance: continue
            if len(columns)>=rank_limit: raise RuntimeError('Independent shared rank budget')
            columns.append(q/norm)
        added = len(columns)-before
        if not added: break
    record = dict(schema='numivivo.org/shared-residual-oracle/v1', sourceRequestSHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  referenceEnergiesHartree=truth.tolist(), levels=levels,
                  scope='same shared finite-sector Galerkin/residual approximation; no optimized reaction or rate')
    write(out/'reference.json', record)
    return record


def global_reference(path, out):
    import numpy as np
    from scipy.linalg import eigh
    from pyscf import gto, ao2mo, fci
    from pyscf.solvent import pcm
    r = load(path, 'globalEmbedding'); molecule = r['molecule']; system = molecule['system']
    if r['gasResidualRounds'] != 0 or any(f.get('orbitalRotation') is not None for f in r['fragments']):
        raise ValueError('This independent fixture adapter requires unrotated fragment frames and zero gas enrichment')
    if molecule['basis']['representation']!='normalized-cartesian' or system['pointCharges']:
        raise ValueError('This oracle does not substitute a method for a different AO/embedding contract')
    atoms, basis = [], {}
    for i, atom in enumerate(system['nuclei']):
        label = gto.mole._atom_symbol(atom['atomicNumber'])+str(i)
        atoms.append((label,atom['positionBohr']))
        basis[label] = [[s['angularMomentum'], *[[p['exponent'],p['coefficient']] for p in s['primitives']]]
                        for s in molecule['basis']['shells'] if s['nucleusIndex']==i]
    ne = (system['alphaElectrons'],system['betaElectrons'])
    m = gto.M(atom=atoms,basis=basis,unit='Bohr',cart=True,verbose=0,
              charge=sum(a['atomicNumber'] for a in system['nuclei'])-sum(ne),spin=ne[0]-ne[1])
    overlap = m.intor('int1e_ovlp'); scale = 1/np.sqrt(np.diag(overlap))
    normalized_s = overlap*scale[:,None]*scale[None,:]
    if molecule.get('coefficients') is not None:
        c = molecule['coefficients']; c = np.array(c['values']).reshape(c['rows'],c['columns'])
    else:
        e, c = eigh(normalized_s); c = c/np.sqrt(e)
    if np.linalg.norm(c.T@normalized_s@c-np.eye(m.nao))>1e-8: raise ValueError('Nonorthonormal oracle molecular frame')
    raw_c = scale[:,None]*c
    h = raw_c.T@(m.intor('int1e_kin')+m.intor('int1e_nuc'))@raw_c
    eri = ao2mo.restore(1,ao2mo.kernel(m,raw_c),m.nao)
    gas, alpha, beta = full_matrix(h,eri,m.nao,ne,m.energy_nuc())
    all_indices = []
    for fragment in r['fragments']: all_indices += seed_indices(fragment['partition'],alpha,beta,m.nao)
    indices = sorted(set(all_indices)); redundancy = len(all_indices)-len(indices)
    if len(indices)>r['space']['maximumDimension']: raise RuntimeError('Independent union exceeds rank budget')
    solvent = pcm.PCM(m); sc=molecule['solvent']
    if sc['angularPoints']!=50 or sc['extraCavitySpheres']: raise ValueError('Unmatched oracle cavity')
    solvent.method='C-PCM';solvent.eps=sc['dielectricConstant'];solvent.lebedev_order=11
    radii=np.ones(119)
    for z,radius in sc['radiiAngstrom'].items():radii[int(z)]=radius
    solvent.radii_table=radii*sc['radiusScale']/0.529177210544
    shape=(len(alpha),len(beta))
    def solve(matrix):
        e, v=eigh(matrix[np.ix_(indices,indices)])
        vector=np.zeros(gas.shape[0]);vector[indices]=v[:,0]
        p=fci.direct_spin1.make_rdm1(vector.reshape(shape),m.nao,ne)
        return vector,p,c@p@c.T,raw_c@p@raw_c.T
    vector,p,density,density_raw=solve(gas); initial=density.copy();previous=None
    for iteration in range(1,257):
        _,old_v=solvent._get_vind(density_raw)
        potential=raw_c.T@old_v@raw_c
        biased,_,_=full_matrix(h+potential,eri,m.nao,ne,m.energy_nuc())
        vector,p,next_density,next_raw=solve(biased)
        pol,new_v=solvent._get_vind(next_raw)
        energy=float(vector@gas@vector+pol)
        error=np.linalg.norm(next_density-density)
        if previous is not None and abs(energy-previous)<1e-12 and error<1e-11: break
        density=next_density;density_raw=next_raw;previous=energy
    else: raise RuntimeError('Independent coherent fragment/PCM did not reach equilibrium')
    final_h,_,_=full_matrix(h+raw_c.T@new_v@raw_c,eri,m.nao,ne,m.energy_nuc())
    residual=final_h@vector-float(vector@final_h@vector)*vector
    projected=residual[indices]
    outside=residual.copy();outside[indices]=0
    record=dict(schema='numivivo.org/global-fragment-solvent-oracle/v1',sourceRequestSHA256=hashlib.sha256(path.read_bytes()).hexdigest(),
                variationalDimension=len(indices),fullSectorDimension=gas.shape[0],redundantSeedColumns=redundancy,
                energyHartree=energy,gasEnergyHartree=float(vector@gas@vector),polarizationEnergyHartree=float(pol),
                densityAO=next_density.tolist(),occupations=eigh(p)[0].tolist(),
                selfConsistentProjectedResidualHartree=float(np.linalg.norm(projected)),externalResidualHartree=float(np.linalg.norm(outside)),
                densityChangeFromGas=float(np.linalg.norm(next_density-initial)),iterations=iteration,
                scope='same coherent global fragment subspace and equilibrium solvent; not democratic ECC-DMET or full-CI approximation accuracy')
    write(out/'reference.json',record)
    return record


def export(directory,out):
    import pyscf
    if pyscf.__version__!='2.8.0': raise RuntimeError('Expected pinned PySCF 2.8.0')
    pyscf.lib.num_threads(1)
    residual_reference(directory/'residual-h3-631g.json',out/'residual')
    global_reference(directory/'global-h2-overlap.json',out/'h2')
    global_reference(directory/'global-lih-overlap.json',out/'lih')
    manifest={name:hashlib.sha256((out/name/'reference.json').read_bytes()).hexdigest() for name in ['residual','h2','lih']}
    write(out/'manifest.json',dict(schema='numivivo.org/scientific-closure-oracle/v1',version=pyscf.__version__,references=manifest))
    print('Independent shared-residual and overlapping-fragment solvent references generated')


def check(binary,examples,oracle,out):
    import copy
    import subprocess
    if out.exists() and any(out.iterdir()): raise RuntimeError('Use an empty output directory; previous results are not deleted')
    out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((oracle/'manifest.json').read_text())
    if manifest['schema']!='numivivo.org/scientific-closure-oracle/v1' or manifest['version']!='2.8.0' or set(manifest['references'])!={'residual','h2','lih'}:
        raise RuntimeError('Incomplete or mismatched independent reference protocol')
    checks=[]
    def require(condition,label,error=None):
        item=dict(label=label,passed=bool(condition))
        if error is not None: item['absoluteError']=error
        checks.append(item);write(out/'checks.json',dict(checks=checks,passed=all(x['passed'] for x in checks)))
        if not condition: raise RuntimeError(label)
    def compare(a,b,label,tolerance=1e-8):
        error=abs(a-b);require(math.isfinite(error) and error<=tolerance,label,error)
    def execute(request,name,result_case,fail=False):
        source=out/(name+'.request.json');destination=out/(name+'.result.json');write(source,request)
        p=subprocess.run([str(binary),'reaction-run',str(source),'--store',str(out/'store'),'--output',str(destination)],
                         capture_output=True,text=True,timeout=900)
        (out/(name+'.log')).write_text(p.stdout+'\n'+p.stderr)
        if fail:
            require(p.returncode!=0 and not destination.exists() and not Path(str(destination)+'.receipt.json').exists(),name)
            return None
        if p.returncode: raise RuntimeError(p.stderr)
        doc=json.loads(destination.read_text())
        return doc[result_case]['result'],json.loads(Path(str(destination)+'.receipt.json').read_text())
    saved={}
    for name,filename,case in [('residual','residual-h3-631g.json','residualBarrier'),('h2','global-h2-overlap.json','globalEmbedding'),('lih','global-lih-overlap.json','globalEmbedding')]:
        data=(oracle/name/'reference.json').read_bytes();reference=json.loads(data)
        if hashlib.sha256(data).hexdigest()!=manifest['references'][name] or hashlib.sha256((examples/filename).read_bytes()).hexdigest()!=reference['sourceRequestSHA256']:
            raise RuntimeError('Reference or declared input bytes changed')
        request=json.loads((examples/filename).read_text());actual,receipt=execute(request,name,case);saved[name]=(request,actual,receipt)
        require(not receipt['reused'],name+' calculation actually executes')
        repeated,again=execute(request,name+'-cached',case)
        require(again['reused'] and repeated==actual and again['result']==receipt['result'],name+' verified cache identity and complete reconstruction')
        if name=='residual':
            require(len(actual['levels'])==len(reference['levels'])==7,'all shared-residual rounds retained')
            require(actual['reducedAccuracyEstablished'] and actual['acceptedRound']==6,'two final distinct reduced ranks meet unchanged accuracy and stability criteria')
            require(actual['baseline']['assessment']=='reducedAccuracyNotEstablished','old inadequate orbital-only hierarchy is not relabeled as successful')
            require(actual['hamiltonianOperatorApplications']<=request['calculation'][case]['request']['baseline']['budget']['maximumOperatorApplications'],'physical Hamiltonian work stays bounded')
            for level,ref in zip(actual['levels'],reference['levels']):
                i=ref['round'];require(level['round']==i and level['variationalDimension']==ref['variationalDimension'] and level['fullSectorDimension']==ref['fullSectorDimension'],f'residual rank/full-sector {i}')
                require(level['stateSpaceIsReduced'] and level['variationalDimension']<level['fullSectorDimension'],f'genuine eigenproblem reduction {i}')
                require(len(level['points'])==len(ref['energiesHartree'])==9,f'all residual path points {i}')
                for j,(point,energy) in enumerate(zip(level['points'],ref['energiesHartree'])):
                    compare(point['energyHartree'],energy,f'residual point energy {i}:{j}')
                    compare(point['externalResidualHartree'],ref['externalResidualHartree'][j],f'external residual {i}:{j}',1e-7)
                    require(point['state']['orbitalCount']==6,f'all six orbital modes disclosed {i}:{j}')
                for metric in ['maximumBarrierErrorHartree','maximumProfileErrorHartree','maximumAbsoluteEnergyErrorHartree']:
                    compare(level[metric],ref[metric],f'{metric} round {i}')
                if i>0: compare(level['maximumSuccessiveChangeHartree'],ref['maximumSuccessiveChangeHartree'],f'successive profile change {i}')
            require('complete determinant-sector vectors retained' in actual['meaning'],'residual result does not claim smaller orbital memory or scalable DMRG')
        else:
            require(actual['variationalDimension']==reference['variationalDimension'] and actual['fullSectorDimension']==reference['fullSectorDimension'],name+' global projector rank')
            require(actual['redundantSeedColumns']==reference['redundantSeedColumns']==4,name+' overlapping directions are counted once')
            for metric in ['energyHartree','gasEnergyHartree','externalResidualHartree']:
                compare(actual[metric],reference[metric],name+' '+metric,2e-8)
            compare(actual['equilibriumField']['polarizationEnergyHartree'],reference['polarizationEnergyHartree'],name+' polarization',2e-8)
            density=actual['densityAO'];expected=reference['densityAO'];flat=[x for row in expected for x in row]
            require(density['rows']==len(expected) and density['columns']==len(expected[0]) and len(density['values'])==len(flat),name+' density matrix cardinality')
            error=math.sqrt(sum((a-b)**2 for a,b in zip(density['values'],flat)))
            require(math.isfinite(error) and error<=2e-7,name+' independent global density',error)
            require(len(actual['occupations'])==len(reference['occupations']),name+' occupation cardinality')
            for i,(a,b) in enumerate(zip(actual['occupations'],reference['occupations'])): compare(a,b,name+f' occupation {i}',2e-7)
            require(actual['selfConsistentProjectedResidualHartree']<=request['calculation'][case]['request']['stationarityToleranceHartree'],name+' actual returned-field stationarity')
            if name=='lih':
                require(actual['externalResidualHartree']>1e-4 and actual['externalResidualHartree']>actual['selfConsistentProjectedResidualHartree'],name+' omitted correlation is not hidden in a projected residual')
            else:
                require(actual['externalResidualHartree']<1e-8,name+' symmetric reduced subspace contains the full ground state')
            if name=='lih': require(reference['densityChangeFromGas']>1e-4,'polar global density has nontrivial equilibrium solvent feedback')
            require('not democratic ECC-DMET' in actual['method'],'global closure keeps its distinct energy-functional meaning '+name)
    residual=saved['residual'][0]
    bad=copy.deepcopy(residual);bad['calculation']['residualBarrier']['request']['space']['maximumDimension']=10
    execute(bad,'rank-budget-fails','residualBarrier',True)
    global_input=saved['h2'][0]
    bad=copy.deepcopy(global_input);bad['calculation']['globalEmbedding']['request']['fragments'][0]['partition']['doublyOccupiedCore']=[0]
    execute(bad,'incompatible-fragment-core-fails','globalEmbedding',True)
    # Negative preflight reports are valid data. They must never masquerade as
    # successful paper execution or authenticated author-provided structures.
    for template,count in [('paper-michael-inputs',5),('paper-btk-inputs',8)]:
        pth=out/(template+'.template.json')
        p=subprocess.run([str(binary),'reaction-template',template,'--output',str(pth)],capture_output=True,text=True,timeout=60)
        if p.returncode: raise RuntimeError(p.stderr)
        data=json.loads(pth.read_text());report,_=execute(data,template,'reproductionPreflight')
        require(len(report['missingRoles'])==count and not report['suppliedInputPackageConsistent'] and not report['reproductionExecuted'] and not report['sourceAuthenticityVerified'],template+' missing assets cannot yield false qualification')
        bad=copy.deepcopy(data);bad['calculation']['reproductionPreflight']['package']['assets']=[dict(role='precomplex',origin='authorSupplied',sourceIdentifier='unverified',sha256='0'*64,payload='e30=')]
        execute(bad,template+'-corrupt-payload','reproductionPreflight',True)
    # Common reaction routing must recognize connectivity even when validation
    # rejects a malformed request. A result from another schema is not accepted.
    digest=bytes(saved['residual'][2]['result']['bytes']).hex();object_path=out/'store'/'objects'/'sha256'/digest[:2]/digest[2:4]/digest
    original=object_path.read_bytes()
    if hashlib.sha256(original).hexdigest()!=digest: raise RuntimeError('Unknown artifact layout')
    try:
        object_path.write_bytes(original+b'\n');execute(residual,'corrupt-residual-cache','residualBarrier',True)
    finally: object_path.write_bytes(original)
    write(out/'checks.json',dict(schema='numivivo.org/scientific-closure-checks/v1',checks=checks,passed=True,
        executableSHA256=hashlib.sha256(binary.read_bytes()).hexdigest(),referenceManifestSHA256=hashlib.sha256((oracle/'manifest.json').read_bytes()).hexdigest(),
        scope='bounded shared residual eigenproblem and coherent overlapping Fock projectors with correlated solvent; matched inputs, no paper reproduction or kinetic-rate certificate'))
    print('PASS',len(checks),'independent scientific-closure and production CLI checks')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('export');p.add_argument('examples',type=Path);p.add_argument('out',type=Path)
    p=sub.add_parser('check');p.add_argument('binary',type=Path);p.add_argument('examples',type=Path);p.add_argument('reference',type=Path);p.add_argument('out',type=Path)
    args=parser.parse_args()
    if args.command=='export':export(args.examples,args.out)
    else:check(args.binary.resolve(),args.examples,args.reference,args.out)
