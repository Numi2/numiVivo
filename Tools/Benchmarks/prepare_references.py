#!/usr/bin/env python3
"""Generate independent OpenMM observations; never executes NumiVivo physics.

The exact OpenMM System XML, source files, rounded coordinates, reference forces,
and package versions are retained. This exporter rejects unrepresented forces.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import urllib.request

import numpy as np
import scipy
import openmm as mm
from openmm import app, unit as u

REVISION = "f6ef22a8b9f66e582df2ffa62f3bb6516de43536"
BASE = f"https://raw.githubusercontent.com/choderalab/openmmtools/{REVISION}/openmmtools/data/"
CASES = {
    "water": ("waterbox/watbox216.prmtop", "waterbox/watbox216.crd"),
    "alanine": ("alanine-dipeptide-explicit/alanine-dipeptide.prmtop", "alanine-dipeptide-explicit/alanine-dipeptide.crd"),
    "protein": ("dhfr/JAC.prmtop", "dhfr/JAC.inpcrd"),
    "complex": ("T4-lysozyme-L99A-implicit/complex.prmtop", "T4-lysozyme-L99A-implicit/complex.crd"),
    "dna": ("dna_dodecamer_explicit/minimized_dna_dodecamer.pdb",),
    "membrane": (),
    "ions": (),
    "water-orthogonal": (),
}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    with path.open("x") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")


def vec(value):
    return dict(zip(("x", "y", "z"), map(float, value)))


def fp(data):
    return {"bytes": list(hashlib.sha256(data).digest())}


def scalar(value, units):
    return float(value.value_in_unit(units))


def load_case(name, root):
    sources = []
    files = []
    for relative in CASES[name]:
        data = urllib.request.urlopen(BASE + relative, timeout=60).read()
        path = root / Path(relative).name
        path.write_bytes(data)
        sources.append({"url": BASE + relative, "sha256": sha(data), "file": path.name})
        files.append(path)
    kwargs = dict(constraints=app.HBonds, rigidWater=True, removeCMMotion=False)
    if name in ("water", "alanine", "protein", "complex"):
        top = app.AmberPrmtopFile(str(files[0]))
        state = app.AmberInpcrdFile(str(files[1]))
        periodic = name != "complex"
        system = top.createSystem(nonbondedMethod=app.PME if periodic else app.NoCutoff,
                                  nonbondedCutoff=0.8*u.nanometer, **kwargs)
        if periodic:
            system.setDefaultPeriodicBoxVectors(*state.boxVectors)
        positions = state.positions
    else:
        forcefields = ["amber14-all.xml", "amber14/tip3p.xml"]
        ff = app.ForceField(*forcefields)
        if name == "dna":
            pdb = app.PDBFile(str(files[0]))
            model = app.Modeller(pdb.topology, pdb.positions)
        elif name == "membrane":
            source = Path(app.__file__).parent / "data" / "POPC.pdb"
            data = source.read_bytes(); (root/source.name).write_bytes(data)
            sources.append({"package": "openmm", "file": source.name, "sha256": sha(data)})
            pdb = app.PDBFile(str(source))
            model = app.Modeller(pdb.topology, pdb.positions)
        else:
            # Deterministic prepackaged water; replace two widely separated
            # water residues with one sodium and one chloride below.
            source = Path(app.__file__).parent / "data" / "tip3p.pdb"
            data = source.read_bytes(); (root/source.name).write_bytes(data)
            sources.append({"package": "openmm", "file": source.name, "sha256": sha(data)})
            pdb = app.PDBFile(str(source))
            model = app.Modeller(pdb.topology, pdb.positions)
            if name == "ions":
                waters=list(model.topology.residues())
                chosen=[waters[0],waters[len(waters)//2]]
                coordinates=[model.positions[next(r.atoms()).index] for r in chosen]
                model.delete(chosen)
                top=app.Topology(); chain=top.addChain()
                for symbol in ["Na","Cl"]:
                    residue=top.addResidue(symbol.upper(),chain)
                    top.addAtom(symbol,app.Element.getBySymbol(symbol),residue)
                model.add(top,u.Quantity([v.value_in_unit(u.nanometer) for v in coordinates],u.nanometer))
        system=ff.createSystem(model.topology,nonbondedMethod=app.PME,nonbondedCutoff=0.8*u.nanometer,**kwargs)
        positions=model.positions; periodic=True
        sources.append({"forcefields": forcefields, "package": "openmm", "version": mm.__version__})
    return system,positions,periodic,sources


def export_system(system, name, xml):
    supported=(mm.NonbondedForce,mm.HarmonicBondForce,mm.HarmonicAngleForce,mm.PeriodicTorsionForce)
    for force in system.getForces():
        if not isinstance(force,supported):
            raise ValueError(f"unrepresented force {type(force).__name__}")
    nonbonded=[f for f in system.getForces() if isinstance(f,mm.NonbondedForce)]
    if len(nonbonded)!=1 or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
        raise ValueError("requires one nonbonded force and physical particles")
    nb=nonbonded[0]
    if nb.getNumGlobalParameters() or nb.getNumParticleParameterOffsets() or nb.getNumExceptionParameterOffsets():
        raise ValueError("unrepresented nonbonded parameter offsets")
    result=dict(schema="numivivo.org/classical-system/v1",identifier=name,structureFingerprint=fp(xml),
                parameterSourceFingerprints=[fp(xml)],mixingRule="lorentzBerthelot",particles=[],
                bonds=[],angles=[],torsions=[],constraints=[],nonbondedExceptions=[],metadata={"reference":"OpenMM serialized System"})
    types={}
    charges=[]
    for i in range(system.getNumParticles()):
        q,sigma,epsilon=nb.getParticleParameters(i)
        q=scalar(q,u.elementary_charge);sigma=scalar(sigma,u.nanometer);epsilon=scalar(epsilon,u.kilojoule_per_mole)
        key=(sigma,epsilon);types.setdefault(key,f"lj-{len(types)}");charges.append(q)
        result["particles"].append(dict(index=i,atomIndex=i,typeIdentifier=types[key],role="atom",
            massDa=scalar(system.getParticleMass(i),u.dalton),chargeE=q,sigmaNM=sigma,epsilonKJPerMol=epsilon))
    for force in system.getForces():
        if isinstance(force,mm.HarmonicBondForce):
            for i in range(force.getNumBonds()):
                a,b,length,k=force.getBondParameters(i)
                result["bonds"].append(dict(a=a,b=b,lengthNM=scalar(length,u.nanometer),forceConstant=scalar(k,u.kilojoule_per_mole/u.nanometer**2)))
        elif isinstance(force,mm.HarmonicAngleForce):
            for i in range(force.getNumAngles()):
                a,b,c,theta,k=force.getAngleParameters(i)
                result["angles"].append(dict(a=a,b=b,c=c,angleRadians=scalar(theta,u.radian),forceConstant=scalar(k,u.kilojoule_per_mole/u.radian**2)))
        elif isinstance(force,mm.PeriodicTorsionForce):
            for i in range(force.getNumTorsions()):
                a,b,c,d,n,phase,k=force.getTorsionParameters(i)
                result["torsions"].append(dict(a=a,b=b,c=c,d=d,periodicity=n,phaseRadians=scalar(phase,u.radian),barrierKJPerMol=scalar(k,u.kilojoule_per_mole),improper=False))
    for i in range(system.getNumConstraints()):
        a,b,d=system.getConstraintParameters(i)
        result["constraints"].append(dict(a=a,b=b,distanceNM=scalar(d,u.nanometer)))
    for i in range(nb.getNumExceptions()):
        a,b,q,sigma,epsilon=nb.getExceptionParameters(i)
        product=charges[a]*charges[b];q=scalar(q,u.elementary_charge**2)
        if product==0 and q!=0:raise ValueError("unrepresented independent exception charge")
        result["nonbondedExceptions"].append(dict(a=a,b=b,coulombScale=q/product if product else 0,
            lennardJonesScale=1,sigmaOverrideNM=scalar(sigma,u.nanometer),epsilonOverrideKJPerMol=scalar(epsilon,u.kilojoule_per_mole)))
    return result


def prepare(name, root, steps, project_constraints=False):
    root.mkdir()
    system,positions,periodic,sources=load_case(name,root)
    nb=next(f for f in system.getForces() if isinstance(f,mm.NonbondedForce))
    nb.setUseDispersionCorrection(False)
    nb.setUseSwitchingFunction(False)
    if periodic:nb.setEwaldErrorTolerance(1e-7)
    xml=mm.XmlSerializer.serialize(system).encode();(root/"reference-system.xml").write_bytes(xml)
    native=export_system(system,name,xml)
    box=np.asarray([v.value_in_unit(u.nanometer) for v in system.getDefaultPeriodicBoxVectors()],dtype=np.float32).astype(float)
    # The comparison uses the exact native FP32 cell in both engines.
    if periodic:system.setDefaultPeriodicBoxVectors(*[mm.Vec3(*v)*u.nanometer for v in box])
    xml=mm.XmlSerializer.serialize(system).encode();(root/"reference-system.xml").write_bytes(xml)
    native["structureFingerprint"]=fp(xml);native["parameterSourceFingerprints"]=[fp(xml)]
    integrator=mm.VerletIntegrator(0.001*u.picosecond)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
    xyz=np.asarray(positions.value_in_unit(u.nanometer),dtype=np.float32).astype(float)
    cutoff=0.8 if periodic else float(max(1,np.linalg.norm(np.ptp(xyz,axis=0))+1))
    config=dict(schema="numivivo.org/md-configuration/v1",timeStepPS=0.001,cutoffNM=cutoff,
        neighborSkinNM=0.1,electrostatics="pme" if periodic else "cutoff",relativeDielectric=1,
        reactionFieldDielectric=78.3,pmeTolerance=1e-7,pmeGridSpacingNM=0.06,ensemble="nvt",
        thermostat="langevinMiddle",targetTemperatureK=300,frictionPerPS=1,barostat="none",
        barostatInterval=25,barostatMaximumLogVolumeStep=0.01,constraintTolerance=1e-6,
        maximumConstraintIterations=128,neighborRebuildInterval=10,neighborListEnabled=periodic,
        maximumNeighborsPerParticle=768,randomSeed=1729)
    cell=dict(zip(("a","b","c"),map(vec,box))) if periodic else None
    references=[]
    # Same coordinates, shifted coordinates, and a deterministic small distortion
    # expose force terms and image handling without running an equilibration protocol.
    for label,coords in [("source",xyz),("translated",xyz+np.array([0.125,-0.0625,0.03125])),
                         ("perturbed",xyz+0.0001*np.sin(np.arange(xyz.size).reshape(xyz.shape)))]:
        coords=coords.astype(np.float32).astype(float)
        context.setPositions(coords*u.nanometer)
        state=context.getState(getEnergy=True,getForces=True)
        references.append(dict(identifier=label,geometry=dict(particlePositionsNM=list(map(vec,coords)),periodicCell=cell),
            energyKJPerMol=scalar(state.getPotentialEnergy(),u.kilojoule_per_mole),
            forcesKJPerMolNM=list(map(vec,state.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole/u.nanometer)))))
    request=dict(schema="numivivo.org/md-benchmark-request/v1",identifier=name,system=native,configuration=config,
        references=references,referenceProvenance=dict(engine="OpenMM",version=mm.__version__,platform="Reference",
        serializedSystemSHA256=sha(xml),sourceRevision=REVISION,precision="FP64 reference at exact FP32 coordinates and cell",
        scope="static Hamiltonian comparison; prepared parameters; no experimental or ensemble qualification"),
        limits=dict(energyAbsolutePerParticleKJPerMol=0.002,forceNormalizedRMS=0.001,forceNormalizedMaximum=0.01,forceNormalizationFloor=1),
        dynamicsSteps=steps)
    if project_constraints:request["dynamicsPreparation"]="projectConstraints"
    write(root/"request.json",request)
    write(root/"sources.json",dict(sources=sources,python=platform.python_version(),openmm=mm.__version__,numpy=np.__version__,scipy=scipy.__version__,
        preparation="original supplied geometry; no minimization; no dispersion correction; no COM removal"))
    return dict(identifier=name,particles=system.getNumParticles(),requestSHA256=sha((root/"request.json").read_bytes()),
                serializedSystemSHA256=sha(xml),status="prepared")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--cases",nargs="+",choices=list(CASES),default=list(CASES))
    parser.add_argument("--steps",type=int,default=100)
    parser.add_argument("--project-initial-constraints",action="store_true",help="explicitly project imported state before native dynamics; static references remain unchanged")
    args=parser.parse_args()
    if not 0<=args.steps<=100_000:parser.error("steps must be 0...100000")
    args.out.mkdir(parents=True,exist_ok=False)
    records=[]
    for name in args.cases:
        try:record=prepare(name,args.out/name,args.steps,args.project_initial_constraints)
        except Exception as e:record=dict(identifier=name,status="preparation-failed",error=f"{type(e).__name__}: {e}")
        records.append(record);print(json.dumps(record),flush=True)
    write(args.out/"manifest.json",dict(schema="numivivo.org/md-reference-campaign/v1",cases=records,
        generatorSHA256=sha(Path(__file__).read_bytes()),sourceRevision=REVISION))
    return 0 if all(r["status"]=="prepared" for r in records) else 1


if __name__=="__main__":raise SystemExit(main())
