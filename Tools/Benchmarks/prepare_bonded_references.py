#!/usr/bin/env python3
"""Independent isolated bonded-force regression references with gradient checks."""
import argparse
import json
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit as u
from prepare_references import export_system,fp,sha,vec,write


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--template",type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    config=json.loads(args.template.read_text())["configuration"]
    config.update(electrostatics="cutoff",neighborListEnabled=False,cutoffNM=10)
    rows=[]
    xyz=np.array([[0,0.25,0.125],[0.25,0.125,0],[0.5,0.25,0.125],[0.75,0.125,0.375]],dtype=float)
    for name in ["bond","angle","torsion-zero","torsion-phase","torsion-pi"]:
        root=args.out/name;root.mkdir();system=mm.System();nb=mm.NonbondedForce()
        for _ in range(4):system.addParticle(12);nb.addParticle(0,0,0)
        system.addForce(nb)
        if name=="bond":force=mm.HarmonicBondForce();force.addBond(0,1,0.23,100)
        elif name=="angle":force=mm.HarmonicAngleForce();force.addAngle(0,1,2,1.3,50)
        else:
            force=mm.PeriodicTorsionForce()
            phase={"torsion-zero":0,"torsion-phase":0.37,"torsion-pi":np.pi}[name]
            force.addTorsion(0,1,2,3,3,phase,2.2)
        system.addForce(force)
        xml=mm.XmlSerializer.serialize(system).encode();(root/"reference-system.xml").write_bytes(xml)
        native=export_system(system,name,xml)
        integrator=mm.VerletIntegrator(0.001)
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
        context.setPositions(xyz)
        state=context.getState(getEnergy=True,getForces=True)
        energy=state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        forces=state.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole/u.nanometer)
        fd=np.zeros_like(xyz)
        for i in range(4):
            for axis in range(3):
                plus=xyz.copy();minus=xyz.copy();plus[i,axis]+=1e-5;minus[i,axis]-=1e-5
                context.setPositions(plus);ep=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                context.setPositions(minus);em=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
                fd[i,axis]=-(ep-em)/(2e-5)
        error=float(np.max(np.abs(fd-forces)))
        assert error<1e-5,(name,error)
        request=dict(schema="numivivo.org/md-benchmark-request/v1",identifier=name,system=native,configuration=config,
            references=[dict(identifier="source",geometry=dict(particlePositionsNM=list(map(vec,xyz))),energyKJPerMol=energy,forcesKJPerMolNM=list(map(vec,forces)))],
            referenceProvenance=dict(engine="OpenMM",version=mm.__version__,platform="Reference",serializedSystemSHA256=sha(xml)),
            limits=dict(energyAbsolutePerParticleKJPerMol=0.002,forceNormalizedRMS=0.001,forceNormalizedMaximum=0.01,forceNormalizationFloor=1),dynamicsSteps=0)
        write(root/"request.json",request)
        write(root/"finite-difference.json",dict(maximumErrorKJPerMolNM=error,stepNM=1e-5,forces=list(map(vec,fd))))
        rows.append(dict(identifier=name,particles=4,requestSHA256=sha((root/"request.json").read_bytes()),serializedSystemSHA256=sha(xml),status="prepared"))
    write(args.out/"manifest.json",dict(schema="numivivo.org/md-reference-campaign/v1",cases=rows,generatorSHA256=sha(Path(__file__).read_bytes())))


if __name__=="__main__":main()
