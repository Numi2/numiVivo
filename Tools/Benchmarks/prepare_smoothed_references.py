#!/usr/bin/env python3
"""Create a separately identified smooth-LJ reference panel for NVE refinement.

The original sharp-cutoff panel is retained. All particle/bonded/constraint
parameters and supplied geometries are preserved; both reference and native
models explicitly switch periodic LJ interactions from 0.7 to 0.8 nm. New
OpenMM Reference energies and forces are calculated for the changed model.
"""
import argparse
import copy
import json
from pathlib import Path
import shutil
import openmm as mm
from openmm import unit as u
from prepare_references import export_system,scalar,sha,vec,write


def prepare(case,source,out):
    root=source/case["identifier"];original=json.loads((root/"request.json").read_text())
    assert sha((root/"request.json").read_bytes())==case["requestSHA256"]
    xml=(root/"reference-system.xml").read_bytes();assert sha(xml)==case["serializedSystemSHA256"]
    system=mm.XmlSerializer.deserialize(xml.decode())
    forces=[f for f in system.getForces() if isinstance(f,mm.NonbondedForce)]
    assert len(forces)==1;nb=forces[0];periodic=nb.getNonbondedMethod()!=mm.NonbondedForce.NoCutoff
    identifier=case["identifier"]+ ("-smooth" if periodic else "-vacuum")
    directory=out/identifier;shutil.copytree(root,directory)
    request=copy.deepcopy(original);request["identifier"]=identifier
    if periodic:
        assert nb.getNonbondedMethod()==mm.NonbondedForce.PME and abs(scalar(nb.getCutoffDistance(),u.nanometer)-0.8)<1e-12
        assert not nb.getUseSwitchingFunction() and not nb.getUseDispersionCorrection()
        nb.setUseSwitchingFunction(True);nb.setSwitchingDistance(0.7*u.nanometer)
        request["configuration"]["lennardJonesSwitchOnNM"]=0.7
    xml=mm.XmlSerializer.serialize(system).encode();native=export_system(system,identifier,xml)
    excluded={"identifier","structureFingerprint","parameterSourceFingerprints","metadata"}
    assert {k:v for k,v in native.items() if k not in excluded}=={k:v for k,v in original["system"].items() if k not in excluded},"unexpected physical parameter change"
    request["system"]=native
    integrator=mm.VerletIntegrator(0.001*u.picosecond)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
    for reference in request["references"]:
        positions=[mm.Vec3(p["x"],p["y"],p["z"]) for p in reference["geometry"]["particlePositionsNM"]]
        context.setPositions(positions*u.nanometer)
        state=context.getState(getEnergy=True,getForces=True)
        reference["energyKJPerMol"]=scalar(state.getPotentialEnergy(),u.kilojoule_per_mole)
        reference["forcesKJPerMolNM"]=list(map(vec,state.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole/u.nanometer)))
    request["referenceProvenance"].update(serializedSystemSHA256=sha(xml),parentRequestSHA256=case["requestSHA256"],
        modelVariant="periodic LJ switched 0.7 to 0.8 nm" if periodic else "unchanged vacuum NoCutoff model",
        version=mm.__version__)
    (directory/"reference-system.xml").write_bytes(xml);write(directory/"request.json",request)
    write(directory/"variant.json",dict(parentIdentifier=case["identifier"],parentRequestSHA256=case["requestSHA256"],
        parentReferenceSystemSHA256=case["serializedSystemSHA256"],scope=request["referenceProvenance"]["modelVariant"]))
    return dict(identifier=identifier,particles=system.getNumParticles(),status="prepared",
        requestSHA256=sha((directory/"request.json").read_bytes()),serializedSystemSHA256=sha(xml),parentIdentifier=case["identifier"])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--references",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    source=a.references.resolve(strict=True);parent=json.loads((source/"manifest.json").read_text());rows=[]
    for case in parent["cases"]:
        if case["status"]!="prepared":rows.append(case);continue
        try:row=prepare(case,source,a.out)
        except Exception as e:row=dict(identifier=case["identifier"],status="preparation-failed",error=f"{type(e).__name__}: {e}")
        rows.append(row);print(json.dumps(row),flush=True)
    write(a.out/"manifest.json",dict(schema=parent["schema"],cases=rows,generatorSHA256=sha(Path(__file__).read_bytes()),
        parentManifestSHA256=sha((source/"manifest.json").read_bytes()),openmm=mm.__version__,
        scope="separate smooth-LJ model; original sharp-cutoff panel and failures retained"))
    return 0 if all(r["status"]=="prepared" for r in rows) else 1


if __name__=="__main__":raise SystemExit(main())
