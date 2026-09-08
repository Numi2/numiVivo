#!/usr/bin/env python3
"""Re-evaluate retained MD endpoints with independent FP64 energy accounting.

OpenMM Reference is an offline auditor only. It never advances or forces the
native trajectory. Original preparation/execution failures remain visible.
"""
import argparse
import json
import math
from pathlib import Path
import platform
from importlib import metadata
from run_campaign import digest, write, verify_checkpoint


def read(path):return json.loads(path.read_text())


def kinetic_check(system,checkpoint,observation):
    particles=system["particles"];velocities=checkpoint["velocitiesNMPerPS"]
    if len(particles)!=len(velocities) or not particles:raise ValueError("kinetic particle count mismatch")
    terms=[]
    for i,(particle,velocity) in enumerate(zip(particles,velocities)):
        mass=particle["massDa"]
        if particle["index"]!=i or particle["role"]!="atom" or not math.isfinite(mass) or mass<=0:
            raise ValueError("endpoint kinetic audit requires positive-mass physical atoms")
        if not all(math.isfinite(velocity[k]) for k in ("x","y","z")):raise ValueError("nonfinite velocity")
        terms.append(.5*mass*math.fsum(velocity[k]**2 for k in ("x","y","z")))
    expected=math.fsum(terms);actual=observation["kineticEnergyKJPerMol"]
    if not math.isfinite(actual) or actual<0:raise ValueError("invalid native kinetic energy")
    # Nonnegative pairwise FP32 accumulation, plus eight rounded operations
    # for mass conversion and each particle's squared-speed/energy term.
    operations=math.ceil(math.log2(len(particles)))+8;unit_roundoff=2**-24
    gamma=operations*unit_roundoff/(1-operations*unit_roundoff)
    bound=gamma*expected+1e-8
    return dict(referenceKineticEnergyKJPerMol=expected,nativeKineticEnergyKJPerMol=actual,
        absoluteErrorKJPerMol=abs(actual-expected),roundoffBoundKJPerMol=bound,
        roundoffOperationBound=operations,passed=abs(actual-expected)<=bound)


def reference_potential(xml,checkpoint):
    import openmm as mm
    from openmm import unit
    system=mm.XmlSerializer.deserialize(xml);integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
    cell=checkpoint.get("periodicCell")
    if cell:context.setPeriodicBoxVectors(*[mm.Vec3(*(cell[k][axis] for axis in ("x","y","z")))*unit.nanometer for k in ("a","b","c")])
    context.setPositions([mm.Vec3(*(v[k] for k in ("x","y","z"))) for v in checkpoint["positionsNM"]]*unit.nanometer)
    result=float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    del context,integrator
    if not math.isfinite(result):raise ValueError("nonfinite independent potential energy")
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("references","campaign","out"):parser.add_argument("--"+name,type=Path,required=True)
    a=parser.parse_args();manifest=read(a.references/"manifest.json");scorecard=read(a.campaign/"scorecard.json")
    if scorecard["referenceManifestSHA256"]!=digest(a.references/"manifest.json"):raise ValueError("changed reference manifest")
    if a.out.exists():raise FileExistsError(a.out)
    tool_hash=digest(Path(__file__));rows=[]
    for case in manifest["cases"]:
        name=case["identifier"]
        if not name or not all(c.isalnum() or c in "-_" for c in name):raise ValueError("unsafe case identifier")
        row=dict(identifier=name,outcome="failed")
        if case["status"]!="prepared":
            row.update(outcome="preparation-failed",error=case.get("error"));rows.append(row);continue
        try:
            native=next(r for r in scorecard["cases"] if r["identifier"]==name)
            if native["outcome"]!="passed":raise ValueError("native execution did not pass: "+native["outcome"])
            request_path=a.references/name/"request.json";report_path=a.campaign/name/"report.json";xml_path=a.references/name/"reference-system.xml"
            if digest(request_path)!=case["requestSHA256"] or digest(request_path)!=native["requestSHA256"] or digest(report_path)!=native["reportSHA256"]:
                raise ValueError("changed native request/report")
            if digest(xml_path)!=case["serializedSystemSHA256"]:raise ValueError("changed independent reference system")
            request=read(request_path);report=read(report_path);dynamics=report["dynamics"]
            if report["numericalContract"]!="numivivo.org/md-metal-numerics/v6" or report["outcome"]!="passed":
                raise ValueError("endpoint roundoff audit requires a passed native v6 report")
            kinetics=[]
            for key,observation in (("preparedCheckpoint","start"),("finalCheckpoint","end")):
                checkpoint=dynamics[key];sample=dynamics[observation]
                if checkpoint["acceptedStep"]!=sample["stepIndex"] or checkpoint["timePS"]!=sample["timePS"]:raise ValueError("endpoint clock mismatch")
                if not verify_checkpoint(request,checkpoint)["passed"]:raise ValueError("independent endpoint constraints failed")
                kinetics.append(dict(endpoint=observation,**kinetic_check(request["system"],checkpoint,sample)))
            reference=reference_potential(xml_path.read_text(),dynamics["finalCheckpoint"])
            native_potential=dynamics["end"]["potentialEnergyKJPerMol"]
            if not math.isfinite(native_potential):raise ValueError("nonfinite native potential energy")
            error=abs(native_potential-reference)/len(request["system"]["particles"])
            limit=request["limits"]["energyAbsolutePerParticleKJPerMol"]
            potential=dict(referenceEnergyKJPerMol=reference,nativeEnergyKJPerMol=native_potential,
                errorPerParticleKJPerMol=error,limitPerParticleKJPerMol=limit,passed=math.isfinite(error) and error<=limit)
            row.update(outcome="passed" if potential["passed"] and all(k["passed"] for k in kinetics) else "failed",
                kinetic=kinetics,finalPotential=potential,requestSHA256=digest(request_path),reportSHA256=digest(report_path),serializedSystemSHA256=digest(xml_path))
        except Exception as error:row["error"]=f"{type(error).__name__}: {error}"
        rows.append(row);print(json.dumps(row,allow_nan=False),flush=True)
    if digest(Path(__file__))!=tool_hash:raise ValueError("endpoint auditor changed during execution")
    prepared=[r for r in rows if r["outcome"]!="preparation-failed"]
    result=dict(schema="numivivo.org/md-endpoint-energy-audit/v1",auditorSHA256=tool_hash,
        scorecardSHA256=digest(a.campaign/"scorecard.json"),referenceManifestSHA256=digest(a.references/"manifest.json"),
        environment=dict(python=platform.python_version(),openmm=metadata.version("openmm")),cases=rows,
        preparedPassed=bool(prepared) and all(r["outcome"]=="passed" for r in prepared),passed=bool(rows) and all(r["outcome"]=="passed" for r in rows),
        scope="Offline FP64 final potential and start/end kinetic agreement, using saved exact checkpoint words. No independent trajectory, ensemble or experimental claim.")
    write(a.out,result);return 0 if result["passed"] else 1


if __name__=="__main__":raise SystemExit(main())
