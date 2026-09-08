#!/usr/bin/env python3
"""Re-evaluate retained MD endpoints with independent FP64 energy accounting.

OpenMM Reference is an offline auditor only. It never advances or forces the
native trajectory. Original preparation/execution failures remain visible.
"""
import argparse
import json
import hashlib
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


def whole_molecule_positions(system,checkpoint):
    """Select whole finite molecular images without changing minimum-image geometry.

    Native v6 wraps atoms independently and evaluates bonded displacements by
    minimum image. OpenMM's serialized bonded forces/exceptions require whole
    molecules. Reject winding cycles, ambiguous edges and long exception images.
    See https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#periodic.
    """
    xyz=[[v[k] for k in ("x","y","z")] for v in checkpoint["positionsNM"]];n=len(xyz)
    if n!=len(system["particles"]) or not all(math.isfinite(v) for p in xyz for v in p):raise ValueError("invalid coordinate count or value")
    cell=checkpoint.get("periodicCell");shifts=[[0,0,0] for _ in xyz];components=n
    if cell:
        matrix=[[cell[k][a] for a in ("x","y","z")] for k in ("a","b","c")]
        if any(not math.isfinite(v) or (i==j and v<=0) or (i!=j and v!=0) for i,row in enumerate(matrix) for j,v in enumerate(row)):
            raise ValueError("endpoint image audit currently requires an orthogonal positive cell")
        lengths=[matrix[i][i] for i in range(3)];adj=[set() for _ in xyz]
        for key,indices in (("bonds","ab"),("constraints","ab"),("angles","abc"),("torsions","abcd")):
            for term in system.get(key,[]):
                for ka,kb in zip(indices,indices[1:]):
                    a,b=term[ka],term[kb]
                    if type(a) is not int or type(b) is not int or not 0<=a<n or not 0<=b<n or a==b:raise ValueError("invalid molecular edge")
                    adj[a].add(b);adj[b].add(a)
        shifts=[None]*n;groups=[None]*n;components=0
        for root in range(n):
            if shifts[root] is not None:continue
            shifts[root]=[0,0,0];groups[root]=components;components+=1;stack=[root]
            while stack:
                a=stack.pop()
                for b in sorted(adj[a]):
                    fractional=[(xyz[b][i]-xyz[a][i])/lengths[i] for i in range(3)]
                    if any(abs(abs(x-round(x))-.5)<=1e-10 for x in fractional):raise ValueError("ambiguous half-box molecular edge")
                    image=[shifts[a][i]-round(fractional[i]) for i in range(3)]
                    if shifts[b] is None:shifts[b]=image;groups[b]=groups[a];stack.append(b)
                    elif shifts[b]!=image:raise ValueError("molecular cycle winds through the periodic cell")
        for term in system.get("nonbondedExceptions",[]):
            a,b=term["a"],term["b"]
            if groups[a]!=groups[b]:raise ValueError("exception joins disconnected molecular components")
            for i in range(3):
                required=-round((xyz[b][i]-xyz[a][i])/lengths[i])
                if shifts[b][i]-shifts[a][i]!=required:raise ValueError("exception image differs from native minimum image")
        xyz=[[v+shift[i]*lengths[i] for i,v in enumerate(position)] for position,shift in zip(xyz,shifts)]
    digest_value=lambda value:hashlib.sha256(json.dumps(value,separators=(",",":"),allow_nan=False).encode()).hexdigest()
    return xyz,dict(convention="whole finite molecular images; integer cell translations only" if cell else "unchanged vacuum coordinates",
        components=components,shiftedParticles=sum(any(s) for s in shifts),imageShiftsSHA256=digest_value(shifts),positionsSHA256=digest_value(xyz))


def reference_potential(xml,checkpoint,native_system):
    import openmm as mm
    from openmm import unit
    system=mm.XmlSerializer.deserialize(xml);integrator=mm.VerletIntegrator(.001)
    context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
    cell=checkpoint.get("periodicCell")
    if cell:context.setPeriodicBoxVectors(*[mm.Vec3(*(cell[k][axis] for axis in ("x","y","z")))*unit.nanometer for k in ("a","b","c")])
    positions,images=whole_molecule_positions(native_system,checkpoint)
    context.setPositions([mm.Vec3(*v) for v in positions]*unit.nanometer)
    result=float(context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
    del context,integrator
    if not math.isfinite(result):raise ValueError("nonfinite independent potential energy")
    return result,images


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
            reference,images=reference_potential(xml_path.read_text(),dynamics["finalCheckpoint"],request["system"])
            native_potential=dynamics["end"]["potentialEnergyKJPerMol"]
            if not math.isfinite(native_potential):raise ValueError("nonfinite native potential energy")
            error=abs(native_potential-reference)/len(request["system"]["particles"])
            limit=request["limits"]["energyAbsolutePerParticleKJPerMol"]
            potential=dict(referenceEnergyKJPerMol=reference,nativeEnergyKJPerMol=native_potential,
                errorPerParticleKJPerMol=error,limitPerParticleKJPerMol=limit,referenceImages=images,passed=math.isfinite(error) and error<=limit)
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
