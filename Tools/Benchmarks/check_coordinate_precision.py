#!/usr/bin/env python3
"""Measure constraint error introduced solely by FP32 coordinate rounding.

This is a reference diagnostic, not a proof that no representable constrained
geometry exists. It does not change native tolerances or trajectories.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import openmm as mm
from openmm import unit as u
from prepare_references import sha,write


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--cases",nargs="+",default=["water-orthogonal","protein"])
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    records=[]
    for name in args.cases:
        assert name and all(c.isalnum() or c in "-_" for c in name)
        root=args.references/name;raw=(root/"reference-system.xml").read_bytes()
        request=json.loads((root/"request.json").read_text())
        assert sha(raw)==request["referenceProvenance"]["serializedSystemSHA256"]
        system=mm.XmlSerializer.deserialize(raw.decode());integrator=mm.VerletIntegrator(0.001)
        context=mm.Context(system,integrator,mm.Platform.getPlatformByName("Reference"))
        xyz=np.array([[v[k] for k in ["x","y","z"]] for v in request["references"][0]["geometry"]["particlePositionsNM"]])
        context.setPositions(xyz);context.applyConstraints(1e-8)
        exact=context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.nanometer)
        rounded=exact.astype(np.float32).astype(float)
        box=np.array([v.value_in_unit(u.nanometer) for v in system.getDefaultPeriodicBoxVectors()])
        assert np.max(np.abs(box-np.diag(np.diag(box))))<1e-12,"diagnostic requires orthogonal cell"
        def errors(coords):
            result=[]
            for i in range(system.getNumConstraints()):
                a,b,d=system.getConstraintParameters(i);d=d.value_in_unit(u.nanometer)
                delta=coords[a]-coords[b];delta-=np.diag(box)*np.rint(delta/np.diag(box))
                result.append(float(abs(np.linalg.norm(delta)-d)/d))
            return result
        reference=errors(exact);quantized=errors(rounded)
        record=dict(identifier=name,referenceSystemSHA256=sha(raw),constraints=len(reference),
            referenceMaximumRelativeError=max(reference),roundedMaximumRelativeError=max(quantized),
            roundedAboveRequestedTolerance=sum(x>request["configuration"]["constraintTolerance"] for x in quantized))
        write(args.out/f"{name}.json",dict(summary=record,referenceRelativeErrors=reference,roundedRelativeErrors=quantized,
            projectedPositionsNM=exact.tolist(),roundedPositionsNM=rounded.tolist()))
        records.append(record);print(json.dumps(record),flush=True)
    write(args.out/"summary.json",dict(generatorSHA256=sha(Path(__file__).read_bytes()),openmm=mm.__version__,cases=records,
        scope="rounding diagnostic only; no impossibility proof and no native acceptance-tolerance change"))


if __name__=="__main__":main()
