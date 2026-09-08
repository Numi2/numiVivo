#!/usr/bin/env python3
"""Derive a named native execution variant without recalculating reference physics.

Only declared native execution controls and the requested case subset may
change. Geometry, force-field models, observations and acceptance limits remain
byte-equivalent JSON values. Parent preparation failures remain in the panel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--references",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--position-precision",choices=["fp32","compensated"])
    p.add_argument("--ensemble",choices=["nve","nvt"])
    p.add_argument("--time-step-ps",type=float)
    p.add_argument("--temperature-k",type=float)
    p.add_argument("--friction-per-ps",type=float)
    p.add_argument("--seed",type=int)
    p.add_argument("--constraint-iterations",type=int)
    p.add_argument("--cases",nargs="+")
    p.add_argument("--steps",type=int)
    p.add_argument("--observe-every",type=int)
    p.add_argument("--project-initial-constraints",action="store_true")
    a=p.parse_args()
    if a.time_step_ps is not None and not 0<a.time_step_ps<=0.01:p.error("timestep must be 0...0.01 ps")
    if a.steps is not None and not 0<=a.steps<=100_000:p.error("steps must be 0...100000")
    if a.observe_every is not None and a.observe_every<=0:p.error("observation interval must be positive")
    if a.temperature_k is not None and not 0<a.temperature_k<=10000:p.error("temperature must be 0...10000 K")
    if a.friction_per_ps is not None and not 0<a.friction_per_ps<=10000:p.error("friction must be 0...10000 per ps")
    if a.seed is not None and not 0<=a.seed<2**64:p.error("seed must fit UInt64")
    if a.constraint_iterations is not None and not 1<=a.constraint_iterations<=16384:p.error("constraint iterations must be 1...16384")
    if a.ensemble=="nve" and (a.temperature_k is not None or a.friction_per_ps is not None):p.error("temperature/friction overrides require NVT")
    source=a.references.resolve(strict=True)
    manifest=json.loads((source/"manifest.json").read_text())
    assert manifest["schema"]=="numivivo.org/md-reference-campaign/v1"
    if a.cases is not None:
        if len(a.cases)!=len(set(a.cases)) or not set(a.cases)<=set(c["identifier"] for c in manifest["cases"]):p.error("case subset has duplicate or unknown identifiers")
        manifest["excludedPreparedCases"]=[c for c in manifest["cases"] if c["status"]=="prepared" and c["identifier"] not in a.cases]
        manifest["cases"]=[c for c in manifest["cases"] if c["status"]!="prepared" or c["identifier"] in a.cases]
    a.out.mkdir(parents=True,exist_ok=False)
    changes={k:v for k,v in vars(a).items() if k not in ("references","out") and v is not None}
    for case in manifest["cases"]:
        name=case["identifier"]
        assert name and all(c.isalnum() or c in "-_" for c in name)
        if case["status"]!="prepared":continue
        root=source/name
        assert digest(root/"request.json")==case["requestSHA256"]
        assert digest(root/"reference-system.xml")==case["serializedSystemSHA256"]
        shutil.copytree(root,a.out/name)
        request=json.loads((root/"request.json").read_text());config=request["configuration"]
        if a.position_precision is not None:config["positionPrecision"]=a.position_precision
        if a.time_step_ps is not None:config["timeStepPS"]=a.time_step_ps
        if a.ensemble is not None:
            config.update(ensemble=a.ensemble,barostat="none")
            if a.ensemble=="nve":config.update(thermostat="none",targetTemperatureK=None,frictionPerPS=None)
            else:config.update(thermostat="langevinMiddle",targetTemperatureK=300,frictionPerPS=1)
        if a.temperature_k is not None or a.friction_per_ps is not None:
            if config["ensemble"]!="nvt":p.error("temperature/friction overrides require an NVT case")
            if a.temperature_k is not None:config["targetTemperatureK"]=a.temperature_k
            if a.friction_per_ps is not None:config["frictionPerPS"]=a.friction_per_ps
        if a.seed is not None:config["randomSeed"]=a.seed
        if a.constraint_iterations is not None:config["maximumConstraintIterations"]=a.constraint_iterations
        if a.steps is not None:request["dynamicsSteps"]=a.steps
        if a.observe_every is not None:request["dynamicsObserveEvery"]=a.observe_every
        if a.project_initial_constraints:request["dynamicsPreparation"]="projectConstraints"
        path=a.out/name/"request.json"
        path.write_text(json.dumps(request,indent=2,allow_nan=False)+"\n")
        case["parentRequestSHA256"]=case["requestSHA256"];case["requestSHA256"]=digest(path)
    manifest["derivation"]=dict(parentManifestSHA256=digest(source/"manifest.json"),
        generatorSHA256=digest(Path(__file__)),nativeExecutionOverrides=changes)
    (a.out/"manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    print(json.dumps(dict(out=str(a.out),cases=len(manifest["cases"]),changes=changes)))


if __name__=="__main__":main()
