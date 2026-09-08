#!/usr/bin/env python3
"""Execute a fixed NVT matrix with immutable inputs and complete failure retention."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
from ensemble_campaign import read, digest, write, validate_policy, rigid_triangle_degrees_of_freedom, campaign_grid, analysis_environment, verify_runtime


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("references","binary","binary-manifest","out"):parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--policy",type=Path,default=Path(__file__).with_name("ensemble_policy.json"))
    a=parser.parse_args();policy=read(a.policy);validate_policy(policy)
    refs=a.references.resolve(strict=True);binary=a.binary.resolve(strict=True);tools=Path(__file__).parent
    parent=read(refs/"manifest.json");byid={r["identifier"]:r for r in parent["cases"]}
    for name in policy["cases"]:
        if name not in byid or byid[name]["status"]!="prepared":raise ValueError("requested ensemble case is not prepared: "+name)
        request=refs/name/"request.json"
        if digest(request)!=byid[name]["requestSHA256"]:raise ValueError("altered parent request")
        native_request=read(request);rigid_triangle_degrees_of_freedom(native_request["system"])
        if native_request["configuration"]["constraintTolerance"]!=policy["constraintTolerance"]:raise ValueError("parent constraint tolerance differs from policy")
    manifest=read(a.binary_manifest);environment=analysis_environment()
    verify_runtime(binary,manifest);a.out.mkdir(parents=True,exist_ok=False)
    for source,name in ((a.policy,"policy.json"),(a.binary_manifest,"binary-manifest.json"),(refs/"manifest.json","parent-manifest.json")):
        with (a.out/name).open("xb") as output:output.write(source.read_bytes())
    if read(a.out/"policy.json")!=policy or read(a.out/"binary-manifest.json")!=manifest or read(a.out/"parent-manifest.json")!=parent:
        raise ValueError("campaign input changed while taking its snapshot")
    identities=dict(binarySHA256=digest(binary),binaryManifestSHA256=digest(a.binary_manifest),nativeSource=manifest["buildSource"],analysisEnvironment=environment,
        parentManifestSHA256=digest(refs/"manifest.json"),policySHA256=digest(a.policy),tools={p.name:digest(p) for p in
        [tools/n for n in ("run_ensemble.py","analyze_ensemble.py","ensemble_campaign.py","ensemble_statistics.py","derive_campaign.py","run_campaign.py")]})
    write(a.out/"identities.json",identities)
    failures=[r for r in parent["cases"] if r["status"]!="prepared"]
    excluded=[r for r in parent["cases"] if r["status"]=="prepared" and r["identifier"] not in policy["cases"]]
    rows=[]
    for ordinal,(dt,temperature,seed) in enumerate(campaign_grid(policy)):
        directory=a.out/f"run-{ordinal:02d}";directory.mkdir();derived=directory/"references";output=directory/"native"
        steps=round(policy["durationPS"]/dt);interval=round(policy["observationIntervalPS"]/dt)
        derive=[sys.executable,str(tools/"derive_campaign.py"),"--references",str(refs),"--out",str(derived),
                "--cases",*policy["cases"],"--ensemble","nvt","--temperature-k",str(temperature),"--friction-per-ps",str(policy["frictionPerPS"]),
                "--time-step-ps",str(dt),"--seed",str(seed),"--constraint-iterations",str(policy["maximumConstraintIterations"]),
                "--position-precision","compensated","--steps",str(steps),"--observe-every",str(interval),"--project-initial-constraints"]
        run=[sys.executable,str(tools/"run_campaign.py"),"--references",str(derived),"--binary",str(binary),"--out",str(output),"--timeout","1800"]
        write(directory/"commands.json",dict(derive=derive,run=run));begin=time.monotonic()
        row=dict(directory=directory.name,timeStepPS=dt,temperatureK=temperature,seed=seed)
        try:
            with (directory/"derive.log").open("x") as log:subprocess.run(derive,stdout=log,stderr=subprocess.STDOUT,check=True)
            with (directory/"run.log").open("x") as log:completed=subprocess.run(run,stdout=log,stderr=subprocess.STDOUT)
            score=read(output/"scorecard.json")
            row.update(exitCode=completed.returncode,scorecardSHA256=digest(output/"scorecard.json"),referenceManifestSHA256=digest(derived/"manifest.json"),
                       cases=[dict(identifier=r["identifier"],outcome=r["outcome"]) for r in score["cases"]])
        except (OSError,ValueError,subprocess.CalledProcessError) as error:row.update(error=f"{type(error).__name__}: {error}",cases=[])
        row["wallSeconds"]=time.monotonic()-begin;rows.append(row);write(directory/"result.json",row)
        print(json.dumps(row),flush=True)
    verify_runtime(binary,manifest)
    for name,value in identities["tools"].items():
        if digest(tools/name)!=value:raise ValueError("campaign implementation changed during execution")
    result=dict(schema="numivivo.org/md-ensemble-campaign/v1",identities=identities,policy=policy,runs=rows,inputFailures=failures,excludedPreparedCases=excluded)
    write(a.out/"campaign.json",result)
    return 0 if not failures and all(r.get("cases") and all(c["outcome"]=="passed" for c in r["cases"]) for r in rows) else 1


if __name__=="__main__":raise SystemExit(main())
