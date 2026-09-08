#!/usr/bin/env python3
"""Measure a preregistered, matched-duration constrained NVE refinement panel."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


def read(path):return json.loads(path.read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write(path,value):
    with path.open("x") as f:json.dump(value,f,indent=2,allow_nan=False);f.write("\n")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--references",type=Path,required=True)
    p.add_argument("--binary",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--policy",type=Path,default=Path(__file__).with_name("nve_policy.json"))
    a=p.parse_args();policy=read(a.policy);tools=Path(__file__).parent
    assert policy["schema"]=="numivivo.org/md-nve-refinement-policy/v1"
    assert policy["timeStepsPS"]==sorted(policy["timeStepsPS"],reverse=True) and len(policy["timeStepsPS"])==3
    assert policy["initialTemperatureK"]==300
    a.out.mkdir(parents=True,exist_ok=False)
    write(a.out/"preregistered-policy.json",policy)
    identities=dict(policySHA256=digest(a.policy),generatorSHA256=digest(tools/"derive_campaign.py"),
        runnerSHA256=digest(tools/"run_campaign.py"),qualifierSHA256=digest(Path(__file__)),
        binarySHA256=digest(a.binary),parentManifestSHA256=digest(a.references/"manifest.json"))
    write(a.out/"identities.json",identities)
    runs=[]
    for index,dt in enumerate(policy["timeStepsPS"]):
        steps=round(policy["durationPS"]/dt);interval=round(policy["observationIntervalPS"]/dt)
        assert 0<steps<=100000 and interval>0 and steps%interval==0
        assert math.isclose(steps*dt,policy["durationPS"],abs_tol=1e-14)
        assert math.isclose(interval*dt,policy["observationIntervalPS"],abs_tol=1e-14)
        refs=a.out/f"references-{index}";out=a.out/f"run-{index}"
        derive=[sys.executable,str(tools/"derive_campaign.py"),"--references",str(a.references),"--out",str(refs),
            "--position-precision","compensated","--ensemble","nve","--time-step-ps",str(dt),"--steps",str(steps),
            "--observe-every",str(interval),"--project-initial-constraints"]
        command=[sys.executable,str(tools/"run_campaign.py"),"--references",str(refs),"--binary",str(a.binary),"--out",str(out)]
        write(a.out/f"commands-{index}.json",dict(derive=derive,run=command))
        with (a.out/f"run-{index}.log").open("w") as log:
            subprocess.run(derive,stdout=log,stderr=subprocess.STDOUT,check=True)
            completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        runs.append(dict(dt=dt,out=out,exitCode=completed.returncode,scorecard=read(out/"scorecard.json")))
    results=[]
    for case in read(a.references/"manifest.json")["cases"]:
        row=dict(identifier=case["identifier"],outcome="failed",scope=policy["scope"])
        if case["status"]!="prepared":row.update(outcome="preparation-failed",error=case.get("error"));results.append(row);continue
        try:
            measurements=[];prepared=[]
            for run in runs:
                score=next(r for r in run["scorecard"]["cases"] if r["identifier"]==case["identifier"])
                assert score["outcome"]=="passed",score
                report=read(run["out"]/case["identifier"]/"report.json");dynamics=report["dynamics"]
                samples=dynamics["observations"];assert len(samples)==round(policy["durationPS"]/policy["observationIntervalPS"])+1
                assert samples[0]==dynamics["start"] and samples[-1]==dynamics["end"],"series boundaries differ from start/end observations"
                cp=dynamics["preparedCheckpoint"]
                prepared.append({k:cp[k] for k in ("positionsNM","positionHighNM","positionCorrectionsNM","velocitiesNMPerPS","acceptedStep","timePS")})
                for i,sample in enumerate(samples):
                    assert math.isclose(sample["timePS"]-samples[0]["timePS"],i*policy["observationIntervalPS"],abs_tol=1e-12)
                deviations=[(s["totalEnergyKJPerMol"]-samples[0]["totalEnergyKJPerMol"])/case["particles"] for s in samples]
                measurements.append(dict(timeStepPS=run["dt"],committedSteps=dynamics["committedSteps"],
                    maximumAbsoluteDeviationPerParticle=max(map(abs,deviations)),
                    rmsDeviationPerParticle=math.sqrt(math.fsum(x*x for x in deviations)/len(deviations)),
                    energyDeviationsPerParticle=deviations,reportSHA256=digest(run["out"]/case["identifier"]/"report.json")))
            assert prepared[0]==prepared[1]==prepared[2],"timestep variants do not share exact prepared position and velocity words"
            ratios=[];refinement=[]
            for coarse,fine in zip(measurements,measurements[1:]):
                c,f=coarse["rmsDeviationPerParticle"],fine["rmsDeviationPerParticle"]
                ratio=c/f if f>0 else None;ratios.append(ratio)
                refinement.append(c<=policy["refinementFloorPerParticleKJPerMol"] or f<=policy["refinementFloorPerParticleKJPerMol"] or ratio>=policy["minimumRMSRefinementRatio"])
            conservation=all(m["maximumAbsoluteDeviationPerParticle"]<=policy["maximumAbsoluteEnergyDeviationPerParticleKJPerMol"] for m in measurements)
            row.update(outcome="passed" if conservation and all(refinement) else "failed",conservationPassed=conservation,
                refinementPassed=all(refinement),rmsRefinementRatios=ratios,measurements=measurements)
        except (AssertionError,KeyError,TypeError,ValueError) as error:row["error"]=f"{type(error).__name__}: {error}"
        results.append(row);print(json.dumps(row),flush=True)
    write(a.out/"qualification.json",dict(schema="numivivo.org/md-nve-refinement-result/v1",identities=identities,
        policy=policy,cases=results,passed=all(r["outcome"]=="passed" for r in results)))
    return 0 if all(r["outcome"]=="passed" for r in results) else 1


if __name__=="__main__":raise SystemExit(main())
