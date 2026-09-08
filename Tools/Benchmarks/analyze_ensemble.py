#!/usr/bin/env python3
"""Independently qualify an immutable NVT matrix; never overwrite an analysis."""
import argparse
import itertools
import json
from pathlib import Path
from ensemble_campaign import read, digest, write, validate_policy, verify_series, campaign_grid, temperature_seeds, analysis_environment
from ensemble_statistics import kinetic_assessment, slope_assessment


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--campaign",type=Path,required=True);parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args();root=a.campaign.resolve(strict=True);campaign=read(root/"campaign.json");policy=campaign["policy"];validate_policy(policy)
    if campaign["schema"]!="numivivo.org/md-ensemble-campaign/v1":raise ValueError("campaign schema")
    if analysis_environment()!=campaign["identities"]["analysisEnvironment"]:raise ValueError("analysis environment differs from the recorded campaign")
    for name,key in (("policy.json","policySHA256"),("binary-manifest.json","binaryManifestSHA256"),("parent-manifest.json","parentManifestSHA256")):
        if digest(root/name)!=campaign["identities"][key]:raise ValueError("changed campaign input snapshot")
    if read(root/"policy.json")!=policy:raise ValueError("campaign policy differs from its original snapshot")
    tools=Path(__file__).parent
    for name in ("analyze_ensemble.py","ensemble_campaign.py","ensemble_statistics.py","run_campaign.py"):
        if digest(tools/name)!=campaign["identities"]["tools"][name]:raise ValueError("analysis differs from the recorded campaign implementation")
    expected=set(campaign_grid(policy))
    indexed={}
    for row in campaign["runs"]:
        key=(row["timeStepPS"],row["temperatureK"],row["seed"])
        if key not in expected or key in indexed:raise ValueError("duplicate or unrequested matrix cell")
        if Path(row["directory"]).name!=row["directory"]:raise ValueError("unsafe run directory")
        indexed[key]=row
    results=[];traces={};errors={};contracts=set()
    for case in policy["cases"]:
        for key in sorted(expected):
            try:
                row=indexed[key];directory=root/row["directory"]
                score=directory/"native/scorecard.json";manifest=directory/"references/manifest.json"
                if digest(score)!=row["scorecardSHA256"] or digest(manifest)!=row["referenceManifestSHA256"]:raise ValueError("changed run scorecard or reference manifest")
                sc=read(score);reference=next(c for c in read(manifest)["cases"] if c["identifier"]==case)
                item=next(c for c in sc["cases"] if c["identifier"]==case)
                if sc["binarySHA256"]!=campaign["identities"]["binarySHA256"] or item["outcome"]!="passed":raise ValueError("native case did not pass with the recorded binary")
                request=directory/"references"/case/"request.json";report=directory/"native"/case/"report.json"
                if digest(request)!=reference["requestSHA256"] or digest(request)!=item["requestSHA256"] or digest(report)!=item["reportSHA256"]:raise ValueError("changed request or report")
                if digest(directory/"references"/case/"reference-system.xml")!=reference["serializedSystemSHA256"]:raise ValueError("changed independent reference system")
                value=read(report);contracts.add(value["numericalContract"])
                traces[(case,*key)]=verify_series(read(request),value,policy,key[1],key[2],key[0])
            except (KeyError,StopIteration,OSError,ValueError,AssertionError) as error:errors[(case,*key)]=f"{type(error).__name__}: {error}"
    if len(contracts)>1:raise ValueError("matrix mixed native numerical contracts")
    for case in policy["cases"]:
        if len({v["hamiltonianIdentity"] for key,v in traces.items() if key[0]==case})>1:
            raise ValueError("matrix changed the Hamiltonian, reference geometry or scientific limits")
    for case in policy["cases"]:
        for dt_index,dt in enumerate(policy["timeStepsPS"]):
            row=dict(identifier=case,timeStepPS=dt,outcome="failed",kinetic=[])
            try:
                for temperature_index,temperature in enumerate(policy["temperaturesK"]):
                    for seed in temperature_seeds(policy,temperature):
                        key=(case,dt,temperature,seed)
                        if key in errors:raise ValueError(errors[key])
                    values=[traces[(case,dt,temperature,seed)] for seed in temperature_seeds(policy,temperature)]
                    dofs={v["degreesOfFreedom"] for v in values}
                    if len(dofs)!=1:raise ValueError("inconsistent independent degree count")
                    row["kinetic"].append(kinetic_assessment([v["kinetic"] for v in values],dofs.pop(),temperature,policy,dt_index*10+temperature_index))
                low,high=policy["temperaturesK"]
                row["configurational"]=slope_assessment([traces[(case,dt,low,s)]["potential"] for s in temperature_seeds(policy,low)],
                    [traces[(case,dt,high,s)]["potential"] for s in temperature_seeds(policy,high)],low,high,policy,100+dt_index*10)
                outcomes=[r["outcome"] for r in row["kinetic"]]+[row["configurational"]["outcome"]]
                row["outcome"]="failed" if "failed" in outcomes else ("inconclusive" if "inconclusive" in outcomes else "passed")
            except (KeyError,ValueError,AssertionError) as error:row["error"]=f"{type(error).__name__}: {error}"
            results.append(row);print(json.dumps(row),flush=True)
    prepared_equal=True
    for case,temperature in itertools.product(policy["cases"],policy["temperaturesK"]):
        for seed in temperature_seeds(policy,temperature):
            pair=[traces.get((case,dt,temperature,seed)) for dt in policy["timeStepsPS"]]
            if any(v is None for v in pair) or pair[0]["prepared"]!=pair[1]["prepared"]:prepared_equal=False
    prepared_outcome="passed" if prepared_equal and all(r["outcome"]=="passed" for r in results) else (
        "failed" if not prepared_equal or any(r["outcome"]=="failed" for r in results) else "inconclusive")
    result=dict(schema="numivivo.org/md-ensemble-qualification/v1",campaignSHA256=digest(root/"campaign.json"),identities=campaign["identities"],policy=policy,
                preparedWordsEqualAcrossTimeSteps=prepared_equal,preparedOutcome=prepared_outcome,cases=results,inputFailures=campaign["inputFailures"],
                excludedPreparedCases=campaign["excludedPreparedCases"],passed=prepared_outcome=="passed" and not campaign["inputFailures"],scope=policy["scope"])
    write(a.out,result);return 0 if result["passed"] else 1


if __name__=="__main__":raise SystemExit(main())
