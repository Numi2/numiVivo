#!/usr/bin/env python3
"""Run the published reference panel and retain every failure and raw result."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def digest(path):
    result=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):result.update(block)
    return result.hexdigest()


def write(path,value):
    with path.open("x") as f:
        json.dump(value,f,indent=2,allow_nan=False);f.write("\n")


def verify_comparisons(request,report):
    refs=request["references"];limits=request["limits"]
    assert len(report["evaluations"])==len(refs)==len(report["comparisons"]),"missing comparisons"
    outcomes=[]
    for ref,value,comparison in zip(refs,report["evaluations"],report["comparisons"]):
        r=[v[k] for v in ref["forcesKJPerMolNM"] for k in ("x","y","z")]
        a=[v[k] for v in value["physicalParticleForcesKJPerMolNM"] for k in ("x","y","z")]
        assert len(r)==len(a)>0
        scale=max(limits["forceNormalizationFloor"],math.sqrt(math.fsum(x*x for x in r)/len(r)))
        errors=[x-y for x,y in zip(a,r)]
        rms=math.sqrt(math.fsum(x*x for x in errors)/len(r))/scale
        maximum=max(map(abs,errors))/scale
        energy=abs(value["energyKJPerMol"]-ref["energyKJPerMol"])/(len(r)/3)
        # Swift optional encoding omits nil keys. Normalize both representations.
        same={k:v for k,v in value["evaluatedGeometry"].items() if v is not None}=={k:v for k,v in ref["geometry"].items() if v is not None}
        passed=same and energy<=limits["energyAbsolutePerParticleKJPerMol"] and rms<=limits["forceNormalizedRMS"] and maximum<=limits["forceNormalizedMaximum"]
        assert comparison["identifier"]==ref["identifier"]
        assert comparison["outcome"]==("passed" if passed else "failed")
        for key,expected in [("energyErrorPerParticleKJPerMol",energy),("forceNormalizedRMS",rms),("forceNormalizedMaximum",maximum)]:
            assert math.isclose(comparison[key],expected,rel_tol=1e-10,abs_tol=1e-12),(key,expected,comparison[key])
        outcomes.append(passed)
    return all(outcomes)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references",type=Path,required=True)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--timeout",type=int,default=900)
    args=parser.parse_args()
    if not 1<=args.timeout<=86400:parser.error("invalid timeout")
    binary=args.binary.resolve(strict=True);references=args.references.resolve(strict=True)
    binary_hash=digest(binary)
    manifest=json.loads((references/"manifest.json").read_text())
    assert manifest["schema"]=="numivivo.org/md-reference-campaign/v1"
    args.out.mkdir(parents=True,exist_ok=False)
    rows=[]
    for case in manifest["cases"]:
        name=case["identifier"]
        assert name and all(c.isalnum() or c in "-_" for c in name),"unsafe case identifier"
        row=dict(identifier=name,referenceStatus=case["status"],ensemble="inconclusive",experimentalAccuracy="not-evaluated")
        if case["status"]!="prepared":
            row.update(outcome="preparation-failed",error=case.get("error"));rows.append(row);continue
        case_dir=references/name;request_path=case_dir/"request.json"
        assert digest(request_path)==case["requestSHA256"],"altered reference request"
        assert digest(case_dir/"reference-system.xml")==case["serializedSystemSHA256"],"altered reference system"
        request=json.loads(request_path.read_text())
        assert request["referenceProvenance"]["serializedSystemSHA256"]==case["serializedSystemSHA256"]
        out=args.out/name;out.mkdir();report_path=out/"report.json"
        command=[str(binary),"md-benchmark",str(request_path),"--output",str(report_path)]
        write(out/"command.json",command)
        begin=time.monotonic();code=None
        with (out/"stdout.log").open("wb") as stdout,(out/"stderr.log").open("wb") as stderr:
            try:
                code=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=args.timeout,env={**os.environ,"MTL_DEBUG_LAYER":"1"}).returncode
            except subprocess.TimeoutExpired:row["error"]="native command exceeded fixed wall-time limit"
        row.update(exitCode=code,wallSeconds=time.monotonic()-begin,particles=case["particles"],requestSHA256=digest(request_path))
        if report_path.exists():
            report=json.loads(report_path.read_text())
            assert report["identifier"]==name
            if report["outcome"]=="unsupported":
                assert not report["capability"]["executable"] and not report["comparisons"]
                row["blockers"]=report["capability"]["blockers"]
            else:
                static_pass=verify_comparisons(request,report)
                row["hamiltonianAgreement"]="passed" if static_pass else "failed"
                row["comparisons"]=report["comparisons"]
                dynamics=report.get("dynamics")
                if dynamics:
                    assert dynamics["ensembleOutcome"]=="inconclusive"
                    row["committedSteps"]=dynamics["committedSteps"]
                    row["requestedSteps"]=dynamics["requestedSteps"]
                expected=static_pass and not report.get("executionError") and (not dynamics or not dynamics.get("rejected"))
                assert report["outcome"]==("passed" if expected else "failed")
            row.update(outcome=report["outcome"],device=report.get("deviceName"),reportSHA256=digest(report_path))
            assert (code==0)==(report["outcome"]=="passed")
        else:row["outcome"]="execution-failed"
        rows.append(row);print(json.dumps(row),flush=True)
    assert digest(binary)==binary_hash,"executable changed during campaign"
    summary=dict(schema="numivivo.org/md-benchmark-campaign/v1",binarySHA256=binary_hash,
        referenceManifestSHA256=digest(references/"manifest.json"),runnerSHA256=digest(Path(__file__)),
        cases=rows,passed=all(r["outcome"]=="passed" for r in rows),
        scope="static Hamiltonian agreement and short execution smoke only; no ensemble or speed leadership claim")
    write(args.out/"scorecard.json",summary)
    return 0 if summary["passed"] else 1


if __name__=="__main__":raise SystemExit(main())
