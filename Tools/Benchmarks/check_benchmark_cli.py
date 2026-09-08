#!/usr/bin/env python3
"""Exercise benchmark document admission, no-clobber publication and false references."""
import argparse
import copy
import json
from pathlib import Path
import subprocess


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary",type=Path,required=True)
    parser.add_argument("--reference",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    original=json.loads(args.reference.read_text());checks=[]
    def run(name,request,expected,alias=False,existing=False):
        directory=args.out/name;directory.mkdir();source=directory/"request.json";target=directory/"report.json"
        source.write_text(json.dumps(request,allow_nan=False));before=source.read_bytes()
        if alias:target.symlink_to(source.resolve())
        elif existing:target.write_text("preserve me")
        command=[str(args.binary.resolve()),"md-benchmark",str(source),"--output",str(target)]
        result=subprocess.run(command,capture_output=True,timeout=60)
        (directory/"stdout.log").write_bytes(result.stdout);(directory/"stderr.log").write_bytes(result.stderr)
        assert result.returncode==expected,(name,result.returncode,result.stderr.decode())
        assert source.read_bytes()==before
        if existing:assert target.read_text()=="preserve me"
        elif expected==65 and not alias:assert not target.exists()
        elif expected==75:assert json.loads(target.read_text())["outcome"]=="failed"
        checks.append(dict(name=name,passed=True,exitCode=result.returncode))
    run("existing-output",original,65,existing=True)
    run("symlink-input-alias",original,65,alias=True)
    mutations=[("schema",lambda r:r.update(schema="unknown")),
               ("work-bound",lambda r:r.update(dynamicsSteps=100001)),
               ("observation-interval",lambda r:r.update(dynamicsSteps=10,dynamicsObserveEvery=0)),
               ("observation-budget",lambda r:r.update(dynamicsSteps=1001,dynamicsObserveEvery=1)),
               ("unknown-precision",lambda r:r["configuration"].update(positionPrecision="fp128")),
               ("missing-force",lambda r:r["references"][0].update(forcesKJPerMolNM=[])),
               ("negative-limit",lambda r:r["limits"].update(forceNormalizedRMS=-1)),
               ("non-fp32-coordinate",lambda r:r["references"][0]["geometry"]["particlePositionsNM"][0].update(x=0.1)),
               ("duplicate-observation",lambda r:r["references"].append(copy.deepcopy(r["references"][0])))]
    for name,change in mutations:
        request=copy.deepcopy(original);change(request);run(name,request,65)
    request=copy.deepcopy(original);request["references"][0]["energyKJPerMol"]+=100
    run("altered-reference-fails",request,75)
    (args.out/"checks.json").write_text(json.dumps(dict(passed=True,checks=checks),indent=2)+"\n")
    print(f"{len(checks)} benchmark CLI checks passed")


if __name__=="__main__":main()
