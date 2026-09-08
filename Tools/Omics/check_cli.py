#!/usr/bin/env python3
"""Host-only checks of the real native CLI. No production calculations run in Python."""
import argparse
import copy
import hashlib
import json
import pathlib
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True, type=pathlib.Path)
    parser.add_argument("--out", required=True, type=pathlib.Path)
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    root = pathlib.Path(__file__).resolve().parents[2]
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    fixture = out / "inputs"
    shutil.copytree(root / "Examples/singlecell", fixture)
    store = out / "store"
    commands = []
    checks = []

    def run(*arguments, success=True):
        command = [str(binary)] + [str(x) for x in arguments]
        result = subprocess.run(command, text=True, capture_output=True, timeout=180)
        index = len(commands)
        (out / f"command-{index:03d}.stdout").write_text(result.stdout)
        (out / f"command-{index:03d}.stderr").write_text(result.stderr)
        commands.append({"arguments": command, "exitCode": result.returncode, "expectedSuccess": success})
        if (result.returncode == 0) != success:
            raise RuntimeError(f"Unexpected exit {result.returncode}: {command}\n{result.stderr[-8000:]}")
        return result.stdout

    def check(name, condition):
        if not condition:
            raise RuntimeError("Check failed: " + name)
        checks.append(name)

    def load(path):
        return json.loads(path.read_text())

    def save(path, value):
        path.write_text(json.dumps(value, indent=2) + "\n")
        return path

    def hex_id(fingerprint):
        return bytes(fingerprint["bytes"]).hex()

    def object_path(fingerprint, selected_store=store):
        descriptor = json.loads(run("artifact-show", hex_id(fingerprint), "--store", selected_store))
        return selected_store / descriptor["objectPath"]

    passed = False
    try:
        check("help exposes supported workflow", "singlecell-run" in run("singlecell-help"))
        manifest = fixture / "manifest.json"
        receipt_path = out / "receipt.json"
        run("singlecell-run", manifest, "--store", store, "--output", receipt_path)
        receipt = load(receipt_path)
        verification = json.loads(run("singlecell-verify", receipt_path, "--store", store))
        check("verified native reconstruction", verification["status"] == "verified-native-reconstruction-not-biological-validation")
        check("source evidence not promoted", verification["evidence"] == "synthetic")
        check("expected dimensions and replicate groups", (verification["cells"], verification["features"], verification["pseudobulkGroups"]) == (6, 3, 2))
        report_path = out / "report.json"
        run("singlecell-export", receipt_path, "--store", store, "--output", report_path)
        report = load(report_path)
        check("exact raw source counts", report["dataset"]["matrix"]["counts"] == [2, 3, 1, 7, 4, 6, 2, 14])
        check("exact per-cell totals", [x["totalCounts"] for x in report["quality"]] == [5, 8, 0, 10, 16, 0])
        check("empty-cell mitochondrial fraction remains missing", report["quality"][2].get("mitochondrialFraction") is None)
        check("pseudobulk retains replicate counts", report["pseudobulk"]["matrix"]["counts"] == [2, 4, 7, 4, 8, 14])
        check("normalized view separate and scale invariant", report["normalized"]["values"][:4] == report["normalized"]["values"][4:])
        check("cell identities retained", report["dataset"]["cells"][0]["sampleID"] != report["dataset"]["cells"][3]["sampleID"])
        second = out / "repeat.json"
        run("singlecell-run", manifest, "--store", store, "--output", second)
        check("repeat run has stable content identities", load(second) == receipt)
        original = receipt_path.read_bytes()
        run("singlecell-run", manifest, "--store", store, "--output", receipt_path, success=False)
        check("receipt not overwritten", receipt_path.read_bytes() == original)
        run("singlecell-run", manifest, "--store", store, "--normalise", "1", success=False)
        run("singlecell-run", manifest, "--store", store, "--output", store / "forbidden.json", success=False)
        alias = out / "store-alias"
        alias.symlink_to(store, target_is_directory=True)
        run("singlecell-run", manifest, "--store", store, "--output", alias / "forbidden.json", success=False)
        check("output cannot enter store through symlink", not (store / "forbidden.json").exists())
        unknown = load(manifest); unknown["normalisationTarget"] = 500
        run("singlecell-run", save(fixture / "unknown.json", unknown), "--store", store, success=False)
        wrong = copy.deepcopy(receipt); wrong["implementation"]["bytes"][0] ^= 1
        run("singlecell-verify", save(out / "wrong-implementation.json", wrong), "--store", store, success=False)
        wrong = copy.deepcopy(receipt); wrong["input"]["bytes"][0] ^= 1
        run("singlecell-verify", save(out / "wrong-input.json", wrong), "--store", store, success=False)
        # Store a correctly hashed but numerically false result; hash integrity
        # alone must not let it become an accepted result.
        record = load(object_path(receipt["result"]))
        record["report"]["quality"][0]["totalCounts"] += 1
        false_result = save(out / "false-result.json", record)
        descriptor = json.loads(run("artifact-put", false_result, "--kind", "vivo.singlecell-stored-result-v1",
                                    "--media-type", "application/json", "--store", store))
        wrong = copy.deepcopy(receipt); wrong["result"] = descriptor["fingerprint"]
        run("singlecell-verify", save(out / "false-receipt.json", wrong), "--store", store, success=False)
        # Exercise the pre-existing workflow scheduler with a stored typed input,
        # so UInt64 counts never pass through the generic JSON number type.
        recipe = {
            "schema": "numivivo.org/workflow-recipe/v1", "identifier": "singlecell-cli-conformance",
            "artifacts": [{"identifier": "counts", "source": {"stored": {
                "kind": "vivo.singlecell-input-bundle-v1", "fingerprint": receipt["input"]}}}],
            "nodes": [{"identifier": "analyze", "operation": "vivo.platform.singlecell", "version": "1",
                       "configuration": {}, "inputs": {"input": {"artifact": {"identifier": "counts"}}},
                       "resources": {"budget": {"maximumBasisFunctions": 64, "maximumBytes": 1073741824,
                                                "maximumDeterminants": 512, "maximumOperatorApplications": 100000000},
                                     "maximumInputBytes": 134217728, "maximumOutputBytes": 536870912,
                                     "numericalBackend": "cpu-fp64"}}],
            "outputs": [{"name": "report", "node": "analyze", "port": "report"}],
            "policy": {"maximumConcurrentMetalTasks": 0, "maximumConcurrentTasks": 1,
                       "maximumInlineBytes": 134217728, "maximumNodes": 256, "maximumReservedBytes": 2147483648}}
        recipe_path = save(out / "workflow.json", recipe)
        run("workflow-plan", recipe_path)
        dag_path = out / "workflow-run.json"
        run("workflow-run", recipe_path, "--store", store, "--output", dag_path)
        check("general workflow executed", load(dag_path)["allTasksSucceeded"] is True)
        dag_receipt = load(pathlib.Path(str(dag_path) + ".receipt.json"))
        exported_path = out / "workflow-report.json"
        run("workflow-export", hex_id(dag_receipt["reportArtifact"]), "--name", "report",
            "--store", store, "--output", exported_path)
        check("general workflow and dedicated command agree", load(exported_path) == report)
        replay_path = out / "workflow-repeat.json"
        run("workflow-run", recipe_path, "--store", store, "--output", replay_path)
        replay = load(replay_path)
        check("validated cache reuse", replay["allTasksSucceeded"] and all(n["outcome"]["succeeded"]["reused"] for n in replay["nodes"]))
        # Changing current files does not invalidate their already snapshotted
        # historical bytes; a new run must get a different identity.
        changed = fixture / "control-1.mtx"
        changed.write_text(changed.read_text().replace("1 1 2", "1 1 9"))
        run("singlecell-verify", receipt_path, "--store", store)
        changed_receipt = out / "changed.json"
        run("singlecell-run", manifest, "--store", store, "--output", changed_receipt)
        check("changed source cannot reuse old input identity", load(changed_receipt)["input"] != receipt["input"])
        damaged_store = out / "damaged-store"
        shutil.copytree(store, damaged_store)
        damaged_input = object_path(receipt["input"], damaged_store)
        damaged_input.write_bytes(damaged_input.read_bytes() + b" ")
        run("singlecell-verify", receipt_path, "--store", damaged_store, success=False)
        check("primary stored source preserved", hashlib.sha256(object_path(receipt["input"]).read_bytes()).hexdigest() == hex_id(receipt["input"]))
        passed = True
    finally:
        save(out / "observations.json", {"schema": "numivivo.org/singlecell-cli-checks/v1",
             "passed": passed, "checks": checks, "commands": commands,
             "binarySHA256": hashlib.sha256(binary.read_bytes()).hexdigest(),
             "scope": "Synthetic native CLI, artifact integrity, reconstruction and existing workflow composition; not biological validation."})
    print(f"Single-cell CLI checks passed: {len(checks)} assertions, {len(commands)} commands")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
