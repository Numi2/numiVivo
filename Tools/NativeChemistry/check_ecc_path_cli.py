#!/usr/bin/env python3
"""Exercise real shared-path CLI, primitive cache reuse and integrity failures.

Uses only the standard library. The executable must contain the actual production
command classes, numerics and artifact store; no interface substitutes are made.
"""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def fingerprint(value):
    data = value["bytes"]
    if len(data) != 32 or any(not isinstance(x, int) or x < 0 or x > 255 for x in data):
        raise RuntimeError("Malformed artifact fingerprint in actual CLI receipt")
    return bytes(data).hex()


def run(binary, reference_dir, out):
    if out.exists() and any(out.iterdir()):
        raise RuntimeError("Use an empty output directory; the check never deletes previous user artifacts")
    out.mkdir(parents=True, exist_ok=True)
    checks = []
    counter = 0

    def require(condition, label):
        checks.append({"label": label, "passed": bool(condition)})
        write(out / "checks.json", {"schema": "numivivo.org/ecc-path-cli-checks/v1", "checks": checks,
              "passed": all(x["passed"] for x in checks),
              "executableSHA256": hashlib.sha256(binary.read_bytes()).hexdigest()})
        if not condition:
            raise RuntimeError(label)
        print("PASS", label, flush=True)

    def execute(*arguments, expect_failure=False):
        nonlocal counter
        counter += 1
        process = subprocess.run([str(binary), *map(str, arguments)], capture_output=True, text=True, timeout=600)
        (out / f"command-{counter:02d}.log").write_text(
            f"args: {json.dumps(list(map(str, arguments)))}\nexit: {process.returncode}\n"
            + process.stdout + "\n" + process.stderr)
        if (process.returncode == 0) == expect_failure:
            raise RuntimeError(f"Unexpected exit {process.returncode}: {arguments}\n{process.stderr}")
        return process

    request_path, store = out / "path.json", out / "artifacts"
    execute("chemistry-path-template", "h2-stretch", "--output", request_path)
    request = json.loads(request_path.read_text())

    def path_run(input_path, name):
        destination = out / f"{name}.json"
        execute("chemistry-path-run", input_path, "--store", store, "--output", destination)
        receipt = json.loads(Path(str(destination) + ".receipt.json").read_text())
        return json.loads(destination.read_text()), receipt

    result, receipt = path_run(request_path, "first")
    require(result["path"]["converged"] and len(result["path"]["pointResults"]) == 5,
            "five-geometry shared ECC executes through the production path command")
    require(len(receipt["nodes"]) == 11 and all(not node["reused"] for node in receipt["nodes"]),
            "first run materializes ten primitive stages and one shared-path stage")
    repeated, repeated_receipt = path_run(request_path, "repeated")
    require(all(node["reused"] for node in repeated_receipt["nodes"]) and repeated == result,
            "unchanged path reuses all verified stages and exports identical numerical content")
    require(receipt["result"] == repeated_receipt["result"], "unchanged path retains its content-addressed result identity")
    require(result["requestFingerprint"] == receipt["input"] and len(result["snapshotFingerprints"]) == 5,
            "export binds the exact request and all five snapshot identities")
    require(result["coordinates"] == [1.2, 1.4, 1.6, 1.8, 2.0]
            and result["energyMeaning"] == "electronic-profile; no transition-state characterization, thermal correction, standard-state correction or kinetic export",
            "ordered electronic profile has no implied barrier or kinetic qualification")
    reference_raw = (reference_dir / "reference.json").read_bytes()
    manifest = json.loads((reference_dir / "manifest.json").read_text())
    if hashlib.sha256(reference_raw).hexdigest() != manifest["sha256"]:
        raise RuntimeError("Independent reference fixture hash differs")
    oracle = json.loads(reference_raw)
    expected = oracle["points"]
    errors = [abs(actual["frame"]["energyHartree"] - reference["fci"])
              for actual, reference in zip(result["path"]["pointResults"], expected)]
    require(len(expected) == 5 and all(math.isfinite(x) and x <= 1e-8 for x in errors),
            "actual CLI profile energies agree with independent PySCF at every geometry")

    # Canonical implementation and primitive resource/configuration identities
    # must also agree with the existing isolated electronic workflow.
    isolated = {"schema": "numivivo.org/electronic-workflow/v1", "system": request["snapshots"][0]["system"],
                "basis": request["basis"], "budget": request["budget"],
                "calculation": {"hartreeFock": {"configuration": request["reference"]}}}
    isolated_path, isolated_output = out / "isolated.json", out / "isolated.result.json"
    write(isolated_path, isolated)
    execute("chemistry-run", isolated_path, "--store", store, "--output", isolated_output)
    isolated_receipt = json.loads(Path(str(isolated_output) + ".receipt.json").read_text())
    require(isolated_receipt["implementation"] == receipt["implementation"],
            "path and isolated commands have the same actual-executable implementation identity")
    require(len(isolated_receipt["nodes"]) == 2 and all(node["reused"] for node in isolated_receipt["nodes"]),
            "isolated RHF reuses the identical AO and reference artifacts produced by the path")

    changed = copy.deepcopy(request)
    changed["snapshots"][2]["coordinate"] = 1.65
    for nucleus, z in zip(changed["snapshots"][2]["system"]["nuclei"], [-0.825, 0.825]):
        nucleus["positionBohr"][2] = z
    changed_path = out / "one-geometry-changed.json"
    write(changed_path, changed)
    changed_result, changed_receipt = path_run(changed_path, "one-geometry-result")
    misses = {node["identifier"] for node in changed_receipt["nodes"] if not node["reused"]}
    require(misses == {"path-2-integrals", "path-2-reference", "path-profile"},
            "editing one geometry invalidates only its AO/RHF stages and the shared path")
    require(changed_result["snapshotFingerprints"][2] != result["snapshotFingerprints"][2]
            and all(changed_result["snapshotFingerprints"][i] == result["snapshotFingerprints"][i] for i in [0, 1, 3, 4]),
            "snapshot provenance identifies exactly which geometry changed")

    reweighted = copy.deepcopy(request)
    reweighted["configuration"]["pointWeights"] = [0.3, 0.1, 0.2, 0.2, 0.2]
    reweighted_path = out / "reweighted.json"
    write(reweighted_path, reweighted)
    _, weighted_receipt = path_run(reweighted_path, "reweighted.result")
    require({node["identifier"] for node in weighted_receipt["nodes"] if not node["reused"]} == {"path-profile"},
            "changing shared objective weights retains all geometry-specific preparation artifacts")

    def reject_request(value, stem):
        path, destination = out / f"{stem}.request.json", out / f"{stem}.output.json"
        write(path, value)
        execute("chemistry-path-run", path, "--store", store, "--output", destination, expect_failure=True)
        return not destination.exists() and not Path(str(destination) + ".receipt.json").exists()

    bad = copy.deepcopy(request); bad["schema"] = "numivivo.org/unknown-path/v99"
    require(reject_request(bad, "unknown-schema"), "unknown path schema fails without publishing a fallback result")
    bad = copy.deepcopy(request); bad["snapshots"][2]["system"]["alphaElectrons"] = 0
    require(reject_request(bad, "wrong-electrons"), "changed electron sector is rejected before result publication")
    bad = copy.deepcopy(request); bad["configuration"]["maximumPointEvaluations"] = 5
    require(reject_request(bad, "unfinished"), "insufficient shared-optimization budget cannot publish success")
    original = request_path.read_bytes()
    execute("chemistry-path-run", request_path, "--store", store, "--output", request_path, expect_failure=True)
    require(request_path.read_bytes() == original, "request/output aliasing is rejected before overwrite")
    alias = out / "request-symlink.json"; alias.symlink_to(request_path.resolve())
    execute("chemistry-path-run", request_path, "--store", store, "--output", alias, expect_failure=True)
    require(request_path.read_bytes() == original and alias.is_symlink(), "symlink alias cannot overwrite a request")
    alias = out / "request-hardlink.json"; os.link(request_path, alias)
    execute("chemistry-path-run", request_path, "--store", store, "--output", alias, expect_failure=True)
    require(request_path.read_bytes() == original and os.path.samefile(request_path, alias),
            "hard-link alias cannot replace request contents")
    destination = store / "user-result.json"
    execute("chemistry-path-run", request_path, "--store", store, "--output", destination, expect_failure=True)
    require(not destination.exists(), "output cannot overwrite a file inside the artifact store")
    execute("chemistry-path-run", request_path, "--solver", "fci", expect_failure=True)
    require(request_path.read_bytes() == original, "unrecognized path options are not interpreted as another method")

    digest = fingerprint(receipt["result"])
    object_path = store / "objects" / "sha256" / digest[:2] / digest[2:4] / digest
    stored = object_path.read_bytes()
    if hashlib.sha256(stored).hexdigest() != digest:
        raise RuntimeError("Unexpected artifact object layout or initial digest")
    try:
        object_path.write_bytes(stored + b"\n")
        rejected_output = out / "corrupt-cache-output.json"
        execute("chemistry-path-run", request_path, "--store", store, "--output", rejected_output, expect_failure=True)
        require(not rejected_output.exists(), "corrupted cached path object fails SHA-256 verification without an exported result")
    finally:
        object_path.write_bytes(stored)
    report = {"schema": "numivivo.org/ecc-path-cli-checks/v1", "checks": checks, "passed": True,
              "executableSHA256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "maximumPySCFProfileEnergyErrorHartree": max(errors),
              "scope": "real scoped production path/electronic command classes, numerical sources and artifact store; not the full biology/Metal product"}
    write(out / "checks.json", report)
    print(f"PASS {len(checks)} production shared-path CLI/store checks", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.binary.resolve(), args.reference.resolve(), args.output.resolve())
