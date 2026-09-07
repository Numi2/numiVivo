#!/usr/bin/env python3
"""Exercise the actual built native CLI; Python is an offline test driver only."""
from __future__ import annotations
import hashlib
import json
import pathlib
import subprocess
import sys


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit("usage: check_refinement_cli.py /path/to/numivivo-chemistry output-directory")
    binary = pathlib.Path(sys.argv[1]).resolve(strict=True)
    out = pathlib.Path(sys.argv[2]).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit("Use an empty output directory so stale artifacts cannot satisfy checks.")
    out.mkdir(parents=True, exist_ok=True)
    store = out / "artifacts"
    checks: list[str] = []
    runs: list[dict[str, object]] = []

    def execute(name: str, args: list[str], expected: int = 0) -> None:
        completed = subprocess.run([str(binary), *args], capture_output=True, text=True, timeout=180)
        (out / f"{name}.stdout.log").write_text(completed.stdout)
        (out / f"{name}.stderr.log").write_text(completed.stderr)
        runs.append({"name": name, "arguments": args, "exitCode": completed.returncode})
        if completed.returncode != expected:
            raise AssertionError(f"{name}: expected exit {expected}, got {completed.returncode}: {completed.stderr}")

    def path(name: str) -> str:
        return str(out / name)

    def read(name: str) -> dict:
        return json.loads((out / name).read_text())

    def write(name: str, value: object) -> None:
        (out / name).write_text(json.dumps(value, sort_keys=True))

    def require(condition: bool, label: str) -> None:
        if not condition:
            raise AssertionError(label)
        checks.append(label)

    execute("template", ["chemistry-template", "algebraic-space-refinement", "--output", path("request.json")])
    execute("refine", ["chemistry-refine", path("request.json"), "--store", str(store), "--output", path("result.json")])
    result = read("result.json")
    require(result["finalSpace"]["active"] == [0, 2], "CLI computes property-directed space")
    require(result["confirmation"]["passed"], "CLI retains confirmation result")
    execute("reuse", ["chemistry-refine", path("request.json"), "--store", str(store), "--output", path("reused.json")])
    require(all(node["reused"] for node in read("reused.json.receipt.json")["nodes"]), "validated cache reuse")
    require(read("reused.json") == result, "cache replay preserves numerical result")
    execute("export", ["chemistry-export-space", path("result.json"), "--point", "barrier-point", "--store", str(store), "--output", path("anchor.json")])
    anchor = read("anchor.json")
    require(anchor["orbitalIdentifiers"] == ["required", "reaction-sensitive"], "export materializes chosen active space")
    require(anchor["provenance"]["refinementEvidenceSHA256"] == hashlib.sha256((out / "result.json").read_bytes()).hexdigest(), "export binds exact canonical evidence bytes")
    execute("solver-template", ["chemistry-template", "solver-selected-ci", "--output", path("solver.json")])
    execute("solve", ["chemistry-solve", path("anchor.json"), "--solver", path("solver.json"), "--store", str(store), "--output", path("electronic.json")])
    electronic = read("electronic.json")
    require(electronic["converged"], "exported Hamiltonian runs through existing solve command")
    expected = next(p["energyHartree"] for p in result["confirmation"]["chosen"]["points"] if p["pointIdentifier"] == "barrier-point")
    require(abs(electronic["variationalEnergyHartree"] - expected) < 1e-9, "exported physical scalar and profile energy agree")
    write("state.json", electronic["state"])
    execute("pairs-template", ["chemistry-template", "orbital-information-selection", "--output", path("pairs.json")])
    execute("information", ["chemistry-correlations", path("state.json"), "--selection", path("pairs.json"), "--store", str(store), "--output", path("information.json")])
    require(read("information.json")["completePairCoverage"], "explicit two-orbital query computes the pair")
    execute("singles", ["chemistry-correlations", path("state.json"), "--store", str(store), "--output", path("singles.json")])
    require(not read("singles.json")["completePairCoverage"], "default single-orbital query does not assert pair coverage")
    request = read("request.json"); request["maximumPointEvaluations"] = 1
    write("capped.json", request)
    execute("cap", ["chemistry-refine", path("capped.json"), "--store", str(store), "--output", path("capped-result.json")], 2)
    require(read("capped-result.json")["termination"] == "incompleteExecution", "incomplete exploration retained with nonzero status")
    execute("reject-export", ["chemistry-export-space", path("capped-result.json"), "--point", "barrier-point", "--output", path("rejected-anchor.json")], 1)
    require(not (out / "rejected-anchor.json").exists(), "incomplete report cannot create a Hamiltonian")
    execute("reject-unseen", ["chemistry-export-space", path("result.json"), "--point", "new-geometry", "--output", path("unseen.json")], 1)
    require(not (out / "unseen.json").exists(), "unseen geometry is not implicitly qualified")
    solver = read("solver.json"); solver["configuration"]["maximumDeterminants"] = 1
    solver["configuration"]["selectionBatchSize"] = 1; solver["requireConverged"] = False
    write("probe.json", solver)
    execute("probe", ["chemistry-solve", path("anchor.json"), "--solver", path("probe.json"), "--store", str(store), "--output", path("probe-result.json")], 2)
    require(not read("probe-result.json")["converged"], "bounded probe remains explicitly unconverged")
    solver["requireConverged"] = True; write("strict.json", solver)
    execute("strict", ["chemistry-solve", path("anchor.json"), "--solver", path("strict.json"), "--store", str(store), "--output", path("strict-result.json")], 1)
    require(not (out / "strict-result.json").exists(), "accepted-solver policy does not publish a failed eigenpair")
    solver["schema"] = "unknown/schema"; write("unknown.json", solver)
    execute("unknown", ["chemistry-solve", path("anchor.json"), "--solver", path("unknown.json"), "--output", path("unknown-result.json")], 1)
    before = (out / "request.json").read_bytes()
    execute("overwrite", ["chemistry-refine", path("request.json"), "--output", path("request.json")], 1)
    require((out / "request.json").read_bytes() == before, "output cannot replace its request")
    execute("molecular-template", ["chemistry-template", "molecular-h2-space-preparation", "--output", path("molecular.json")])
    execute("prepare-molecule", ["chemistry-prepare-space", path("molecular.json"), "--store", str(store), "--output", path("molecular-request.json")])
    require(read("molecular-request.json")["initialSpace"]["active"] == [0, 1], "native molecular projection seeds complete H2 space")
    require(len(read("molecular-request.json.receipt.json")["nodes"]) == 11, "molecular CLI reuses separately addressable integral and HF nodes")
    execute("prepare-molecule-reuse", ["chemistry-prepare-space", path("molecular.json"), "--store", str(store), "--output", path("molecular-reused.json")])
    require(all(node["reused"] for node in read("molecular-reused.json.receipt.json")["nodes"]), "molecular primitive and projection caches revalidate before reuse")
    execute("molecular-refine", ["chemistry-refine", path("molecular-request.json"), "--store", str(store), "--output", path("molecular-result.json")])
    require(read("molecular-result.json")["confirmation"]["passed"], "prepared molecular request runs existing refinement")
    execute("multistate-template", ["chemistry-template", "algebraic-multistate-refinement", "--output", path("multistate.json")])
    execute("multistate-refine", ["chemistry-refine", path("multistate.json"), "--store", str(store), "--output", path("multistate-result.json")])
    require(read("multistate-result.json")["finalSpace"]["active"] == [0, 1, 2, 3], "multistate CLI rejects mean-energy cancellation")
    write("integration-results.json", {"scope": "actual native scoped CLI; not full Apple/Metal product or chemical qualification", "binarySHA256": hashlib.sha256(binary.read_bytes()).hexdigest(), "checks": checks, "runs": runs})
    print(f"{len(checks)} native refinement CLI checks passed")


if __name__ == "__main__":
    main()
