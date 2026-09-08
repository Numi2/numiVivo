#!/usr/bin/env python3
"""Extend every finer-step replica through exact checkpoints and bounded CLI runs.

The original benchmark matrix is immutable. Continuation uses the same owning
Metal runtime via md-run --restore. It never rethermalizes, discards failed
evidence, changes scientific acceptance limits, or fabricates a protocol cursor.
"""
import argparse
import copy
import math
from pathlib import Path
import subprocess
import time
from importlib import metadata

from ensemble_campaign import (read, digest, write, validate_policy, verify_series,
    campaign_grid, analysis_environment, verify_runtime, rigid_triangle_degrees_of_freedom)
from run_campaign import verify_checkpoint
from audit_endpoints import kinetic_check, reference_potential
from ensemble_statistics import kinetic_assessment, slope_assessment

CONTRACT = "numivivo.org/md-metal-numerics/v6"
SOURCES = ("continue_ensemble.py", "ensemble_campaign.py", "ensemble_statistics.py",
           "run_campaign.py", "audit_endpoints.py")


def valid_component(value):
    return isinstance(value, str) and bool(value) and all(c.isalnum() or c in "-_" for c in value)


def continuation_schedule(base, policy):
    validate_policy(base)
    if policy.get("schema") != "numivivo.org/md-ensemble-continuation-policy/v1":
        raise ValueError("continuation policy schema")
    if set(policy) != {"schema", "timeStepPS", "totalDurationPS", "chunkDurationPS", "bootstrapBlockPS", "reason", "scope"}:
        raise ValueError("continuation permits duration/block changes only")
    if policy["timeStepPS"] != base["timeStepsPS"][-1]:
        raise ValueError("continuation must include the original finer time step")
    for key in ("totalDurationPS", "chunkDurationPS", "bootstrapBlockPS"):
        value = policy[key]
        if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
            raise ValueError("invalid continuation duration: " + key)
    dt = policy["timeStepPS"]
    remaining = policy["totalDurationPS"] - base["durationPS"]
    chunk = round(policy["chunkDurationPS"] / dt)
    every = round(base["observationIntervalPS"] / dt)
    steps = round(remaining / dt)
    if not 0 < steps <= 10_000_000 or not 0 < chunk <= 100_000 or chunk // every > 1000:
        raise ValueError("bounded continuation capacity")
    for value in (remaining, policy["chunkDurationPS"], policy["bootstrapBlockPS"]):
        n = round(value / base["observationIntervalPS"])
        if not math.isclose(n * base["observationIntervalPS"], value, abs_tol=1e-12):
            raise ValueError("continuation must end on observation boundaries")
    production = round((policy["totalDurationPS"] - base["equilibrationPS"]) / base["observationIntervalPS"])
    block = round(policy["bootstrapBlockPS"] / base["observationIntervalPS"])
    if block <= 0 or production % block or production // block < base["minimumBlocksPerReplica"]:
        raise ValueError("continuation cannot meet complete-block minimum")
    chunks = [min(chunk, steps - start) for start in range(0, steps, chunk)]
    if any(n % every for n in chunks):
        raise ValueError("partial continuation observation interval")
    return chunks, every


def check_clock(actual, step, dt):
    expected = step * dt
    if type(step) is not int or not math.isfinite(actual) or abs(actual - expected) > max(1e-12, step * math.ulp(expected)):
        raise ValueError("continuation physical clock mismatch")


def verify_chunk(request, start, report, steps, every, policy):
    """Validate clocks, exact checkpoint words, constraints and stored observables."""
    if report.get("schema") != "numivivo.org/md-run-report/v1":
        raise ValueError("continuation report schema")
    if report.get("rejected") or report["requestedSteps"] != steps or report["committedSteps"] != steps:
        raise ValueError("incomplete continuation dynamics")
    final = report["finalCheckpoint"]
    dt = request["configuration"]["timeStepPS"]
    if report["startStep"] != start["acceptedStep"] or report["startTimePS"] != start["timePS"]:
        raise ValueError("continuation reset or changed the input clock")
    if report["endStep"] != final["acceptedStep"] or report["endTimePS"] != final["timePS"] or final["acceptedStep"] != start["acceptedStep"] + steps:
        raise ValueError("continuation final checkpoint clock mismatch")
    for checkpoint in (start, final):
        if checkpoint["numericalContract"] != CONTRACT:
            raise ValueError("continuation requires the v6 checkpoint contract")
        for key in ("systemFingerprint", "configurationFingerprint"):
            if report[key] != start[key] or checkpoint[key] != start[key]:
                raise ValueError("continuation state identity changed")
        check_clock(checkpoint["timePS"], checkpoint["acceptedStep"], dt)
        if checkpoint.get("periodicCell") != request["references"][0]["geometry"].get("periodicCell"):
            raise ValueError("continuation changed the NVT cell")
        if not verify_checkpoint(request, checkpoint)["passed"]:
            raise ValueError("independent continuation checkpoint constraints failed")
    samples = report["samples"]
    if len(samples) != steps // every:
        raise ValueError("missing or duplicate continuation observations")
    if steps == 0:
        if final != start:
            raise ValueError("zero-step restore did not preserve the exact checkpoint")
        return dict(kinetic=[], potential=[], endpointKinetic=None)
    dof = rigid_triangle_degrees_of_freedom(request["system"])
    kinetic, potential = [], []
    for ordinal, sample in enumerate(samples, 1):
        expected_step = start["acceptedStep"] + ordinal * every
        observation = sample["observables"]
        if sample["stepIndex"] != expected_step or observation["stepIndex"] != expected_step:
            raise ValueError("continuation observation step mismatch")
        if sample["timePS"] != observation["timePS"]:
            raise ValueError("continuation snapshot and observable clocks disagree")
        check_clock(sample["timePS"], expected_step, dt)
        if sample.get("periodicCell") != start.get("periodicCell"):
            raise ValueError("continuation snapshot cell changed")
        positions = sample["positionsNM"]
        if len(positions) != len(request["system"]["particles"]) or not all(math.isfinite(p[k]) for p in positions for k in ("x", "y", "z")):
            raise ValueError("invalid continuation snapshot positions")
        if observation["degreesOfFreedom"] != dof or any(observation[k] != start[k] for k in ("systemFingerprint", "configurationFingerprint")):
            raise ValueError("continuation observable identity or degree count changed")
        k, u = observation["kineticEnergyKJPerMol"], observation["potentialEnergyKJPerMol"]
        total, temperature = observation["totalEnergyKJPerMol"], observation["temperatureK"]
        if not all(math.isfinite(v) for v in (k, u, total, temperature)) or k < 0:
            raise ValueError("nonfinite or negative continuation observable")
        if not math.isclose(total, k + u, rel_tol=1e-12, abs_tol=1e-9) or not math.isclose(temperature, 2*k/(dof*policy["boltzmannKJPerMolK"]), rel_tol=1e-12, abs_tol=1e-9):
            raise ValueError("continuation energy accounting or temperature mismatch")
        kinetic.append(k)
        potential.append(u)
    if samples[-1]["positionsNM"] != final["positionsNM"] or samples[-1]["timePS"] != final["timePS"]:
        raise ValueError("final continuation snapshot differs from the checkpoint")
    audit = kinetic_check(request["system"], final, samples[-1]["observables"])
    if not audit["passed"]:
        raise ValueError("independent continuation endpoint kinetic energy failed")
    return dict(kinetic=kinetic, potential=potential, endpointKinetic=audit)


def load_parent(root, qualification_path):
    campaign = read(root / "campaign.json")
    if campaign["schema"] != "numivivo.org/md-ensemble-campaign/v1":
        raise ValueError("parent campaign schema")
    base = campaign["policy"]
    validate_policy(base)
    if analysis_environment() != campaign["identities"]["analysisEnvironment"]:
        raise ValueError("parent analysis environment changed")
    for name, key in (("policy.json", "policySHA256"), ("binary-manifest.json", "binaryManifestSHA256"), ("parent-manifest.json", "parentManifestSHA256")):
        if digest(root / name) != campaign["identities"][key]:
            raise ValueError("parent input snapshot changed")
    if read(root / "policy.json") != base:
        raise ValueError("parent policy changed")
    qualification = read(qualification_path)
    if qualification["campaignSHA256"] != digest(root / "campaign.json"):
        raise ValueError("parent must have a completed, bound qualification, including failures")
    for name in ("ensemble_campaign.py", "ensemble_statistics.py", "run_campaign.py"):
        if digest(Path(__file__).with_name(name)) != campaign["identities"]["tools"][name]:
            raise ValueError("parent validation implementation changed")
    grid = campaign_grid(base)
    keys = [(r["timeStepPS"], r["temperatureK"], r["seed"]) for r in campaign["runs"]]
    if sorted(keys) != sorted(grid):
        raise ValueError("parent matrix incomplete or duplicated")
    directories = [r.get("directory") for r in campaign["runs"]]
    if any(not valid_component(name) for name in directories) or len(set(directories)) != len(directories):
        raise ValueError("unsafe or duplicated parent run directory")
    return campaign, qualification


def parent_cell(root, campaign, row, case):
    if not valid_component(row["directory"]):
        raise ValueError("unsafe parent directory")
    directory = root / row["directory"]
    score_path = directory / "native/scorecard.json"
    manifest_path = directory / "references/manifest.json"
    if digest(score_path) != row["scorecardSHA256"] or digest(manifest_path) != row["referenceManifestSHA256"]:
        raise ValueError("parent scorecard or manifest changed")
    score = read(score_path)
    reference = next(c for c in read(manifest_path)["cases"] if c["identifier"] == case)
    native = next(c for c in score["cases"] if c["identifier"] == case)
    if score["binarySHA256"] != campaign["identities"]["binarySHA256"] or native["outcome"] != "passed":
        raise ValueError("parent native execution did not pass with the bound runtime")
    paths = dict(request=directory / "references" / case / "request.json",
                 report=directory / "native" / case / "report.json",
                 xml=directory / "references" / case / "reference-system.xml")
    if digest(paths["request"]) != reference["requestSHA256"] or digest(paths["request"]) != native["requestSHA256"] or digest(paths["report"]) != native["reportSHA256"] or digest(paths["xml"]) != reference["serializedSystemSHA256"]:
        raise ValueError("parent native or reference artifact changed")
    request, report = read(paths["request"]), read(paths["report"])
    if report["numericalContract"] != CONTRACT:
        raise ValueError("parent numerical contract changed")
    trace = verify_series(request, report, campaign["policy"], row["temperatureK"], row["seed"], row["timeStepPS"])
    return paths, request, report, trace


def file_manifest(directory):
    return {str(p.relative_to(directory)): digest(p) for p in directory.rglob("*") if p.is_file()}


def verify_files(directory, files):
    for name, sha in files.items():
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts or digest(directory / relative) != sha:
            raise ValueError("changed continuation artifact: " + name)


def run(args):
    parent = args.parent.resolve(strict=True)
    campaign, qualification = load_parent(parent, args.parent_qualification)
    base, policy = campaign["policy"], read(args.policy)
    chunks, every = continuation_schedule(base, policy)
    binary = args.binary.resolve(strict=True)
    manifest = read(parent / "binary-manifest.json")
    verify_runtime(binary, manifest)
    if digest(binary) != campaign["identities"]["binarySHA256"]:
        raise ValueError("continuation changed the parent executable")
    args.out.mkdir(parents=True, exist_ok=False)
    tools = Path(__file__).parent
    identities = dict(parentCampaignSHA256=digest(parent / "campaign.json"),
        parentQualificationSHA256=digest(args.parent_qualification),
        environment=dict(**analysis_environment(), openmm=metadata.version("openmm")),
        tools={n: digest(tools/n) for n in SOURCES}, binarySHA256=digest(binary),
        nativeSource=manifest["buildSource"], binaryManifest=manifest)
    write(args.out / "policy.json", policy)
    (args.out / "parent-campaign.json").write_bytes((parent / "campaign.json").read_bytes())
    (args.out / "parent-qualification.json").write_bytes(args.parent_qualification.read_bytes())
    if read(args.out/"parent-campaign.json") != campaign or read(args.out/"parent-qualification.json") != qualification:
        raise ValueError("parent changed during its continuation snapshot")
    write(args.out / "identities.json", identities)
    snapshots = file_manifest(args.out)
    rows = []
    for source in campaign["runs"]:
        if source["timeStepPS"] != policy["timeStepPS"]:
            continue
        for case in base["cases"]:
            directory = args.out / (source["directory"] + "-" + case)
            directory.mkdir()
            row = dict(directory=directory.name, identifier=case, timeStepPS=source["timeStepPS"],
                temperatureK=source["temperatureK"], seed=source["seed"], outcome="failed", chunks=[])
            begin = time.monotonic()
            try:
                paths, request, original, _ = parent_cell(parent, campaign, source, case)
                for key, path in paths.items():
                    (directory / ("parent-" + path.name)).write_bytes(path.read_bytes())
                if read(directory/"parent-request.json") != request or read(directory/"parent-report.json") != original:
                    raise ValueError("parent changed while copying its inputs")
                if digest(directory/"parent-reference-system.xml") != digest(paths["xml"]):
                    raise ValueError("reference changed while copying its input")
                write(directory / "system.json", request["system"])
                write(directory / "configuration.json", request["configuration"])
                start = original["dynamics"]["finalCheckpoint"]
                checkpoint_path = directory / "start-checkpoint.json"
                write(checkpoint_path, start)
                inputs = file_manifest(directory)
                # A zero-step round trip checks exact high/low words, velocities,
                # cell, fingerprints and stochastic step counter before extension.
                for index, steps in enumerate([0] + chunks):
                    chunk_dir = directory / f"chunk-{index:03d}"
                    chunk_dir.mkdir()
                    report_path = chunk_dir / "report.json"
                    next_checkpoint = chunk_dir / "checkpoint.json"
                    command = [str(binary), "md-run", str(directory / "system.json"),
                        "--config", str(directory / "configuration.json"), "--restore", str(checkpoint_path),
                        "--steps", str(steps), "--sample-every", str(every), "--observables-every", str(every),
                        "--max-samples", "1000", "--checkpoint", str(next_checkpoint), "--output", str(report_path)]
                    command_record = dict(command=command, restoreSHA256=digest(checkpoint_path), requestedSteps=steps)
                    write(chunk_dir / "command.json", command_record)
                    with (chunk_dir / "native.log").open("x") as log:
                        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
                    if completed.returncode != 0:
                        raise ValueError("native continuation exit " + str(completed.returncode))
                    report = read(report_path)
                    verified = verify_chunk(request, start, report, steps, every, base)
                    if read(next_checkpoint) != report["finalCheckpoint"]:
                        raise ValueError("separate final checkpoint differs from report")
                    audit = dict(endpointKinetic=verified["endpointKinetic"])
                    if steps:
                        energy, images = reference_potential((directory/"parent-reference-system.xml").read_text(), report["finalCheckpoint"], request["system"])
                        actual = report["samples"][-1]["observables"]["potentialEnergyKJPerMol"]
                        error = abs(energy-actual) / len(request["system"]["particles"])
                        limit = request["limits"]["energyAbsolutePerParticleKJPerMol"]
                        audit["endpointPotential"] = dict(referenceEnergyKJPerMol=energy, nativeEnergyKJPerMol=actual,
                            errorPerParticleKJPerMol=error, limitPerParticleKJPerMol=limit, images=images, passed=error <= limit)
                    write(chunk_dir / "audit.json", audit)
                    if steps and not audit["endpointPotential"]["passed"]:
                        raise ValueError("independent continuation endpoint potential energy failed")
                    row["chunks"].append(dict(directory=chunk_dir.name, requestedSteps=steps,
                        restoreSHA256=command_record["restoreSHA256"], files=file_manifest(chunk_dir)))
                    checkpoint_path, start = next_checkpoint, report["finalCheckpoint"]
                    print(f"{directory.name} {chunk_dir.name}: {start['acceptedStep']} steps, {start['timePS']:.3f} ps", flush=True)
                    del report, verified
                verify_files(directory, inputs)
                row["outcome"] = "passed"
            except Exception as error:
                row["error"] = f"{type(error).__name__}: {error}"
            row["wallSeconds"] = time.monotonic() - begin
            row["files"] = file_manifest(directory)
            write(directory / "result.json", row)
            rows.append(row)
    verify_runtime(binary, manifest)
    for name, sha in identities["tools"].items():
        if digest(tools/name) != sha:
            raise ValueError("continuation implementation changed during execution")
    verify_files(args.out, snapshots)
    result = dict(schema="numivivo.org/md-ensemble-continuation/v1", policy=policy,
        identities=identities, snapshots=snapshots, parentPolicy=base, runs=rows,
        inputFailures=campaign["inputFailures"], excludedPreparedCases=campaign["excludedPreparedCases"],
        parentPreparedOutcome=qualification["preparedOutcome"], scope=policy["scope"])
    write(args.out / "campaign.json", result)
    return 0 if rows and all(r["outcome"] == "passed" for r in rows) else 1


def analyze(args):
    root = args.campaign.resolve(strict=True)
    campaign = read(root / "campaign.json")
    if campaign["schema"] != "numivivo.org/md-ensemble-continuation/v1":
        raise ValueError("continuation campaign schema")
    verify_files(root, campaign["snapshots"])
    base, policy = campaign["parentPolicy"], campaign["policy"]
    chunks, every = continuation_schedule(base, policy)
    if read(root/"policy.json") != policy or read(root/"parent-campaign.json")["policy"] != base:
        raise ValueError("continuation policy changed")
    identities = campaign["identities"]
    if digest(root/"parent-campaign.json") != identities["parentCampaignSHA256"] or digest(root/"parent-qualification.json") != identities["parentQualificationSHA256"]:
        raise ValueError("continuation parent evidence binding changed")
    if read(root/"identities.json") != identities or identities["environment"] != dict(**analysis_environment(), openmm=metadata.version("openmm")):
        raise ValueError("continuation identity or environment changed")
    for name, sha in identities["tools"].items():
        if digest(Path(__file__).with_name(name)) != sha:
            raise ValueError("continuation analysis implementation changed")
    expected = {(case, dt, t, s) for case in base["cases"] for dt, t, s in campaign_grid(base) if dt == policy["timeStepPS"]}
    indexed = {}
    for row in campaign["runs"]:
        key = (row["identifier"], row["timeStepPS"], row["temperatureK"], row["seed"])
        if key not in expected or key in indexed or not valid_component(row["directory"]):
            raise ValueError("invalid continuation matrix cell")
        indexed[key] = row
    if set(indexed) != expected:
        raise ValueError("incomplete continuation matrix")
    traces, errors = {}, []
    for key, row in indexed.items():
        try:
            directory = root / row["directory"]
            verify_files(directory, row["files"])
            if row["outcome"] != "passed" or [r["requestedSteps"] for r in row["chunks"]] != [0] + chunks:
                raise ValueError("native continuation incomplete or failed")
            request, original = read(directory / "parent-request.json"), read(directory / "parent-report.json")
            trace = verify_series(request, original, base, key[2], key[3], key[1])
            checkpoint = directory / "start-checkpoint.json"
            start = original["dynamics"]["finalCheckpoint"]
            if read(checkpoint) != start or read(directory/"system.json") != request["system"] or read(directory/"configuration.json") != request["configuration"]:
                raise ValueError("continuation changed its original state or configuration")
            for index, chunk in enumerate(row["chunks"]):
                if chunk["directory"] != f"chunk-{index:03d}" or digest(checkpoint) != chunk["restoreSHA256"]:
                    raise ValueError("broken continuation checkpoint chain")
                subdir = directory / chunk["directory"]
                verify_files(subdir, chunk["files"])
                report = read(subdir / "report.json")
                values = verify_chunk(request, start, report, chunk["requestedSteps"], every, base)
                checkpoint = subdir / "checkpoint.json"
                if read(checkpoint) != report["finalCheckpoint"]:
                    raise ValueError("changed continuation checkpoint")
                audit = read(subdir/"audit.json")
                if audit["endpointKinetic"] != values["endpointKinetic"] or (index and not audit["endpointPotential"]["passed"]):
                    raise ValueError("continuation endpoint audit failed")
                trace["kinetic"].extend(values["kinetic"])
                trace["potential"].extend(values["potential"])
                start = report["finalCheckpoint"]
                del report, values
            expected_samples = round((policy["totalDurationPS"] - base["equilibrationPS"]) / base["observationIntervalPS"])
            if len(trace["kinetic"]) != expected_samples or start["acceptedStep"] != round(policy["totalDurationPS"] / key[1]):
                raise ValueError("continuation final duration or production sample count mismatch")
            traces[key] = trace
        except Exception as error:
            errors.append(dict(identifier=key[0], timeStepPS=key[1], temperatureK=key[2], seed=key[3], error=f"{type(error).__name__}: {error}"))
    effective = copy.deepcopy(base)
    effective.update(schema="numivivo.org/md-ensemble-continuation-analysis-policy/v1",
        timeStepsPS=[policy["timeStepPS"]], durationPS=policy["totalDurationPS"],
        bootstrapBlockPS=policy["bootstrapBlockPS"], scope=policy["scope"])
    results = []
    for case in base["cases"]:
        row = dict(identifier=case, timeStepPS=policy["timeStepPS"], outcome="failed")
        try:
            selected = {key: traces[key] for key in expected if key[0] == case}
            if len({v["hamiltonianIdentity"] for v in selected.values()}) != 1:
                raise ValueError("continuation matrix changed its Hamiltonian")
            dofs = {v["degreesOfFreedom"] for v in selected.values()}
            if len(dofs) != 1:
                raise ValueError("continuation degree counts disagree")
            dof = dofs.pop()
            rows = [[v for key, v in sorted(selected.items()) if key[2] == t] for t in base["temperaturesK"]]
            row["kinetic"] = [kinetic_assessment([v["kinetic"] for v in values], dof, t, effective, 200+i) for i, (t, values) in enumerate(zip(base["temperaturesK"], rows))]
            row["configurational"] = slope_assessment([v["potential"] for v in rows[0]], [v["potential"] for v in rows[1]], *base["temperaturesK"], effective, 300)
            outcomes = [r["outcome"] for r in row["kinetic"]] + [row["configurational"]["outcome"]]
            row["outcome"] = "failed" if "failed" in outcomes else ("inconclusive" if "inconclusive" in outcomes else "passed")
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
        results.append(row)
    outcome = "failed" if errors or any(r["outcome"] == "failed" for r in results) else ("inconclusive" if any(r["outcome"] == "inconclusive" for r in results) else "passed")
    result = dict(schema="numivivo.org/md-ensemble-continuation-qualification/v1", campaignSHA256=digest(root/"campaign.json"),
        preparedOutcome=outcome, cases=results, validationErrors=errors, policy=policy, effectiveStatisticsPolicy=effective,
        parentPreparedOutcome=campaign["parentPreparedOutcome"], inputFailures=campaign["inputFailures"],
        passed=outcome == "passed" and not campaign["inputFailures"], scope=policy["scope"])
    write(args.out, result)
    return 0 if result["passed"] else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    runner = sub.add_parser("run")
    for name in ("parent", "parent-qualification", "binary", "out"):
        runner.add_argument("--" + name, type=Path, required=True)
    runner.add_argument("--policy", type=Path, default=Path(__file__).with_name("ensemble_continuation_policy.json"))
    analyzer = sub.add_parser("analyze")
    for name in ("campaign", "out"):
        analyzer.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    return run(args) if args.action == "run" else analyze(args)


if __name__ == "__main__":
    raise SystemExit(main())
