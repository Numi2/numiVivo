#!/usr/bin/env python3
"""Offline Arc 2026 evaluator adapter. No metric implementations or training here.

Requires a native NumiVivo submission and a clean pinned cell-eval2 checkout
installed in this Python environment. References/bundles are scoring-only inputs.
A local score is not an official leaderboard submission or acceptance receipt.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import sys
import tempfile

COMMIT = "5e64833518a6603a0301cbe28185d49c30f4a986"
VERSION = "0.16.0"
RULE_VERSION = 3
METRICS = (
    "pds_cosine", "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_direction_fidelity_yield_raw", "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard", "de_wilcoxon_lfc_nmae",
)
MAX_FILE = 64 * 1024**3
MAX_JSON = 16 * 1024**2


def require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def unique_object(pairs):
    obj = {}
    for key, value in pairs:
        require(key not in obj, f"duplicate JSON field: {key}")
        obj[key] = value
    return obj


def reject_constant(value):
    raise ValueError(f"nonfinite JSON constant: {value}")


def read_json(path: Path):
    with open_regular(path, MAX_JSON) as f:
        return json.loads(f.read(), object_pairs_hook=unique_object, parse_constant=reject_constant)


def write_json(path: Path, obj, *, replace: bool = False):
    data = (json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()
    if not replace:
        with path.open("xb") as f:
            f.write(data)
        return
    fd, tmp = tempfile.mkstemp(prefix=".summary-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data); f.flush(); os.fsync(f.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def open_regular(path: Path, limit=MAX_FILE):
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        info = os.fstat(fd)
        require(stat.S_ISREG(info.st_mode) and 0 <= info.st_size <= limit,
                f"not a bounded regular file: {path}")
        return os.fdopen(fd, "rb")
    except BaseException:
        os.close(fd)
        raise


def digest(path: Path, destination: Path | None = None, limit=MAX_FILE) -> str:
    h = hashlib.sha256()
    with open_regular(path, limit) as src:
        before = os.fstat(src.fileno())
        out = destination.open("xb") if destination is not None else None
        try:
            total = 0
            while block := src.read(1024**2):
                total += len(block)
                require(total <= limit, f"file grew beyond limit: {path}")
                h.update(block)
                if out:
                    out.write(block)
            after = os.fstat(src.fileno())
            require((before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                    (after.st_size, after.st_mtime_ns, after.st_ctime_ns) and total == before.st_size,
                    f"source changed during snapshot: {path}")
            if out:
                out.flush(); os.fsync(out.fileno())
        finally:
            if out:
                out.close()
    return h.hexdigest()


def safe_id(value):
    return isinstance(value, str) and 0 < len(value) <= 80 and all(
        c in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in value)


def checked_ids(values):
    require(isinstance(values, list) and len(values) == 3 and
            all(safe_id(x) for x in values) and len(set(values)) == 3,
            "a phase requires exactly three unique context IDs")
    return values


def snapshot_submission(source: Path, target: Path):
    """Copy all inference artifacts before any reference is read; never copy arbitrary paths."""
    target.mkdir()
    sh = digest(source / "submission.json", target / "submission.json", MAX_JSON)
    receipt = read_json(target / "submission.json")
    require(receipt["schemaVersion"] == 1 and receipt["evidenceClass"] == "input-conformance-only" and
            receipt["evaluatorCommit"] == COMMIT, "submission contract/evaluator mismatch")
    require(receipt["phase"] in ("validation", "final-test"), "unknown phase")
    ids = checked_ids([c["id"] for c in receipt["contexts"]])
    (target / "query").mkdir()
    require(digest(source / "query/query.json", target / "query/query.json", MAX_JSON) == receipt["querySHA256"],
            "query manifest hash mismatch")
    query = read_json(target / "query/query.json")
    require(query["schemaVersion"] == 1 and query["phase"] == receipt["phase"] and query["evaluatorCommit"] == COMMIT,
            "query/submission mismatch")
    qids = checked_ids([c["id"] for c in query["contexts"]])
    require(set(qids) == set(ids) == set(receipt["model"]["excludedPerturbedContexts"]) and
            len(receipt["model"]["excludedPerturbedContexts"]) == 3, "zero-shot context declaration mismatch")
    for c in query["contexts"]:
        path = Path("contexts") / c["id"]
        (target / "query" / path).mkdir(parents=True)
        for name, key in (("controls.h5ad", "controlsSHA256"), ("targets.json", "targetsSHA256")):
            require(digest(source / "query" / path / name, target / "query" / path / name) == c[key],
                    f"query payload hash mismatch: {c['id']}/{name}")
        (target / path).mkdir(parents=True)
        prediction = next(r for r in receipt["contexts"] if r["id"] == c["id"])
        require(digest(source / path / "prediction.h5ad", target / path / "prediction.h5ad") == prediction["predictionSHA256"],
                f"prediction hash mismatch: {c['id']}")
    require(digest(source / "training-data.json", target / "training-data.json", MAX_JSON) == receipt["model"]["trainingManifestSHA256"],
            "training-data manifest hash mismatch")
    return receipt, query, sh


def tree_hashes(root: Path, copy_to: Path | None = None):
    require(root.is_dir() and not root.is_symlink(), f"not a real directory: {root}")
    result = {}
    for p in sorted(root.rglob("*")):
        require(not p.is_symlink(), f"symlink in artifact tree: {p}")
        relative = p.relative_to(root)
        if p.is_dir():
            if copy_to:
                (copy_to / relative).mkdir(parents=True, exist_ok=True)
            continue
        require(len(result) < 4096, "artifact tree file-count resource limit")
        if copy_to:
            (copy_to / relative).parent.mkdir(parents=True, exist_ok=True)
        result[relative.as_posix()] = digest(p, copy_to / relative if copy_to else None)
    require(bool(result), "empty artifact tree")
    return result


def evaluator_source(root: Path, copy_to: Path | None = None):
    def git(*args):
        return subprocess.check_output(["git", "-C", str(root), *args])
    require(git("rev-parse", "HEAD").decode().strip() == COMMIT, "wrong cell-eval2 source commit")
    entries = git("ls-tree", "-r", "-z", "HEAD", "src/cell_eval2", "pyproject.toml", "LICENSE").split(b"\0")
    files = {}
    for entry in entries:
        if not entry:
            continue
        meta, path = entry.split(b"\t", 1)
        mode, kind, expected = meta.split()
        name = path.decode()
        require(kind == b"blob" and mode in (b"100644", b"100755"), "unsupported evaluator source entry")
        with open_regular(root / name, MAX_JSON) as f:
            data = f.read()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        require(blob == expected.decode(), f"modified evaluator source: {name}")
        files[name] = hashlib.sha256(data).hexdigest()
        if copy_to:
            destination = copy_to / name; destination.parent.mkdir(parents=True, exist_ok=True)
            with destination.open("xb") as out:
                out.write(data)
    require("src/cell_eval2/competition.py" in files, "evaluator source missing")
    # Reject even ignored/untracked importable source. Bytecode is disabled below.
    for p in (root / "src/cell_eval2").rglob("*"):
        require(not p.is_symlink(), "evaluator source symlink")
        if p.is_file() and "__pycache__" not in p.parts:
            require(p.relative_to(root).as_posix() in files, f"untracked evaluator package payload: {p}")
    return files


def python_command(root: Path, body: str, *args):
    # Isolated interpreter: no caller PYTHONPATH, user site, or working-directory import.
    code = "import sys; sys.dont_write_bytecode=True; sys.path.insert(0," + repr(str(root / "src")) + ");\n" + body
    return [sys.executable, "-I", "-B", "-c", code, *map(str, args)]


def execute(command, log: Path, timeout: int):
    with log.open("xb") as out:
        out.write((json.dumps({"argv": command}) + "\n").encode()); out.flush()
        env = dict(os.environ)
        env["NUMBA_CACHE_DIR"] = str(log.parent / "numba-cache")
        run = subprocess.run(command, cwd=log.parent, stdout=out, stderr=subprocess.STDOUT,
                             timeout=timeout, check=False, env=env)
    require(run.returncode == 0, f"command failed ({run.returncode}); retained log: {log}")


PROBE = '''import json, pathlib, importlib.metadata as md
import cell_eval2
from cell_eval2 import competition
from cell_eval2.config import EvalConfig
assert pathlib.Path(cell_eval2.__file__).resolve().parent == pathlib.Path(sys.path[0]).resolve()/"cell_eval2"
payload = competition.competition_payload()
print("ARC2026_PROBE=" + json.dumps({"version":md.version("cell-eval2"), "ruleVersion":payload["rule_version"],
 "ruleDigest":competition.competition_digest(), "members":list(competition.competition_members()),
 "config":EvalConfig.from_preset("vcc2026").to_dict(),
 "python":sys.version, "dependencies":sorted((d.metadata.get("Name", ""),d.version) for d in md.distributions())}, allow_nan=False))
'''

REFERENCE_AXIS_CHECK = '''import json
import anndata as ad
with open(sys.argv[1]) as f: q=json.load(f)
c=next(c for c in q["contexts"] if c["id"]==sys.argv[2])
a=ad.read_h5ad(sys.argv[3], backed="r")
try:
 assert list(a.var_names)==c["genes"], "reference gene axis mismatch"
 assert a.obs_names.is_unique, "duplicate reference observation identities"
 assert not a.obs[c["perturbationColumn"]].isna().any(), "missing reference labels"
 assert set(a.obs[c["perturbationColumn"]].tolist())==set(c["targets"]+["non-targeting"]), "reference panel mismatch"
 print(json.dumps({"referenceShape":list(a.shape), "axisVerified":True}))
finally:
 a.file.close()
'''


def read_scores(path: Path):
    with open_regular(path, MAX_JSON) as f:
        rows = list(csv.DictReader(f.read().decode().splitlines()))
    require(len(rows) == 7 and {r.get("metric") for r in rows} == set(METRICS) | {"avg_score"},
            "expected exactly six official metrics plus avg_score")
    require(all(r.get("anchor_source") == "real_bundle" and r.get("real_bundle_id") and
                r.get("real_bundle_digest") and r.get("anchor_digest") for r in rows),
            "scores are not real-bundle anchored")
    values = {r["metric"]: float(r["from_replicate"]) for r in rows}
    require(all(math.isfinite(x) for x in values.values()), "missing/nonfinite official score")
    average = next(r for r in rows if r["metric"] == "avg_score")
    require(average.get("from_baseline") == "", "diagnostic/baseline average cannot be a competition score")
    require(math.isclose(values["avg_score"], math.fsum(values[m] for m in METRICS) / 6,
                         rel_tol=1e-10, abs_tol=1e-12), "official six-metric average mismatch")
    # No global clamp: Arc applies each metric's own policy inside its scorer.
    return values


def aggregate(expected, completed):
    checked_ids(expected)
    require(set(completed).issubset(expected), "foreign completed context")
    missing = [x for x in expected if x not in completed]
    for c, scores in completed.items():
        require(set(scores) == set(METRICS) | {"avg_score"} and all(math.isfinite(x) for x in scores.values()),
                f"incomplete metric row: {c}")
        require(math.isclose(scores["avg_score"], math.fsum(scores[m] for m in METRICS) / 6,
                             rel_tol=1e-10, abs_tol=1e-12), f"metric mean mismatch: {c}")
    return {"status": "complete" if not missing else "incomplete", "missingContexts": missing,
            "contextScores": completed,
            "score": math.fsum(completed[x]["avg_score"] for x in expected) / 3 if not missing else None,
            "scope": "local-pinned-evaluator-not-leaderboard-acceptance"}


def run(args):
    output = args.output.absolute()
    output.mkdir(parents=False, exist_ok=False)
    expected, completed, errors = [], {}, {}
    try:
        native = args.native.absolute(); root = args.evaluator.absolute()
        native_hash = digest(native)
        receipt, query, sh = snapshot_submission(args.submission.absolute(), output / "submission")
        expected = [c["id"] for c in query["contexts"]]
        execute([str(native), "arc2026-verify", str(output / "submission")], output / "native-verify.log", args.timeout_seconds)
        require(digest(native) == native_hash, "native executable changed during verification")
        frozen = output / "evaluator-source"; frozen.mkdir()
        source_files = evaluator_source(root, frozen)
        # Execute only copied, tracked bytes. Existing checkout bytecode is never imported.
        execute(python_command(frozen, PROBE), output / "evaluator-probe.log", args.timeout_seconds)
        lines = (output / "evaluator-probe.log").read_text().splitlines()
        probes = [json.loads(x[len("ARC2026_PROBE="):]) for x in lines if x.startswith("ARC2026_PROBE=")]
        require(len(probes) == 1, "missing or ambiguous evaluator probe")
        probe = probes[0]
        require(probe["version"] == VERSION and probe["ruleVersion"] == RULE_VERSION and
                set(probe["members"]) == set(METRICS) and len(probe["members"]) == 6,
                "installed evaluator/version/rule mismatch")
        write_json(output / "runtime.json", {"evaluatorCommit": COMMIT, "sourceFiles": source_files,
                   "probe": probe, "nativeSHA256": native_hash, "submissionSHA256": sh,
                   "adapterSHA256": digest(Path(__file__)),
                   "nativeVerificationLogSHA256": digest(output / "native-verify.log"),
                   "evaluatorProbeLogSHA256": digest(output / "evaluator-probe.log")})
        # Outcomes are opened only AFTER all prediction bytes have been frozen.
        digest(args.references, output / "references-plan.json", MAX_JSON)
        plan = read_json(output / "references-plan.json")
        require(set(plan) == {"schemaVersion", "sourceDescription", "contexts"} and plan["schemaVersion"] == 1
                and isinstance(plan["sourceDescription"], str) and bool(plan["sourceDescription"]), "reference plan schema")
        refs = plan["contexts"]
        require(isinstance(refs, list) and len(refs) <= 3 and
                all(set(r) == {"id", "controlsSHA256", "targetsSHA256", "referenceH5AD",
                               "referenceSHA256", "realBundle", "realBundleManifestSHA256"} for r in refs), "reference contexts schema")
        require(len({r["id"] for r in refs}) == len(refs) and {r["id"] for r in refs}.issubset(expected),
                "duplicate/foreign scoring contexts")
        for ref in refs:
            cid = ref["id"]; work = output / cid; work.mkdir()
            try:
                context = next(c for c in query["contexts"] if c["id"] == cid)
                require(ref["controlsSHA256"] == context["controlsSHA256"] and
                        ref["targetsSHA256"] == context["targetsSHA256"], "reference plan binds a different query context")
                reference_hash = digest(Path(ref["referenceH5AD"]), work / "reference.h5ad")
                require(reference_hash == ref["referenceSHA256"], "reference file differs from the frozen release binding")
                bundle = work / "bundle"; bundle.mkdir()
                bundle_files = tree_hashes(Path(ref["realBundle"]), bundle)
                require(bundle_files.get("manifest.json") == ref["realBundleManifestSHA256"], "wrong scoring bundle for this release binding")
                manifest = read_json(bundle / "manifest.json")
                require(manifest.get("cell_eval2_version") == VERSION and manifest.get("rule_digest") == probe["ruleDigest"],
                        "bundle must match pinned 0.16.0/rule 3; no version restamping or diagnostic fallback")
                device, backend = manifest.get("resolved_device"), manifest.get("resolved_de_backend")
                require(isinstance(device, str) and bool(device) and isinstance(backend, str) and bool(backend),
                        "bundle lacks resolved device/DE backend")
                execute(python_command(frozen, REFERENCE_AXIS_CHECK, output / "submission/query/query.json", cid, work / "reference.h5ad"),
                        work / "reference-axis.log", args.timeout_seconds)
                run_dir = work / "metrics"
                cmd = python_command(frozen, "from cell_eval2.cli import main\nmain()\n", "run",
                    "-ap", output / "submission/contexts" / cid / "prediction.h5ad",
                    "-ar", work / "reference.h5ad", "--preset", "vcc2026",
                    "--pert-col", context["perturbationColumn"],
                    "--set", "device=" + json.dumps(device), "--set", "de.backend=" + json.dumps(backend), "-o", run_dir)
                execute(cmd, work / "metrics.log", args.timeout_seconds)
                score_cmd = python_command(frozen, "from cell_eval2.cli import main\nmain()\n", "score",
                    "--user-agg", run_dir / "agg_results.csv", "--user-meta", run_dir / "run_meta.json",
                    "--real-bundle", bundle, "--comparison-statistic", "mean", "--output", work / "scores.csv")
                execute(score_cmd, work / "score.log", args.timeout_seconds)
                values = read_scores(work / "scores.csv")
                require(tree_hashes(frozen) == source_files, "frozen evaluator changed during scoring")
                require(tree_hashes(bundle) == bundle_files and digest(work / "reference.h5ad") == reference_hash,
                        "scoring inputs changed")
                prediction_hash = next(c["predictionSHA256"] for c in receipt["contexts"] if c["id"] == cid)
                require(digest(output / "submission/contexts" / cid / "prediction.h5ad") == prediction_hash, "prediction changed during scoring")
                write_json(work / "receipt.json", {"schemaVersion": 1, "id": cid, "referenceSHA256": reference_hash,
                    "predictionSHA256": prediction_hash, "bundleFiles": bundle_files, "scores": values,
                    "scoresSHA256": digest(work / "scores.csv"), "metricsFiles": tree_hashes(run_dir),
                    "runtimeSHA256": digest(output / "runtime.json"),
                    "logs": {name: digest(work / name) for name in ("reference-axis.log", "metrics.log", "score.log")}})
                completed[cid] = values
            except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
                errors[cid] = str(exc)
                write_json(work / "failure.json", {"error": str(exc), "score": None})
            summary = aggregate(expected, completed)
            summary.update({"phase": query["phase"], "errors": errors})
            write_json(output / "summary.json", summary, replace=True)
        summary = aggregate(expected, completed)
        summary.update({"phase": query["phase"], "errors": errors})
        write_json(output / "summary.json", summary, replace=True)
        print(json.dumps(summary, indent=2, allow_nan=False))
        return 0 if summary["status"] == "complete" else 1
    except BaseException as exc:
        write_json(output / "failure.json", {"error": str(exc), "score": None})
        # A failure after a partial run must not leave a plausible complete score.
        write_json(output / "summary.json", {"status": "failed", "score": None,
                   "contextScores": completed, "missingContexts": [x for x in expected if x not in completed]}, replace=True)
        raise


def verify(args):
    """Artifact/aggregation replay, not a fresh re-evaluation of the six metrics."""
    root = args.directory
    runtime = read_json(root / "runtime.json")
    require(runtime["evaluatorCommit"] == COMMIT and runtime["probe"]["version"] == VERSION and
            runtime["probe"]["ruleVersion"] == RULE_VERSION and set(runtime["probe"]["members"]) == set(METRICS), "runtime identity mismatch")
    require(tree_hashes(root / "evaluator-source") == runtime["sourceFiles"], "changed evaluator source snapshot")
    require(digest(root / "submission/submission.json") == runtime["submissionSHA256"] and
            digest(root / "native-verify.log") == runtime["nativeVerificationLogSHA256"] and
            digest(root / "evaluator-probe.log") == runtime["evaluatorProbeLogSHA256"], "changed run inputs or native verification log")
    with tempfile.TemporaryDirectory(prefix="numivivo-arc2026-verify-") as temp:
        snapshot_submission(root / "submission", Path(temp) / "submission")
    query = read_json(root / "submission/query/query.json")
    expected = checked_ids([c["id"] for c in query["contexts"]])
    completed = {}
    for cid in expected:
        work = root / cid
        if not (work / "receipt.json").exists():
            continue
        r = read_json(work / "receipt.json")
        require(r["schemaVersion"] == 1 and r["id"] == cid, "context receipt mismatch")
        require(digest(work / "scores.csv") == r["scoresSHA256"] and tree_hashes(work / "bundle") == r["bundleFiles"] and
                tree_hashes(work / "metrics") == r["metricsFiles"] and digest(work / "reference.h5ad") == r["referenceSHA256"] and
                digest(root / "runtime.json") == r["runtimeSHA256"] and
                digest(root / "submission/contexts" / cid / "prediction.h5ad") == r["predictionSHA256"], "changed scoring artifacts")
        require(set(r["logs"]) == {"reference-axis.log", "metrics.log", "score.log"} and
                all(digest(work / name) == sha for name, sha in r["logs"].items()), "changed scoring logs")
        values = read_scores(work / "scores.csv")
        require(values == r["scores"], "changed score receipt")
        completed[cid] = values
    actual = aggregate(expected, completed)
    stored = read_json(root / "summary.json")
    require(stored["phase"] == query["phase"] and
            all(stored[k] == actual[k] for k in ("status", "score", "missingContexts", "contextScores")), "changed aggregate")
    print(json.dumps(actual, indent=2, allow_nan=False))
    return 0 if actual["status"] == "complete" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("run")
    for name in ("submission", "references", "evaluator", "native", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--timeout-seconds", type=int, default=86400)
    p = commands.add_parser("verify", help="verify retained hashes/aggregation; does not rerun metrics")
    p.add_argument("directory", type=Path)
    args = parser.parse_args()
    try:
        if args.command == "run":
            require(0 < args.timeout_seconds <= 604800, "timeout must be between 1 second and 7 days")
            return run(args)
        return verify(args)
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        print(f"Arc 2026 evaluation rejected: {exc}", file=sys.stderr)
        return 65


if __name__ == "__main__":
    sys.exit(main())
