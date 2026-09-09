#!/usr/bin/env python3
"""Managed RegTools + pVACsplice runner for NumiVivo public-reference research.

Consumes the exact splice-job.json emitted by the native CLI after the researcher
fills explicit local paths, SHA-256 identities, versions and predictors. It does
not call variants, align RNA, infer HLA, install predictors, or authorize treatment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from typing import Any

MAX_DOC = 32 * 1024 * 1024
EXPECTED_INPUTS = ["annotatedVCF", "rnaBAM", "rnaBAI", "referenceFASTA", "referenceFAI", "annotationGTF"]
ALLOWED_PREDICTORS = {"NetMHC", "NetMHCpan", "NetMHCcons", "PickPocket", "SMM", "SMMPMBEC"}


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key")
        result[key] = value
    return result


def loads(data: bytes) -> Any:
    return json.loads(data, object_pairs_hook=no_duplicates,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError("Nonfinite JSON number")))


def read(path: Path, limit: int = MAX_DOC) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise ValueError("Expected a regular nonsymlink file")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise ValueError("Input exceeds its size bound")
    return data


def file_hash(path: Path) -> str:
    if path.is_symlink() or not path.is_file():
        raise ValueError("Expected regular nonsymlink file")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_new(path: Path, data: bytes) -> None:
    with path.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def is_digest(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[a-f0-9]{64}", value) is not None


def checked_record(record: dict[str, Any]) -> Path:
    if set(record) != {"path", "sha256"} or not is_digest(record["sha256"]):
        raise ValueError("File records require path and SHA-256")
    original = Path(record["path"])
    if original.is_symlink():
        raise ValueError("Input symlinks are not accepted")
    path = original.resolve(strict=True)
    if file_hash(path) != record["sha256"]:
        raise ValueError("Input digest mismatch")
    return path


def checked_executable(record: dict[str, Any]) -> Path:
    path = checked_record(record)
    if not os.access(path, os.X_OK):
        raise ValueError("Declared tool is not executable")
    return path


def inventory(directory: Path) -> tuple[str, dict[str, str]]:
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Local predictor root must be a regular directory")
    files: dict[str, str] = {}
    for path in sorted(directory.rglob("*")):
        if path.is_symlink():
            raise ValueError("Predictor-resource symlinks are not accepted")
        if path.is_file():
            files[path.relative_to(directory).as_posix()] = file_hash(path)
            if len(files) > 100_000:
                raise ValueError("Predictor inventory exceeds 100,000 files")
    if not files:
        raise ValueError("Predictor inventory is empty")
    return sha(encoded(files)), files


def validate_job(job: dict[str, Any]) -> None:
    required = {"schema", "binding", "inputs", "tools", "resourceFiles", "sampleName", "normalSampleName", "rnaStrand",
                "predictors", "epitopeLengths", "threads", "junctionMinimumReads", "variantDistance", "expressionMinimum",
                "resourceVersions", "sourceCitation"}
    if set(job) != required or job["schema"] != "numivivo.org/splice-job/v1":
        raise ValueError("Unsupported splice-job schema")
    binding = job["binding"]
    if (binding.get("dataClass") != "publicReference" or binding.get("assembly") not in {"GRCh37", "GRCh38"}
            or not is_digest(binding.get("referenceSHA256")) or not is_digest(binding.get("parentReportSHA256"))
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}", binding.get("tumorRNASampleID", ""))):
        raise ValueError("Invalid genomic case binding")
    if set(job["inputs"]) != set(EXPECTED_INPUTS) or set(job["tools"]) != {"regtools", "pvacsplice"}:
        raise ValueError("Prepared splice job requires the exact declared inputs and tools")
    if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", job["sampleName"])
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", job["normalSampleName"])
            or job["rnaStrand"] not in {"RF", "FR", "XS"}
            or type(job["threads"]) is not int or not 1 <= job["threads"] <= 16
            or type(job["junctionMinimumReads"]) is not int or not 0 <= job["junctionMinimumReads"] <= 1_000_000
            or type(job["variantDistance"]) is not int or not 0 <= job["variantDistance"] <= 1_000_000
            or not isinstance(job["expressionMinimum"], (int, float)) or not math_finite_nonnegative(job["expressionMinimum"])):
        raise ValueError("Invalid bounded splice execution settings")
    predictors = job["predictors"]
    if not isinstance(predictors, list) or not predictors or len(set(predictors)) != len(predictors) or not set(predictors) <= ALLOWED_PREDICTORS:
        raise ValueError("Only explicit supported local class-I predictors are accepted")
    lengths = job["epitopeLengths"]
    if not isinstance(lengths, list) or not lengths or len(set(lengths)) != len(lengths) or any(type(v) is not int or not 8 <= v <= 15 for v in lengths):
        raise ValueError("Invalid class-I epitope lengths")
    alleles = binding.get("hlaAlleles", [])
    if not alleles or any(re.fullmatch(r"HLA-[ABC]\*[0-9]{2,3}:[0-9]{2,3}", value) is None for value in alleles):
        raise ValueError("Case binding requires explicit class-I HLA alleles")
    versions = job["resourceVersions"]
    if (not isinstance(versions, dict) or any(not isinstance(k, str) or not isinstance(v, str) or not k or not v or v.lower() == "latest"
                                              for k, v in versions.items())
            or "regtools" not in versions or "pvacsplice" not in versions
            or any(p not in versions for p in predictors)):
        raise ValueError("Exact tool and predictor versions are required")
    if not isinstance(job["sourceCitation"], str) or not job["sourceCitation"] or len(job["sourceCitation"].encode()) > 4096:
        raise ValueError("A bounded public-reference source citation is required")
    if not isinstance(job["resourceFiles"], dict):
        raise ValueError("resourceFiles must be an explicit file-hash map")
    for path, digest in job["resourceFiles"].items():
        if not isinstance(path, str) or not path or not is_digest(digest):
            raise ValueError("resourceFiles maps paths to SHA-256 values")


def math_finite_nonnegative(value: float) -> bool:
    import math
    return math.isfinite(float(value)) and float(value) >= 0


def run_step(tool: str, argv: list[str], executable_hash: str, output_dir: Path, timeout: int) -> dict[str, Any]:
    stdout_path = output_dir / f"{tool}.stdout.log"
    stderr_path = output_dir / f"{tool}.stderr.log"
    started = time.monotonic()
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        process = subprocess.Popen(argv, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except (subprocess.TimeoutExpired, KeyboardInterrupt):
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                process.wait()
            raise
    if code != 0:
        raise RuntimeError(f"{tool} exited nonzero; inspect the bounded local logs")
    if time.monotonic() - started > timeout + 5:
        raise RuntimeError(f"{tool} exceeded its execution deadline")
    return {"tool": tool, "executableSHA256": executable_hash, "arguments": argv, "exitCode": code,
            "stdoutSHA256": file_hash(stdout_path), "stderrSHA256": file_hash(stderr_path)}


def completion(files: dict[str, bytes]) -> bytes:
    return encoded({"schema": "numivivo.org/external-capture-complete/v1",
                    "files": {name: sha(data) for name, data in sorted(files.items())}})


def execute(job_data: bytes, output: Path, timeout: int, predictor_directory: Path) -> dict[str, Any]:
    job = loads(job_data)
    validate_job(job)
    inputs = {name: checked_record(job["inputs"][name]) for name in EXPECTED_INPUTS}
    if job["inputs"]["referenceFASTA"]["sha256"] != job["binding"]["referenceSHA256"]:
        raise ValueError("Reference FASTA digest differs from the case binding")
    tools = {name: checked_executable(job["tools"][name]) for name in ["regtools", "pvacsplice"]}
    tool_hashes = {name: file_hash(path) for name, path in tools.items()}
    predictor_root = predictor_directory
    if predictor_root.is_symlink():
        raise ValueError("Predictor root symlink is not accepted")
    predictor_root = predictor_root.resolve(strict=True)
    inventory_sha, predictor_inventory = inventory(predictor_root)
    for path_text, digest in job["resourceFiles"].items():
        path = Path(path_text).resolve(strict=True)
        if file_hash(path) != digest:
            raise ValueError("Declared resource file changed")
    output.mkdir(mode=0o700)
    write_new(output / "job.json", job_data)
    write_new(output / "predictor-inventory.json", encoded(predictor_inventory))
    started_at = now()
    regtools_out = output / "regtools.raw.tsv"
    regtools_argv = [str(tools["regtools"]), "cis-splice-effects", "identify", "-o", str(regtools_out), "-s", job["rnaStrand"],
                     str(inputs["annotatedVCF"]), str(inputs["rnaBAM"]), str(inputs["referenceFASTA"]), str(inputs["annotationGTF"])]
    regtools_step = run_step("regtools", regtools_argv, tool_hashes["regtools"], output, timeout)
    if not regtools_out.is_file() or regtools_out.is_symlink():
        raise ValueError("RegTools completed without the expected TSV")
    upstream = output / "upstream"
    alleles = ",".join(job["binding"]["hlaAlleles"])
    pvac_argv = [str(tools["pvacsplice"]), "run", str(regtools_out), job["sampleName"], alleles, *job["predictors"],
                 str(upstream), str(inputs["annotatedVCF"]), str(inputs["referenceFASTA"]), str(inputs["annotationGTF"]),
                 "--iedb-install-directory", str(predictor_root), "--normal-sample-name", job["normalSampleName"],
                 "--class-i-epitope-length", ",".join(map(str, job["epitopeLengths"])), "--n-threads", str(job["threads"]),
                 "--junction-score", str(job["junctionMinimumReads"]), "--variant-distance", str(job["variantDistance"]),
                 "--expn-val", str(job["expressionMinimum"])]
    pvac_step = run_step("pvacsplice", pvac_argv, tool_hashes["pvacsplice"], output, timeout)
    reports = list((upstream / "MHC_Class_I").glob("*.all_epitopes.tsv"))
    transcripts = list(upstream.glob("*.transcripts.fa"))
    if len(reports) != 1 or len(transcripts) != 1:
        raise ValueError("Expected exactly one class-I all_epitopes report and one pVACsplice transcripts FASTA")
    report_data = read(reports[0], 16 * 1024 * 1024)
    regtools_data = read(regtools_out, 16 * 1024 * 1024)
    transcripts_data = read(transcripts[0], 16 * 1024 * 1024)
    input_hashes = {name: job["inputs"][name]["sha256"] for name in EXPECTED_INPUTS}
    input_hashes["localPredictorInventory"] = inventory_sha
    output_hashes = {"splice-report.tsv": sha(report_data), "regtools.tsv": sha(regtools_data), "transcripts.fa": sha(transcripts_data)}
    receipt = {"schema": "numivivo.org/splice-execution/v1", "binding": job["binding"], "inputSHA256": input_hashes,
               "outputSHA256": output_hashes, "steps": [regtools_step, pvac_step], "resourceVersions": job["resourceVersions"],
               "startedAt": started_at, "finishedAt": now()}
    receipt_data = encoded(receipt)
    revision = ";".join(f"{key}:{value}" for key, value in sorted(job["resourceVersions"].items()))
    if len(revision.encode()) > 256:
        revision = "resources-sha256:" + sha(encoded(job["resourceVersions"]))
    manifest = {"schema": "numivivo.org/splice-manifest/v1", "binding": job["binding"], "provenanceMode": "externalRecorded",
                "sourceCitation": job["sourceCitation"], "sourceRevision": revision, "reportStage": "allEpitopes",
                "reportSHA256": output_hashes["splice-report.tsv"], "regtoolsSHA256": output_hashes["regtools.tsv"],
                "transcriptsSHA256": output_hashes["transcripts.fa"], "executionReceiptSHA256": sha(receipt_data)}
    manifest_data = encoded(manifest)
    evidence = {"splice-manifest.json": manifest_data, "splice-report.tsv": report_data, "regtools.tsv": regtools_data,
                "transcripts.fa": transcripts_data, "splice-execution.json": receipt_data}
    for name, data in evidence.items():
        write_new(output / name, data)
    write_new(output / "complete.json", completion(evidence))  # written last
    # Detect common concurrent mutation after execution. This is not filesystem isolation.
    if any(file_hash(path) != job["inputs"][name]["sha256"] for name, path in inputs.items()):
        raise RuntimeError("An input changed during execution; discard this bundle")
    if any(file_hash(path) != tool_hashes[name] for name, path in tools.items()):
        raise RuntimeError("A tool executable changed during execution; discard this bundle")
    return {"status": "external-splice-evidence-written", "reportSHA256": output_hashes["splice-report.tsv"],
            "clinicalQualification": "none"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iedb-install-directory", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=3600)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink() or not 1 <= args.timeout_seconds <= 86400:
        parser.error("Output must be new and timeout must be 1–86400 seconds")
    result = execute(read(args.job, 128 * 1024), args.output, args.timeout_seconds, args.iedb_install_directory)
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.TimeoutExpired:
        print("pvacsplice-adapter: external step exceeded its deadline", file=sys.stderr)
        raise SystemExit(2)
    except KeyboardInterrupt:
        print("pvacsplice-adapter: cancelled", file=sys.stderr)
        raise SystemExit(2)
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as error:
        print(f"pvacsplice-adapter: {error}", file=sys.stderr)
        raise SystemExit(1)
