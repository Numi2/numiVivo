#!/usr/bin/env python3
"""One ordered PXD004894 production RAW transaction; Comet and Sage, no adjudication.

Default mode verifies the frozen real inputs without network or search. --execute
requires a sealed exact-catalog Sage resource preflight and an idle Mini owner.
Failures retain the RAW, mzML, outputs, and logs; only a complete dual-engine
receipt permits deletion of exact verified reproducible working files.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path

import run_first_raw_pilot as core


HERE = Path(__file__).resolve().parent
LOCK_PATH = HERE / "PHASE2_PRODUCTION_LOCK.json"
PROTOCOL_PATH = HERE / "PHASE2_DISCOVERY_PROTOCOL.md"
COMET_TEMPLATE = HERE / "comet_txt_only.params.template"
SAGE_TEMPLATE = HERE / "sage_production.json.template"
INPUTS = {
    "pxd_inventory": "PXD004894_HLA_I_RAW_INVENTORY.tsv",
    "pxd_donors": "PXD004894_DONOR_HLA_SUMMARY.tsv",
    "msv_inventory": "MSV000084787_VALIDATION_INVENTORY.tsv",
    "catalog": "catalog-from-repo-final/combined-targets.fasta",
    "provenance": "catalog-from-repo-final/protein-provenance.tsv",
    "catalog_manifest": "catalog-from-repo-final/CATALOG_MANIFEST.json",
    "trc_bed": "Ribo-seq_ORFs.primary.bed",
    "trc_fna": "Ribo-seq_ORFs.primary.fna",
    "trc_faa": "inputs/Ribo-seq_ORFs.primary.faa",
}
SAGE_REQUIRED_COLUMNS = {
    "peptide", "proteins", "filename", "scannr", "rank", "label", "charge",
    "hyperscore", "posterior_error", "spectrum_q", "peptide_q",
}
SAGE_RESULTS = {"results.sage.tsv", "matched_fragments.sage.tsv", "results.json"}
BLOCK = 8 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise core.PilotError(message)


def run_dir_name(index, row):
    stem = Path(row["raw_file"]).stem
    require(re.fullmatch(r"[A-Za-z0-9._-]+", stem) is not None, "unsafe RAW basename")
    return f"production-{index:03d}-{row['pride_file_accession'][:16]}-{stem}"


def resolve_inputs(root):
    require(root.is_dir() and not root.is_symlink(), "input-root must be an existing real directory")
    root = root.resolve(strict=True)
    paths = {key: root / relative for key, relative in INPUTS.items()}
    require(all(root in path.parents for path in paths.values()), "input path escaped root")
    return paths


def static_preflight(input_root, run_index):
    lock = json.loads(LOCK_PATH.read_text())
    require(lock.get("schema") == "numivivo.melanoma.phase2.dual_engine_production.v1" and
            lock.get("analysis_kind") == "production_88_run", "wrong production lock")
    for path, key in ((PROTOCOL_PATH, "protocol_sha256"),
                      (COMET_TEMPLATE, "comet_template_sha256"),
                      (SAGE_TEMPLATE, "sage_template_sha256")):
        core.verify_sha(path, lock[key], key)
    paths = resolve_inputs(input_root)
    command = [sys.executable, str(HERE / "verify_phase2_protocol.py")]
    for name, path in paths.items():
        command += ["--" + name.replace("_", "-"), str(path)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    require(result.returncode == 0, "frozen nine-input preflight failed: " + result.stderr[-1500:])
    report = json.loads(result.stdout)
    require(report.get("status") == "frozen_input_preflight_pass" and
            report.get("pxd_raws") == 88 and report.get("pxd_donors") == 25 and
            report.get("heldout_runs") == 16 and
            report["sha256"]["pxd_inventory"] == lock["inventory_sha256"] and
            report["sha256"]["catalog"] == lock["catalog_fasta_sha256"] and
            report["sha256"]["catalog_manifest"] == lock["catalog_manifest_sha256"] and
            report["comet_template_sha256"] == lock["comet_template_sha256"] and
            report["sage_template_sha256"] == lock["sage_template_sha256"],
            "frozen preflight and production lock disagree")
    with paths["pxd_inventory"].open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(1 <= run_index <= 88 and len(rows) == 88, "run-index must select one of 88 frozen rows")
    row = rows[run_index - 1]
    core.check_url(row["pride_record_url"], "www.ebi.ac.uk")
    core.check_url(row["https_download_url"], "ftp.pride.ebi.ac.uk", "/" + row["raw_file"])
    return lock, paths, rows, row, report


def previous_receipt(path, run_index, rows, lock, source_commit):
    if run_index == 1:
        require(path is None, "run 1 cannot have a previous receipt")
        return None, None
    require(path is not None, "run index >1 requires --previous-run-dir")
    require(path.is_dir() and not path.is_symlink(), "previous run directory is invalid")
    receipt_path = path / "RECEIPT.json"
    hash_path = path / "RECEIPT.sha256"
    status_path = path / "STATUS.complete.json"
    for item in (receipt_path, hash_path, status_path):
        core.exact_regular(item, "previous sealed run file")
    digest = core.sha256(receipt_path)
    require(hash_path.read_text() == digest + "  RECEIPT.json\n", "previous receipt checksum changed")
    receipt = json.loads(receipt_path.read_text())
    status = json.loads(status_path.read_text())
    previous_row = rows[run_index - 2]
    require(receipt.get("schema") == "numivivo.melanoma.phase2.production_transaction.v1" and
            receipt.get("analysis_kind") == lock["analysis_kind"] and
            receipt.get("run_index") == run_index - 1 and
            receipt.get("inventory_row") == previous_row and
            receipt.get("source_git_commit") == source_commit and
            receipt.get("input_hashes", {}).get("catalog_fasta_sha256") == lock["catalog_fasta_sha256"] and
            receipt.get("input_hashes", {}).get("production_lock_sha256") == core.sha256(LOCK_PATH) and
            status.get("state") == "complete_production" and status.get("receipt_sha256") == digest,
            "previous transaction is not the exact completed predecessor")
    require(not (path / "_working").exists() and not (path / "STATUS.failed.json").exists(),
            "previous transaction has failed or unsealed working files")
    retained = receipt.get("retained_files")
    require(isinstance(retained, list), "previous transaction retained-file manifest missing")
    actual = set()
    for directory in (path / "evidence", path / "output"):
        require(directory.is_dir() and not directory.is_symlink(), "previous evidence directory invalid")
        for file in directory.rglob("*"):
            require(not file.is_symlink(), "previous evidence symlink")
            if file.is_file():
                actual.add(file.relative_to(path).as_posix())
    listed = set()
    for item in retained:
        relative = item.get("path")
        require(isinstance(relative, str) and relative.startswith(("evidence/", "output/")) and
                ".." not in Path(relative).parts and relative not in listed,
                "previous retained evidence path invalid")
        listed.add(relative)
        file = path / relative
        core.exact_regular(file, "previous retained evidence")
        require(file.stat().st_size == item.get("bytes") and core.sha256(file) == item.get("sha256"),
                "previous retained evidence checksum changed")
    require(actual == listed, "previous retained evidence file set changed")
    outputs = receipt.get("compressed_outputs", {})
    require(set(outputs) == {"comet_txt", "sage_psm", "sage_fragments", "sage_metadata"},
            "previous dual-engine outputs incomplete")
    for label, item in outputs.items():
        file = path / "output" / (label + ".gz")
        core.exact_regular(file, "previous compressed search output")
        require(file.stat().st_size == item.get("bytes") and core.sha256(file) == item.get("sha256"),
                "previous compressed search hash changed")
        content_hash = hashlib.sha256()
        with gzip.open(file, "rb") as source:
            for block in iter(lambda: source.read(BLOCK), b""):
                content_hash.update(block)
        require(content_hash.hexdigest() == item.get("decompressed_sha256"),
                "previous compressed search content changed")
    return digest, receipt.get("sage_resource_preflight_receipt_sha256")


def sage_preflight(path, lock, catalog, source_commit):
    resource_input = lock.get("real_noncohort_sage_resource_input")
    require(isinstance(resource_input, dict) and resource_input.get("source_type") == "real_hla_i_noncohort" and
            resource_input.get("dataset_accession") not in ("PXD004894", "MSV000084787") and
            isinstance(resource_input.get("sha256"), str) and
            re.fullmatch(r"[0-9a-f]{64}", resource_input["sha256"]) is not None and
            isinstance(resource_input.get("source_record_url"), str) and
            resource_input["source_record_url"].startswith("https://") and
            all(cohort not in resource_input["source_record_url"] for cohort in
                ("PXD004894", "MSV000084787")),
            "real, provenance-verified noncohort Sage resource input is not pinned")
    require(path.is_dir() and not path.is_symlink(), "Sage preflight directory is invalid")
    receipt_path = path / "RECEIPT.json"
    hash_path = path / "RECEIPT.sha256"
    complete_path = path / "STATUS.complete.json"
    for item in (receipt_path, hash_path, complete_path):
        core.exact_regular(item, "Sage preflight evidence")
    digest = core.sha256(receipt_path)
    require(hash_path.read_text() == digest + "  RECEIPT.json\n", "Sage preflight receipt checksum changed")
    receipt = json.loads(receipt_path.read_text())
    complete = json.loads(complete_path.read_text())
    require(receipt.get("schema") == "numivivo.melanoma.phase2.sage_resource_preflight.v1" and
            receipt.get("status") == "pass_resource_only_no_candidate_inspection" and
            complete.get("state") == "complete_resource_preflight" and
            complete.get("receipt_sha256") == digest and
            receipt.get("catalog_fasta_sha256") == core.sha256(catalog) == lock["catalog_fasta_sha256"] and
            receipt.get("sage_template_sha256") == lock["sage_template_sha256"] and
            receipt.get("sage_binary_sha256") == lock["tools"]["sage"]["sha256"] and
            receipt.get("resource_input_mzml_sha256") == resource_input["sha256"] and
            receipt.get("resource_input_mzml_bytes") == resource_input.get("bytes") and
            receipt.get("mzml_validation", {}).get("ms2_markers", 0) >=
            resource_input.get("minimum_ms2_markers", 10 ** 12) and
            receipt.get("source_dataset") == resource_input["dataset_accession"] and
            receipt.get("source_record_url") == resource_input["source_record_url"] and
            receipt.get("source_git_commit") == source_commit and
            receipt.get("runner_sha256") == core.sha256(HERE / "run_sage_resource_preflight.py") and
            receipt.get("production_lock_sha256") == core.sha256(LOCK_PATH) and
            receipt.get("tool_signature", {}).get("sha256") == lock["tools"]["sage"]["sha256"] and
            receipt.get("tool_signature", {}).get("cdhash") == lock["tools"]["sage"]["cdhash"] and
            receipt.get("sage_command", {}).get("exit_code") == 0 and
            receipt.get("sage_validation", {}).get("output_files") == sorted(SAGE_RESULTS) and
            receipt.get("memory_limit_bytes") == lock["sage_maximum_rss_bytes"] and
            0 < receipt.get("maximum_resident_set_size_bytes", 0) <= lock["sage_maximum_rss_bytes"] and
            receipt.get("hardware", {}).get("machine") == platform.node() and
            receipt.get("hardware", {}).get("arm64_memory_bytes") == physical_memory_bytes(),
            "Sage exact-catalog resource preflight missing, stale, or failed")
    require(not (path / "_working").exists() and not (path / "STATUS.failed.json").exists(),
            "Sage resource preflight retained failed working files")
    retained = receipt.get("retained_files")
    require(isinstance(retained, list), "Sage preflight retained-file manifest missing")
    actual = set()
    for directory in (path / "evidence", path / "output"):
        require(directory.is_dir() and not directory.is_symlink(), "Sage preflight evidence directory invalid")
        for file in directory.rglob("*"):
            require(not file.is_symlink(), "Sage preflight evidence symlink")
            if file.is_file():
                actual.add(file.relative_to(path).as_posix())
    listed = set()
    for item in retained:
        relative = item.get("path")
        require(isinstance(relative, str) and relative.startswith(("evidence/", "output/")) and
                ".." not in Path(relative).parts and relative not in listed,
                "Sage preflight retained path invalid")
        listed.add(relative)
        file = path / relative
        core.exact_regular(file, "Sage preflight retained file")
        require(file.stat().st_size == item.get("bytes") and core.sha256(file) == item.get("sha256"),
                "Sage preflight retained evidence hash changed")
    require(actual == listed, "Sage preflight retained file set changed")
    compressed = receipt.get("compressed_outputs", {})
    require(set(compressed) == {"sage_psm", "sage_fragments", "sage_metadata"},
            "Sage preflight compressed output set changed")
    for label, item in compressed.items():
        file = path / "output" / (label + ".gz")
        core.exact_regular(file, "Sage preflight compressed output")
        require(file.stat().st_size == item.get("bytes") and core.sha256(file) == item.get("sha256"),
                "Sage preflight compressed hash changed")
        digest = hashlib.sha256()
        with gzip.open(file, "rb") as source:
            for block in iter(lambda: source.read(BLOCK), b""):
                digest.update(block)
        require(digest.hexdigest() == item.get("decompressed_sha256"),
                "Sage preflight compressed content changed")
    return digest, receipt


def compress_verified_bounded(source, destination, work_root, lock):
    """Lossless gzip with a pre-write 2 GiB free-space reserve on every block."""
    core.exact_regular(source, "compression source")
    floor = lock["minimum_remaining_free_bytes"]
    digest = hashlib.sha256()
    with source.open("rb") as input_file, destination.open("xb") as output:
        with gzip.GzipFile(filename="", fileobj=output, mode="wb", compresslevel=6, mtime=0) as zipped:
            while True:
                block = input_file.read(BLOCK)
                if not block:
                    break
                # Deflate can briefly expand incompressible input. Reserve the
                # full input block plus 1 MiB before writing it.
                require(shutil.disk_usage(work_root).free >= floor + len(block) + 1024 * 1024,
                        "free space lacks the 2 GiB compression reserve")
                zipped.write(block)
                require(shutil.disk_usage(work_root).free >= floor,
                        "free space fell below 2 GiB during compression")
        output.flush(); os.fsync(output.fileno())
    require(shutil.disk_usage(work_root).free >= floor, "free space fell below 2 GiB after compression")
    with gzip.open(destination, "rb") as zipped:
        for block in iter(lambda: zipped.read(BLOCK), b""):
            digest.update(block)
    require(digest.hexdigest() == core.sha256(source), "gzip decompressed content changed")
    core.fsync_directory(destination.parent)
    return {"bytes": destination.stat().st_size, "sha256": core.sha256(destination),
            "decompressed_sha256": digest.hexdigest()}


def retained_output_bytes(work_root):
    total = 0
    for path in work_root.rglob("*"):
        if path.is_symlink():
            raise core.PilotError("symlink in phase-2 output root")
        if path.is_file() and "_working" not in path.relative_to(work_root).parts:
            total += path.stat().st_size
    return total


def physical_memory_bytes():
    result = subprocess.run(["/usr/sbin/sysctl", "-n", "hw.memsize"], capture_output=True, text=True)
    require(result.returncode == 0 and result.stdout.strip().isdigit(), "cannot read physical memory")
    return int(result.stdout.strip())


def process_group_rss_bytes(pgid):
    result = subprocess.run(["/bin/ps", "-axo", "pgid=,rss="], capture_output=True, text=True, timeout=10)
    require(result.returncode == 0, "cannot sample Sage process RSS")
    total_kib = 0
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and fields[1].isdigit() and int(fields[0]) == pgid:
            total_kib += int(fields[1])
    return total_kib * 1024


class SageBoundedFailure(core.PilotError):
    def __init__(self, error, observation):
        super().__init__(str(error))
        self.observation = observation


def run_sage_bounded(command, cwd, stdout_path, stderr_path, work_root, lock):
    started = time.monotonic()
    sampled_peak = 0
    last_sample = 0
    process = None
    phase = "open_sage_logs"
    try:
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            phase = "start_sage"
            process = subprocess.Popen(["/usr/bin/time", "-l"] + command, cwd=cwd,
                                       stdout=stdout, stderr=stderr, start_new_session=True)
            while True:
                phase = "monitor_sage_process"
                try:
                    code = process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    code = None
                phase = "sample_sage_rss"
                last_sample = process_group_rss_bytes(process.pid)
                sampled_peak = max(sampled_peak, last_sample)
                phase = "sage_rss_stop_gate"
                require(sampled_peak <= lock["sage_maximum_rss_bytes"], "Sage exceeded 12 GiB RSS stop gate")
                phase = "sage_working_bytes_stop_gate"
                require(core.working_bytes(cwd.parent) <= lock["maximum_modeled_working_bytes"],
                        "combined RAW/mzML/Comet/Sage working bytes exceeded frozen peak")
                phase = "sage_free_space_stop_gate"
                require(shutil.disk_usage(work_root).free >= lock["minimum_remaining_free_bytes"],
                        "free space fell below 2 GiB during Sage")
                if code is not None:
                    break
            stdout.flush(); os.fsync(stdout.fileno())
            stderr.flush(); os.fsync(stderr.fileno())
        phase = "sage_terminal_report"
        require(code == 0, f"Sage exited {code}; outputs retained")
        log = stderr_path.read_text(errors="replace")
        match = re.search(r"^\s*(\d+)\s+maximum resident set size\s*$", log, re.M)
        require(match is not None, "Sage maximum RSS unavailable in terminal time report")
        reported_peak = int(match.group(1))
        require(reported_peak <= lock["sage_maximum_rss_bytes"] and "finished in" in log,
                "Sage terminal completion or 12 GiB RSS gate failed")
        return {"command": command, "exit_code": code, "elapsed_seconds": round(time.monotonic() - started, 3),
                "sampled_peak_rss_bytes": sampled_peak, "maximum_resident_set_size_bytes": reported_peak,
                "stdout_sha256": core.sha256(stdout_path), "stderr_time_sha256": core.sha256(stderr_path)}
    except BaseException as error:
        try:
            free_bytes = shutil.disk_usage(work_root).free
        except OSError:
            free_bytes = None
        observation = {"runner_phase": phase, "elapsed_seconds": round(time.monotonic() - started, 3),
                       "sampled_peak_rss_bytes": sampled_peak, "last_sampled_rss_bytes": last_sample,
                       "rss_stop_gate_bytes": lock["sage_maximum_rss_bytes"],
                       "free_bytes": free_bytes,
                       "process_returncode_before_stop": process.poll() if process is not None else None,
                       "error": str(error)}
        monitor_error = None
        try:
            core.write_json(stdout_path.parent / "sage-monitor.failed.json", observation)
        except Exception as receipt_error:
            monitor_error = receipt_error
        finally:
            if process is not None and process.poll() is None:
                try:
                    os.killpg(process.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    process.wait()
        if isinstance(error, Exception):
            if monitor_error is not None:
                observation["monitor_sidecar_error"] = str(monitor_error)
            raise SageBoundedFailure(error, observation) from error
        raise


def validate_comet(path, stem, catalog, stdout_path):
    core.exact_regular(path, "production Comet TXT")
    log = stdout_path.read_text(errors="replace")
    require(log.count("Run stats:") == 1 and log.count("Search end:") == 1 and
            ("Input file: " + stem + ".mzML") in log and
            not re.search(r"(^|[^A-Za-z])(error|fatal|segmentation|abort)([^A-Za-z]|$)", log, re.I),
            "Comet terminal completion invalid")
    with path.open(newline="") as handle:
        metadata = handle.readline().rstrip("\r\n").split("\t")
        require(len(metadata) == 4 and metadata[0] == "CometVersion 2026.02 rev. 2 (6edec91)" and
                metadata[1] == stem and metadata[3] == str(catalog), "Comet TXT metadata changed")
        require(tuple(handle.readline().rstrip("\r\n").split("\t")) == core.TXT_COLUMNS,
                "Comet TXT schema changed")
        rows = 0
        for line in handle:
            require(len(line.rstrip("\r\n").split("\t")) == len(core.TXT_COLUMNS),
                    "malformed Comet TXT row")
            rows += 1
    return {"bytes": path.stat().st_size, "sha256": core.sha256(path), "rows_schema_checked": rows}


def validate_sage(directory, mzml):
    require(directory.is_dir() and not directory.is_symlink(), "Sage result directory missing")
    names = {path.name for path in directory.iterdir()}
    require(names == SAGE_RESULTS and all(path.is_file() and not path.is_symlink()
                                          for path in directory.iterdir()),
            f"Sage output set changed: {sorted(names)}")
    table = directory / "results.sage.tsv"
    with table.open(newline="") as handle:
        psm_header_line = handle.readline().rstrip("\r\n")
        header = psm_header_line.split("\t")
        require(len(header) == len(set(header)) and SAGE_REQUIRED_COLUMNS <= set(header),
                "Sage PSM schema changed")
        rows = 0
        for line in handle:
            require(len(line.rstrip("\r\n").split("\t")) == len(header), "malformed Sage PSM row")
            rows += 1
    metadata = json.loads((directory / "results.json").read_text())
    require(isinstance(metadata, dict), "Sage results.json is not an object")
    fragments = directory / "matched_fragments.sage.tsv"
    with fragments.open() as handle:
        fragment_header_line = handle.readline().rstrip("\r\n")
        require(bool(fragment_header_line), "Sage matched fragment header missing")
    return {"psm_rows_schema_checked": rows, "output_files": sorted(names),
            "psm_header_sha256": hashlib.sha256(psm_header_line.encode()).hexdigest(),
            "fragment_header_sha256": hashlib.sha256(fragment_header_line.encode()).hexdigest(),
            "metadata_sha256": core.sha256(directory / "results.json"),
            "mzml_sha256": core.sha256(mzml)}


def render_sage(catalog, mzml, output_dir):
    config = json.loads(SAGE_TEMPLATE.read_text())
    require(config["database"]["fasta"] == "__CATALOG_FASTA__" and
            config["output_directory"] == "__OUTPUT_DIRECTORY__" and
            config["mzml_paths"] == ["__MZML_PATH__"], "Sage template placeholders changed")
    config["database"]["fasta"] = str(catalog)
    config["output_directory"] = str(output_dir)
    config["mzml_paths"] = [str(mzml)]
    return (json.dumps(config, indent=2, sort_keys=True) + "\n").encode()


def execute(input_root, paths, rows, row, index, work_root, previous_dir, preflight_dir, lock):
    require(sys.platform == "darwin" and platform.machine() == "arm64", "production requires macOS arm64")
    require(isinstance(lock.get("real_noncohort_sage_resource_input"), dict),
            "real noncohort Sage resource input is not pinned; production blocked")
    source_commit = core.verify_source_checkout()
    lock_hash = core.sha256(LOCK_PATH)
    previous_hash, previous_preflight_hash = previous_receipt(previous_dir, index, rows, lock, source_commit)
    tools = {name: Path(item["path"]) for name, item in lock["tools"].items()}
    labels = {"comet": "Comet", "thermorawfileparser": "ThermoRawFileParser", "sage": "Sage"}
    tool_evidence = {name: core.verify_tool(tools[name], lock["tools"][name], labels[name])
                     for name in ("comet", "thermorawfileparser", "sage")}
    tool_versions = {
        "comet": core.tool_version([str(tools["comet"]), "--version"], lock["tools"]["comet"]["version"], True),
        "thermorawfileparser": core.tool_version([str(tools["thermorawfileparser"]), "-v"],
                                                 lock["tools"]["thermorawfileparser"]["version"]),
        "sage": core.tool_version([str(tools["sage"]), "--version"], lock["tools"]["sage"]["version"]),
    }
    preflight_hash, preflight = sage_preflight(preflight_dir, lock, paths["catalog"], source_commit)
    require(previous_preflight_hash is None or previous_preflight_hash == preflight_hash,
            "Sage resource preflight changed within the ordered 88-run campaign")
    require(physical_memory_bytes() >= 2 * lock["sage_maximum_rss_bytes"],
            "physical memory cannot safely host Sage RSS gate")
    work_root, free_before, persistent_before = core.check_work_root(work_root, lock)
    name = run_dir_name(index, row)
    run_dir = work_root / name
    run_dir.mkdir(mode=0o750, exist_ok=False)
    working = run_dir / "_working"
    evidence = run_dir / "evidence"
    output = run_dir / "output"
    for directory in (working, evidence, output):
        directory.mkdir(mode=0o750)
    core.write_json(run_dir / "STATUS.started.json", {"state": "started_production", "utc": core.utc_now(),
                                                      "run_index": index, "analysis_kind": lock["analysis_kind"]})
    try:
        core.write_json(evidence / "sage-resource-preflight-receipt.json", preflight)
        core.write_json(evidence / "selected-raw-row.json", row)
        live_record, live_headers = core.live_pride_record(row)
        core.write_json(evidence / "pride-live-record.json", live_record)
        core.write_json(evidence / "pride-live-headers.json", live_headers)
        require(shutil.disk_usage(work_root).free >= lock["minimum_free_bytes"],
                "free space fell below 7 GiB before RAW transfer")
        raw = working / row["raw_file"]
        raw_result = core.download_raw(row, raw, lock["minimum_remaining_free_bytes"])
        core.write_json(evidence / "raw-verification.json", raw_result)
        stem = Path(row["raw_file"]).stem
        mzml = working / (stem + ".mzML")
        core.verify_tool(tools["thermorawfileparser"], lock["tools"]["thermorawfileparser"], "ThermoRawFileParser")
        conversion = core.run_logged([str(tools["thermorawfileparser"]), "-i=" + str(raw), "-b=" + str(mzml),
                                      "-f=1", "-m=2", "-l=2", "-w"], working,
                                     evidence / "thermo.stdout", evidence / "thermo.stderr-time", work_root, lock)
        mzml_result = core.validate_mzml(mzml)
        core.write_json(evidence / "mzml-verification.json", mzml_result)
        core.verify_sha(paths["pxd_inventory"], lock["inventory_sha256"], "inventory before search")
        core.verify_sha(paths["catalog"], lock["catalog_fasta_sha256"], "catalog before search")
        core.verify_sha(tools["comet"], lock["tools"]["comet"]["sha256"], "Comet before search")
        core.verify_sha(tools["sage"], lock["tools"]["sage"]["sha256"], "Sage before search")
        comet_config = evidence / "comet-txt-only.params"
        core.write_durable(comet_config, COMET_TEMPLATE.read_text().replace("__CATALOG_FASTA__", str(paths["catalog"])).encode())
        comet = core.run_logged([str(tools["comet"]), "-P" + str(comet_config), str(mzml)], working,
                                evidence / "comet.stdout", evidence / "comet.stderr-time", work_root, lock)
        txt = working / (stem + ".txt")
        require(not any(path.suffix.lower() in (".pin", ".sqt", ".mzid") or path.name.endswith(".pep.xml")
                        for path in working.iterdir()), "Comet emitted forbidden redundant output")
        comet_txt = validate_comet(txt, stem, paths["catalog"], evidence / "comet.stdout")
        sage_results = working / "sage-results"
        sage_results.mkdir(mode=0o750)
        sage_config = evidence / "sage-production.json"
        core.write_durable(sage_config, render_sage(paths["catalog"], mzml, sage_results))
        core.verify_tool(tools["sage"], lock["tools"]["sage"], "Sage")
        sage_command = [str(tools["sage"]), "--disable-telemetry-i-dont-want-to-improve-sage",
                        "--annotate-matches", "--batch-size", "1", str(sage_config)]
        sage = run_sage_bounded(sage_command, sage_results, evidence / "sage.stdout",
                                evidence / "sage.stderr-time", work_root, lock)
        sage_result = validate_sage(sage_results, mzml)
        require(all(sage_result[key] == preflight["sage_validation"][key]
                    for key in ("output_files", "psm_header_sha256", "fragment_header_sha256")),
                "production Sage outputs differ from the real-resource schema preflight")
        core.write_json(evidence / "sage-validation.json", sage_result)
        require(core.sha256(LOCK_PATH) == lock_hash and core.sha256(PROTOCOL_PATH) == lock["protocol_sha256"],
                "production lock or protocol changed during search")
        for key, expected in (("pxd_inventory", "inventory_sha256"), ("catalog", "catalog_fasta_sha256")):
            core.verify_sha(paths[key], lock[expected], key + " after search")
        for name in ("comet", "thermorawfileparser", "sage"):
            core.verify_sha(tools[name], lock["tools"][name]["sha256"], name + " after search")
        require(core.verify_source_checkout() == source_commit, "production source commit changed during search")
        compressed = {}
        for label, file in (("comet_txt", txt),
                            ("sage_psm", sage_results / "results.sage.tsv"),
                            ("sage_fragments", sage_results / "matched_fragments.sage.tsv"),
                            ("sage_metadata", sage_results / "results.json")):
            destination = output / (label + ".gz")
            compressed[label] = compress_verified_bounded(file, destination, work_root, lock)
            require(compressed[label]["decompressed_sha256"] == core.sha256(file),
                    label + " gzip round-trip mismatch")
        require(shutil.disk_usage(work_root).free >= lock["minimum_remaining_free_bytes"],
                "free space fell below 2 GiB before receipt")
        require(retained_output_bytes(work_root) <= lock["maximum_persistent_output_bytes"],
                "retained phase-2 outputs exceeded 4 GiB before receipt")
        require(core.verify_source_checkout() == source_commit and core.sha256(LOCK_PATH) == lock_hash and
                core.sha256(PROTOCOL_PATH) == lock["protocol_sha256"] and
                core.sha256(COMET_TEMPLATE) == lock["comet_template_sha256"] and
                core.sha256(SAGE_TEMPLATE) == lock["sage_template_sha256"],
                "frozen source/config changed before production receipt")
        core.fsync_directory(evidence)
        core.fsync_directory(output)
        receipt = {
            "schema": "numivivo.melanoma.phase2.production_transaction.v1",
            "analysis_kind": lock["analysis_kind"], "sealed_utc": core.utc_now(),
            "run_index": index, "inventory_row": row, "previous_receipt_sha256": previous_hash,
            "sage_resource_preflight_receipt_sha256": preflight_hash,
            "source_git_commit": source_commit,
            "input_hashes": {"inventory_sha256": core.sha256(paths["pxd_inventory"]),
                             "catalog_fasta_sha256": core.sha256(paths["catalog"]),
                             "catalog_manifest_sha256": core.sha256(paths["catalog_manifest"]),
                             "production_lock_sha256": core.sha256(LOCK_PATH),
                             "protocol_sha256": core.sha256(PROTOCOL_PATH),
                             "comet_template_sha256": core.sha256(COMET_TEMPLATE),
                             "sage_template_sha256": core.sha256(SAGE_TEMPLATE),
                             "rendered_comet_config_sha256": core.sha256(comet_config),
                             "rendered_sage_config_sha256": core.sha256(sage_config),
                             "production_runner_sha256": core.sha256(Path(__file__)),
                             "pilot_helpers_sha256": core.sha256(Path(core.__file__))},
            "tool_signatures": tool_evidence, "tool_versions": tool_versions,
            "disk": {"free_before_bytes": free_before, "persistent_before_bytes": persistent_before,
                     "free_before_cleanup_bytes": shutil.disk_usage(work_root).free},
            "raw": raw_result, "mzml": mzml_result,
            "comet_txt": comet_txt, "sage_results": sage_result,
            "compressed_outputs": compressed,
            "commands": {"conversion": conversion, "comet": comet, "sage": sage},
            "retained_files": core.retained_hashes(run_dir),
        }
        receipt_path = run_dir / "RECEIPT.json"
        core.write_json(receipt_path, receipt)
        receipt_hash = core.sha256(receipt_path)
        core.write_durable(run_dir / "RECEIPT.sha256", (receipt_hash + "  RECEIPT.json\n").encode())
        require(core.sha256(receipt_path) == receipt_hash and
                core.sha256(raw) == raw_result["sha256"] and core.sha256(mzml) == mzml_result["sha256"] and
                core.sha256(txt) == comet_txt["sha256"], "sealed receipt or working files changed")
        for item in receipt["retained_files"]:
            path = run_dir / item["path"]
            require(path.stat().st_size == item["bytes"] and core.sha256(path) == item["sha256"],
                    "retained output changed before cleanup")
        for label, file in (("comet_txt", txt),
                            ("sage_psm", sage_results / "results.sage.tsv"),
                            ("sage_fragments", sage_results / "matched_fragments.sage.tsv"),
                            ("sage_metadata", sage_results / "results.json")):
            require(core.sha256(file) == compressed[label]["decompressed_sha256"],
                    label + " source changed after compression")
        expected_sage = {sage_results / name for name in SAGE_RESULTS}
        require(set(working.iterdir()) == {raw, mzml, txt, sage_results} and
                set(sage_results.iterdir()) == expected_sage, "unexpected working file before cleanup")
        for file in (raw, mzml, txt, *expected_sage):
            core.exact_regular(file, "verified reproducible cleanup input")
        for file in (raw, mzml, txt, *expected_sage):
            file.unlink()
        core.fsync_directory(sage_results)
        sage_results.rmdir()
        core.fsync_directory(working)
        working.rmdir()
        free_after = shutil.disk_usage(work_root).free
        persistent_after = sum(path.stat().st_size for path in work_root.rglob("*") if path.is_file())
        core.write_json(run_dir / "STATUS.complete.json", {"state": "complete_production", "utc": core.utc_now(),
                                                          "run_index": index, "receipt_sha256": receipt_hash,
                                                          "free_after_bytes": free_after,
                                                          "persistent_after_bytes": persistent_after,
                                                          "next_run_storage_gate": free_after >= lock["minimum_free_bytes"] and
                                                          persistent_after <= lock["maximum_persistent_output_bytes"]})
        for directory in (evidence, output):
            for file in directory.rglob("*"):
                if file.is_file():
                    file.chmod(0o444)
            directory.chmod(0o555)
        for file in (receipt_path, run_dir / "RECEIPT.sha256", run_dir / "STATUS.started.json",
                     run_dir / "STATUS.complete.json"):
            file.chmod(0o444)
        run_dir.chmod(0o555)
        print(json.dumps({"status": "complete_production", "run_index": index,
                          "run_dir": str(run_dir), "receipt_sha256": receipt_hash}, sort_keys=True))
    except Exception as error:
        failure = {"state": "failed_preserved", "utc": core.utc_now(),
                   "run_index": index, "error": str(error)}
        if isinstance(error, SageBoundedFailure):
            failure["sage_monitor"] = error.observation
        core.write_json(run_dir / "STATUS.failed.json", failure)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True, help="external verified phase-2 inputs")
    parser.add_argument("--run-index", type=int, required=True, help="one-based frozen inventory row 1..88")
    parser.add_argument("--work-root", type=Path, help="external rolling evidence directory")
    parser.add_argument("--previous-run-dir", type=Path, help="completed predecessor, local or verified archive")
    parser.add_argument("--sage-preflight-dir", type=Path, help="sealed exact-catalog resource preflight")
    parser.add_argument("--execute", action="store_true", help="perform one dual-engine production transaction")
    args = parser.parse_args()
    try:
        lock, paths, rows, row, report = static_preflight(args.input_root, args.run_index)
        if not args.execute:
            print(json.dumps({"status": "real_input_preflight_pass_no_search", "analysis_kind": lock["analysis_kind"],
                              "sage_resource_gate": "unqualified_real_noncohort_input_not_pinned" if
                              lock.get("real_noncohort_sage_resource_input") is None else "pinned_not_executed",
                              "run_index": args.run_index, "raw_file": row["raw_file"],
                              "raw_bytes": int(row["raw_bytes"]), "catalog_sha256": lock["catalog_fasta_sha256"],
                              "sage_template_sha256": lock["sage_template_sha256"],
                              "cohort_counts": {"discovery_raws": report["pxd_raws"],
                                                "discovery_donors": report["pxd_donors"],
                                                "heldout_runs": report["heldout_runs"]}}, sort_keys=True))
            return 0
        require(args.work_root is not None and args.sage_preflight_dir is not None,
                "--execute requires --work-root and --sage-preflight-dir")
        execute(args.input_root, paths, rows, row, args.run_index, args.work_root,
                args.previous_run_dir, args.sage_preflight_dir, lock)
        return 0
    except (OSError, ValueError, TypeError, KeyError, core.PilotError, subprocess.TimeoutExpired) as error:
        print(f"production failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
