#!/usr/bin/env python3
"""One PXD004894 RAW, exploratory phase-2 search. No candidate nomination.

Default mode verifies the frozen real inventory/catalog/config without network or
search. --execute additionally performs one fail-closed RAW transaction on macOS.
Failed transactions retain every working file for inspection.
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
import ssl
import stat
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
LOCK_PATH = HERE / "FIRST_RAW_PILOT_LOCK.json"
PARAMS_PATH = HERE / "comet_txt_only.params.template"
TXT_COLUMNS = (
    "scan", "num", "charge", "exp_neutral_mass", "calc_neutral_mass",
    "e-value", "xcorr", "delta_cn", "sp_score", "ions_matched",
    "ions_total", "plain_peptide", "modified_peptide", "prev_aa",
    "next_aa", "protein", "protein_count", "modifications",
    "retention_time_sec", "sp_rank",
)
USER_AGENT = "NumiVivo-phase2-first-raw-pilot/1"
BLOCK = 8 * 1024 * 1024


class PilotError(RuntimeError):
    pass


def utc_now():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_regular(path, label):
    if path.is_symlink() or not path.is_file() or not stat.S_ISREG(path.stat().st_mode):
        raise PilotError(f"{label} must be a regular file, without a symlink: {path}")


def verify_sha(path, expected, label):
    exact_regular(path, label)
    actual = sha256(path)
    if actual != expected:
        raise PilotError(f"{label} SHA-256 mismatch: {actual}")
    return actual


def fsync_directory(path):
    descriptor = os.open(str(path), os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write_durable(path, payload):
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    fsync_directory(path.parent)


def write_json(path, value):
    write_durable(path, (json.dumps(value, indent=2, sort_keys=True) + "\n").encode())


def check_url(url, host, suffix=None):
    parsed = urllib.parse.urlsplit(url)
    if (parsed.scheme != "https" or parsed.hostname != host or parsed.port not in (None, 443)
            or parsed.username or parsed.password or parsed.query or parsed.fragment):
        raise PilotError(f"unapproved URL: {url}")
    if suffix is not None and not parsed.path.endswith(suffix):
        raise PilotError(f"URL filename differs from frozen RAW: {url}")


class PinnedRedirect(urllib.request.HTTPRedirectHandler):
    def __init__(self, host, suffix=None):
        self.host = host
        self.suffix = suffix

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        check_url(newurl, self.host, self.suffix)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def urlopen_pinned(url, host, suffix=None):
    check_url(url, host, suffix)
    opener = urllib.request.build_opener(PinnedRedirect(host, suffix), urllib.request.HTTPSHandler(context=ssl.create_default_context()))
    response = opener.open(urllib.request.Request(url, headers={"User-Agent": USER_AGENT}), timeout=120)
    check_url(response.geturl(), host, suffix)
    return response


def validate_static(inventory, catalog):
    lock = json.loads(LOCK_PATH.read_text())
    if lock.get("schema") != "numivivo.melanoma.phase2.first_raw_exploratory_pilot.v1" or lock.get("scientific_status") != "exploratory_only_no_candidate_nomination":
        raise PilotError("pilot lock is not the frozen exploratory schema")
    verify_sha(inventory, lock["inventory"]["sha256"], "PXD inventory")
    verify_sha(catalog, lock["catalog_fasta_sha256"], "published CON/UPR/TRC catalog")
    verify_sha(PARAMS_PATH, lock["comet_params_template_sha256"], "Comet parameter template")
    raw_lines = inventory.read_bytes().splitlines(keepends=True)
    if len(raw_lines) - 1 != lock["inventory"]["rows"] or len(raw_lines) < 2:
        raise PilotError("inventory does not contain exactly 88 HLA-I RAW rows")
    if hashlib.sha256(raw_lines[1]).hexdigest() != lock["inventory"]["first_row_tsv_sha256"]:
        raise PilotError("first inventory row changed")
    with inventory.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        records = list(reader)
    row = records[0]
    for key, expected in lock["first_raw"].items():
        if row.get(key) != str(expected):
            raise PilotError(f"first inventory row {key} differs from pilot lock")
    if len({item["raw_file"] for item in records}) != 88 or any(item["dataset_accession"] != "PXD004894" for item in records):
        raise PilotError("inventory has duplicate RAW files or wrong accession")
    text = PARAMS_PATH.read_text()
    if text.count("__CATALOG_FASTA__") != 1:
        raise PilotError("Comet database placeholder must occur once")
    settings = {}
    for line in text.splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, value = line.split("=", 1)
            if key.strip() in settings:
                raise PilotError(f"duplicate Comet setting: {key.strip()}")
            settings[key.strip()] = value.strip()
    required = {"output_txtfile": "1", "output_pepxmlfile": "0", "output_percolatorfile": "0", "output_sqtfile": "0", "output_mzidentmlfile": "0", "decoy_search": "1", "peptide_length_range": "8 15", "database_name": "__CATALOG_FASTA__"}
    if any(settings.get(key) != value for key, value in required.items()):
        raise PilotError("Comet template is not the frozen TXT-only search")
    if any(char.isspace() for char in str(catalog)):
        raise PilotError("catalog path contains whitespace unsupported by the frozen Comet template")
    check_url(row["pride_record_url"], "www.ebi.ac.uk")
    check_url(row["https_download_url"], "ftp.pride.ebi.ac.uk", "/" + row["raw_file"])
    return lock, row, text


def verify_tool(path, expected, label):
    verify_sha(path, expected["sha256"], label)
    result = subprocess.run(["/usr/bin/codesign", "--verify", "--strict", str(path)], capture_output=True, text=True)
    if result.returncode:
        raise PilotError(f"{label} code signature verification failed: {result.stderr}")
    detail = subprocess.run(["/usr/bin/codesign", "-dv", "--verbose=4", str(path)], capture_output=True, text=True)
    if detail.returncode or not re.search(r"^CDHash=" + re.escape(expected["cdhash"]) + r"$", detail.stderr, re.M):
        raise PilotError(f"{label} signed CDHash changed")
    if label == "ThermoRawFileParser":
        verify_sha(path.with_suffix(".dll"), expected["dll_sha256"], "ThermoRawFileParser DLL")
    return {"sha256": expected["sha256"], "cdhash": expected["cdhash"], "codesign_detail": detail.stderr}


def tool_version(command, expected, may_exit_one=False):
    result = subprocess.run(command, capture_output=True, text=True, timeout=30)
    combined = result.stdout + result.stderr
    if (result.returncode not in ((0, 1) if may_exit_one else (0,))) or expected not in combined:
        raise PilotError(f"unexpected tool version or failed launch: {command[0]}")
    return {"command": command, "exit_code": result.returncode, "stdout": result.stdout, "stderr": result.stderr}


def live_pride_record(row):
    with urlopen_pinned(row["pride_record_url"], "www.ebi.ac.uk") as response:
        payload = response.read(16 * 1024 * 1024 + 1)
        if len(payload) > 16 * 1024 * 1024:
            raise PilotError("PRIDE record exceeded 16 MiB")
        headers = dict(response.headers.items())
    record = json.loads(payload)
    if not isinstance(record, dict):
        raise PilotError("PRIDE file endpoint did not return one record")
    locations = set()
    for location in record.get("publicFileLocations", []):
        value = location.get("value", "") if isinstance(location, dict) else ""
        if value.startswith("ftp://ftp.pride.ebi.ac.uk/"):
            locations.add("https://" + value[len("ftp://"):])
    checks = {
        "project": isinstance(record.get("projectAccessions"), list) and "PXD004894" in record["projectAccessions"],
        "accession": record.get("accession") == row["pride_file_accession"],
        "category": isinstance(record.get("fileCategory"), dict) and record["fileCategory"].get("value") == "RAW",
        "filename": record.get("fileName") == row["raw_file"],
        "bytes": record.get("fileSizeBytes") == int(row["raw_bytes"]),
        "sha1": record.get("checksum") == row["sha1"],
        "location": row["https_download_url"] in locations,
    }
    if not all(checks.values()):
        raise PilotError("live PRIDE record changed: " + ", ".join(key for key, valid in checks.items() if not valid))
    return record, headers


def download_raw(row, destination, free_floor):
    partial = destination.with_name(destination.name + ".part")
    if partial.exists() or destination.exists():
        raise PilotError("RAW or partial download already exists")
    expected_bytes = int(row["raw_bytes"])
    sha1, sha256_digest, total = hashlib.sha1(), hashlib.sha256(), 0
    with urlopen_pinned(row["https_download_url"], "ftp.pride.ebi.ac.uk", "/" + row["raw_file"]) as response:
        headers = dict(response.headers.items())
        header_length = response.headers.get("Content-Length")
        if header_length is not None and int(header_length) != expected_bytes:
            raise PilotError("RAW Content-Length differs from inventory")
        with partial.open("xb") as output:
            while True:
                block = response.read(BLOCK)
                if not block:
                    break
                total += len(block)
                if total > expected_bytes:
                    raise PilotError("RAW stream exceeded pinned byte ceiling")
                output.write(block)
                sha1.update(block)
                sha256_digest.update(block)
                if shutil.disk_usage(destination.parent).free < free_floor:
                    raise PilotError("free space fell below 2 GiB during RAW download")
            output.flush()
            os.fsync(output.fileno())
    if total != expected_bytes or sha1.hexdigest() != row["sha1"]:
        raise PilotError("RAW byte count or SHA-1 differs from inventory")
    os.replace(partial, destination)
    fsync_directory(destination.parent)
    return {"bytes": total, "sha1": sha1.hexdigest(), "sha256": sha256_digest.hexdigest(), "response_headers": headers}


def working_bytes(cwd):
    paths = list(cwd.rglob("*"))
    if any(path.is_symlink() for path in paths):
        raise PilotError("tool created a symlink in the working directory")
    return sum(path.stat().st_size for path in paths if path.is_file())


def run_logged(command, cwd, stdout_path, stderr_path, work_root, lock):
    started = time.monotonic()
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        process = subprocess.Popen(["/usr/bin/time", "-l"] + command, cwd=cwd, stdout=stdout, stderr=stderr,
                                   start_new_session=True)
        try:
            while True:
                try:
                    result_code = process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    result_code = None
                if working_bytes(cwd) > lock["maximum_modeled_working_bytes"]:
                    raise PilotError("working files exceeded frozen one-file peak model")
                if shutil.disk_usage(work_root).free < lock["minimum_remaining_free_bytes"]:
                    raise PilotError("free space fell below 2 GiB during conversion or search")
                if result_code is not None:
                    break
        except BaseException:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
            raise
        stdout.flush(); os.fsync(stdout.fileno())
        stderr.flush(); os.fsync(stderr.fileno())
    if result_code:
        raise PilotError(f"command failed with exit {result_code}: {command[0]}")
    return {"command": command, "elapsed_seconds": round(time.monotonic() - started, 3), "exit_code": result_code,
            "stdout_sha256": sha256(stdout_path), "stderr_time_sha256": sha256(stderr_path)}


def validate_mzml(path):
    exact_regular(path, "converted mzML")
    spectra, ms2, root = 0, 0, None
    for event, element in ET.iterparse(path, events=("start", "end")):
        name = element.tag.rsplit("}", 1)[-1]
        if root is None and event == "start":
            root = name
        if event == "end":
            if name == "spectrum":
                spectra += 1
            if name == "cvParam" and element.attrib.get("accession") == "MS:1000511" and element.attrib.get("value") == "2":
                ms2 += 1
            element.clear()
    if root not in ("mzML", "indexedmzML") or spectra == 0 or ms2 == 0:
        raise PilotError(f"invalid mzML or missing MS2: root={root}, spectra={spectra}, ms2={ms2}")
    return {"bytes": path.stat().st_size, "sha256": sha256(path), "spectra": spectra, "ms2_markers": ms2}


def validate_txt(path, stem, catalog, comet_stdout):
    exact_regular(path, "Comet TXT")
    if path.stat().st_size == 0:
        raise PilotError("Comet TXT is empty")
    log = comet_stdout.read_text(errors="replace")
    if log.count("Run stats:") != 1 or log.count("Search end:") != 1 or ("Input file: " + stem + ".mzML") not in log:
        raise PilotError("Comet terminal completion marker is missing or duplicated")
    if re.search(r"(^|[^A-Za-z])(error|fatal|segmentation|abort)([^A-Za-z]|$)", log, re.I):
        raise PilotError("Comet terminal log contains a fatal marker")
    rows = 0
    with path.open(newline="") as handle:
        metadata = handle.readline().rstrip("\r\n").split("\t")
        if len(metadata) != 4 or metadata[0] != "CometVersion 2026.02 rev. 2 (6edec91)" or metadata[1] != stem or metadata[3] != str(catalog):
            raise PilotError("Comet TXT metadata differs from frozen run/catalog")
        if tuple(handle.readline().rstrip("\r\n").split("\t")) != TXT_COLUMNS:
            raise PilotError("Comet TXT column schema changed")
        for line in handle:
            if len(line.rstrip("\r\n").split("\t")) != len(TXT_COLUMNS):
                raise PilotError("Comet TXT contains a malformed result row")
            rows += 1
    if rows == 0:
        raise PilotError("Comet TXT contains no spectrum rows")
    return {"bytes": path.stat().st_size, "sha256": sha256(path), "rows_schema_checked": rows}


def compress_verified(source, destination):
    with source.open("rb") as input_file, destination.open("xb") as output:
        with gzip.GzipFile(filename="", fileobj=output, mode="wb", compresslevel=6, mtime=0) as zipped:
            shutil.copyfileobj(input_file, zipped, BLOCK)
        output.flush(); os.fsync(output.fileno())
    content_hash = hashlib.sha256()
    with gzip.open(destination, "rb") as zipped:
        for block in iter(lambda: zipped.read(BLOCK), b""):
            content_hash.update(block)
    if content_hash.hexdigest() != sha256(source):
        raise PilotError("gzip decompressed content differs from original TXT")
    fsync_directory(destination.parent)
    return {"bytes": destination.stat().st_size, "sha256": sha256(destination), "decompressed_sha256": content_hash.hexdigest()}


def retained_hashes(run_dir):
    result = []
    for root in (run_dir / "evidence", run_dir / "output"):
        for path in sorted(root.rglob("*")):
            if path.is_symlink():
                raise PilotError("symlink in retained evidence")
            if path.is_file():
                result.append({"path": path.relative_to(run_dir).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)})
    return result


def check_work_root(work_root, lock):
    if not work_root.is_dir() or work_root.is_symlink():
        raise PilotError("work-root must already exist as a real directory")
    work_root = work_root.resolve(strict=True)
    if work_root == REPO_ROOT or REPO_ROOT in work_root.parents or work_root in REPO_ROOT.parents:
        raise PilotError("work-root must be outside the Git checkout")
    free = shutil.disk_usage(work_root).free
    if free < lock["minimum_free_bytes"]:
        raise PilotError(f"actual free space {free} is below frozen 7 GiB gate")
    contents = list(work_root.rglob("*"))
    if any(path.is_symlink() for path in contents):
        raise PilotError("work-root contains a symlink")
    persistent = sum(path.stat().st_size for path in contents if path.is_file())
    if persistent > lock["maximum_persistent_output_bytes"]:
        raise PilotError("persistent phase-2 output exceeds 4 GiB stop gate")
    return work_root, free, persistent


def verify_source_checkout():
    command = ["git", "rev-parse", "--show-toplevel", "HEAD"]
    result = subprocess.run(command, cwd=HERE, capture_output=True, text=True)
    lines = result.stdout.splitlines()
    if result.returncode or len(lines) != 2 or Path(lines[0]).resolve() != REPO_ROOT:
        raise PilotError("runner must execute from its exact Git checkout")
    scope = "Tools/Neoantigen/MelanomaPhase2"
    status = subprocess.run(["git", "status", "--porcelain", "--", scope], cwd=REPO_ROOT,
                            capture_output=True, text=True)
    if status.returncode or status.stdout.strip():
        raise PilotError("phase-2 runner source has uncommitted changes")
    return lines[1]


def execute(inventory, catalog, work_root, lock, row, template):
    if sys.platform != "darwin" or platform.machine() != "arm64":
        raise PilotError("execution requires macOS arm64")
    source_commit = verify_source_checkout()
    comet = Path(lock["tools"]["comet"]["path"])
    parser = Path(lock["tools"]["thermorawfileparser"]["path"])
    tool_evidence = {
        "comet": verify_tool(comet, lock["tools"]["comet"], "Comet"),
        "thermorawfileparser": verify_tool(parser, lock["tools"]["thermorawfileparser"], "ThermoRawFileParser"),
    }
    tool_versions = {
        "comet": tool_version([str(comet), "--version"], lock["tools"]["comet"]["version"], True),
        "thermorawfileparser": tool_version([str(parser), "-v"], lock["tools"]["thermorawfileparser"]["version"]),
    }
    work_root, free_before, persistent_before = check_work_root(work_root, lock)
    lock_hash_at_start = sha256(LOCK_PATH)
    stem = Path(row["raw_file"]).stem
    run_dir = work_root / ("exploratory-" + row["pride_file_accession"][:16] + "-" + stem)
    run_dir.mkdir(mode=0o750, exist_ok=False)
    working = run_dir / "_working"
    evidence = run_dir / "evidence"
    output = run_dir / "output"
    for path in (working, evidence, output):
        path.mkdir(mode=0o750)
    write_json(run_dir / "STATUS.started.json", {"state": "started", "utc": utc_now(), "scientific_status": lock["scientific_status"]})
    try:
        live_record, live_headers = live_pride_record(row)
        write_json(evidence / "pride-live-record.json", live_record)
        write_json(evidence / "pride-live-headers.json", live_headers)
        write_json(evidence / "selected-raw-row.json", row)
        if shutil.disk_usage(work_root).free < lock["minimum_free_bytes"]:
            raise PilotError("free space fell below 7 GiB before download")
        raw = working / row["raw_file"]
        raw_result = download_raw(row, raw, lock["minimum_remaining_free_bytes"])
        write_json(evidence / "raw-verification.json", raw_result)
        verify_tool(parser, lock["tools"]["thermorawfileparser"], "ThermoRawFileParser")
        mzml = working / (stem + ".mzML")
        parser_command = [str(parser), "-i=" + str(raw), "-b=" + str(mzml), "-f=1", "-m=2", "-l=2", "-w"]
        parser_run = run_logged(parser_command, working, evidence / "thermo.stdout", evidence / "thermo.stderr-time", work_root, lock)
        mzml_result = validate_mzml(mzml)
        if raw.stat().st_size + mzml.stat().st_size > lock["maximum_modeled_working_bytes"]:
            raise PilotError("RAW+mzML exceeded frozen one-file peak model")
        write_json(evidence / "mzml-verification.json", mzml_result)
        verify_sha(inventory, lock["inventory"]["sha256"], "PXD inventory before search")
        verify_sha(catalog, lock["catalog_fasta_sha256"], "catalog before search")
        verify_sha(comet, lock["tools"]["comet"]["sha256"], "Comet before search")
        verify_sha(parser, lock["tools"]["thermorawfileparser"]["sha256"], "ThermoRawFileParser after conversion")
        verify_sha(parser.with_suffix(".dll"), lock["tools"]["thermorawfileparser"]["dll_sha256"], "ThermoRawFileParser DLL after conversion")
        rendered = template.replace("__CATALOG_FASTA__", str(catalog))
        params = evidence / "comet-txt-only.params"
        write_durable(params, rendered.encode())
        verify_tool(comet, lock["tools"]["comet"], "Comet")
        comet_command = [str(comet), "-P" + str(params), str(mzml)]
        comet_run = run_logged(comet_command, working, evidence / "comet.stdout", evidence / "comet.stderr-time", work_root, lock)
        verify_sha(inventory, lock["inventory"]["sha256"], "PXD inventory after search")
        verify_sha(catalog, lock["catalog_fasta_sha256"], "catalog after search")
        verify_sha(comet, lock["tools"]["comet"]["sha256"], "Comet after search")
        txt = working / (stem + ".txt")
        if any(path.suffix.lower() in (".pin", ".sqt", ".mzid") or path.name.endswith(".pep.xml") for path in working.iterdir()):
            raise PilotError("Comet emitted a forbidden redundant output")
        txt_result = validate_txt(txt, stem, catalog, evidence / "comet.stdout")
        if sum(path.stat().st_size for path in (raw, mzml, txt)) > lock["maximum_modeled_working_bytes"]:
            raise PilotError("RAW+mzML+TXT exceeded frozen one-file peak model")
        compressed = output / (stem + ".txt.gz")
        compressed_result = compress_verified(txt, compressed)
        if compressed_result["decompressed_sha256"] != txt_result["sha256"]:
            raise PilotError("compressed TXT content hash changed")
        if sha256(LOCK_PATH) != lock_hash_at_start:
            raise PilotError("pilot lock changed during transaction")
        if verify_source_checkout() != source_commit:
            raise PilotError("runner Git commit changed during transaction")
        verify_sha(PARAMS_PATH, lock["comet_params_template_sha256"], "parameter template after search")
        fsync_directory(evidence)
        fsync_directory(output)
        receipt = {
            "schema": "numivivo.melanoma.phase2.first_raw_transaction.v1",
            "scientific_status": lock["scientific_status"],
            "sealed_utc": utc_now(),
            "source_git_commit": source_commit,
            "input_hashes": {"inventory_sha256": sha256(inventory), "first_row_tsv_sha256": lock["inventory"]["first_row_tsv_sha256"],
                             "catalog_fasta_sha256": sha256(catalog), "pilot_lock_sha256": sha256(LOCK_PATH),
                             "params_template_sha256": sha256(PARAMS_PATH), "rendered_params_sha256": sha256(params),
                             "runner_sha256": sha256(Path(__file__))},
            "first_raw": lock["first_raw"], "tool_signatures": tool_evidence, "tool_versions": tool_versions,
            "disk": {"free_before_bytes": free_before, "persistent_before_bytes": persistent_before,
                     "free_before_cleanup_bytes": shutil.disk_usage(work_root).free},
            "raw": raw_result, "mzml": mzml_result, "txt": txt_result, "txt_gzip": compressed_result,
            "commands": {"conversion": parser_run, "search": comet_run},
            "retained_files": retained_hashes(run_dir),
        }
        receipt_path = run_dir / "RECEIPT.json"
        write_json(receipt_path, receipt)
        receipt_hash = sha256(receipt_path)
        write_durable(run_dir / "RECEIPT.sha256", (receipt_hash + "  RECEIPT.json\n").encode())
        if sha256(receipt_path) != receipt_hash or sha256(raw) != raw_result["sha256"] or sha256(mzml) != mzml_result["sha256"] or sha256(txt) != txt_result["sha256"] or sha256(compressed) != compressed_result["sha256"]:
            raise PilotError("sealed receipt or working files changed before cleanup")
        for item in receipt["retained_files"]:
            path = run_dir / item["path"]
            if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
                raise PilotError("retained evidence changed before cleanup")
        expected_working = {raw, mzml, txt}
        if set(working.iterdir()) != expected_working:
            raise PilotError("unexpected working file before cleanup; retaining all files")
        for path in expected_working:
            exact_regular(path, "verified cleanup input")
        # Only these exact verified, reproducible paths are removed.
        for path in (raw, mzml, txt):
            path.unlink()
        fsync_directory(working)
        working.rmdir()
        free_after = shutil.disk_usage(work_root).free
        persistent_after = sum(path.stat().st_size for path in work_root.rglob("*") if path.is_file())
        write_json(run_dir / "STATUS.complete.json", {"state": "complete_exploratory", "utc": utc_now(),
                                                       "receipt_sha256": receipt_hash,
                                                       "free_after_bytes": free_after,
                                                       "persistent_after_bytes": persistent_after,
                                                       "next_run_storage_gate": free_after >= lock["minimum_free_bytes"] and persistent_after <= lock["maximum_persistent_output_bytes"]})
        for root in (evidence, output):
            for path in root.rglob("*"):
                if path.is_file():
                    path.chmod(0o444)
            root.chmod(0o555)
        for path in (receipt_path, run_dir / "RECEIPT.sha256", run_dir / "STATUS.started.json", run_dir / "STATUS.complete.json"):
            path.chmod(0o444)
        run_dir.chmod(0o555)
        print(json.dumps({"status": "complete_exploratory", "run_dir": str(run_dir), "receipt_sha256": receipt_hash,
                          "txt_gzip_sha256": compressed_result["sha256"]}, sort_keys=True))
    except Exception as error:
        write_json(run_dir / "STATUS.failed.json", {"state": "failed_preserved", "utc": utc_now(), "error": str(error)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True, help="frozen 88-row PXD HLA-I inventory TSV")
    parser.add_argument("--catalog", type=Path, required=True, help="published combined-targets.fasta")
    parser.add_argument("--work-root", type=Path, help="existing external scratch directory; required with --execute")
    parser.add_argument("--execute", action="store_true", help="perform one exploratory RAW transaction")
    args = parser.parse_args()
    try:
        lock, row, template = validate_static(args.inventory, args.catalog)
        if not args.execute:
            print(json.dumps({"status": "real_input_preflight_pass", "scientific_status": lock["scientific_status"],
                              "raw_file": row["raw_file"], "raw_bytes": int(row["raw_bytes"]),
                              "inventory_sha256": lock["inventory"]["sha256"],
                              "catalog_fasta_sha256": lock["catalog_fasta_sha256"],
                              "params_template_sha256": lock["comet_params_template_sha256"]}, sort_keys=True))
            return 0
        if args.work_root is None:
            raise PilotError("--work-root is required with --execute")
        execute(args.inventory.resolve(strict=True), args.catalog.resolve(strict=True), args.work_root, lock, row, template)
        return 0
    except (OSError, ValueError, KeyError, PilotError, ET.ParseError, subprocess.TimeoutExpired) as error:
        print(f"pilot failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
