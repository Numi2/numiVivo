#!/usr/bin/env python3
"""Verify all 88 sealed phase-2 production transactions without peptide adjudication.

This reads inventory, receipts, metadata, and compressed output byte streams to
check hashes and gzip integrity. It never parses or prints peptide result rows.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import stat
from pathlib import Path

import run_first_raw_pilot as core
import run_production_raw as production


EXPECTED_OUTPUTS = {"comet_txt", "sage_psm", "sage_fragments", "sage_metadata"}
EXPECTED_TOP = {"evidence", "output", "RECEIPT.json", "RECEIPT.sha256",
                "STATUS.started.json", "STATUS.complete.json"}


def require(condition, message):
    if not condition:
        raise core.PilotError(message)


def regular(path):
    require(not path.is_symlink() and path.is_file() and stat.S_ISREG(path.stat().st_mode),
            f"missing or nonregular sealed file: {path}")


def gzip_content_digest(path):
    digest = hashlib.sha256()
    count = 0
    regular(path)
    with gzip.open(path, "rb") as source:
        while True:
            block = source.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
            count += len(block)
    return digest.hexdigest(), count


def gzip_header_lines(path, count):
    lines = []
    with gzip.open(path, "rb") as source:
        for _ in range(count):
            line = source.readline(1024 * 1024)
            require(line.endswith(b"\n"), "compressed search header is missing or oversized")
            lines.append(line.decode("utf-8").rstrip("\r\n"))
    return lines


def verify_one(run_dir, index, row, lock, preceding, expected_commit, expected_preflight):
    require(run_dir.is_dir() and not run_dir.is_symlink(), f"run {index} directory invalid")
    require({entry.name for entry in run_dir.iterdir()} == EXPECTED_TOP,
            f"run {index} has missing or unexpected top-level files")
    for dirname in ("evidence", "output"):
        require((run_dir / dirname).is_dir() and not (run_dir / dirname).is_symlink(),
                f"run {index} {dirname} directory invalid")
    for filename in EXPECTED_TOP - {"evidence", "output"}:
        regular(run_dir / filename)
    receipt_path = run_dir / "RECEIPT.json"
    digest = core.sha256(receipt_path)
    require((run_dir / "RECEIPT.sha256").read_text() == digest + "  RECEIPT.json\n",
            f"run {index} receipt checksum mismatch")
    receipt = json.loads(receipt_path.read_text())
    started = json.loads((run_dir / "STATUS.started.json").read_text())
    complete = json.loads((run_dir / "STATUS.complete.json").read_text())
    require(receipt.get("schema") == "numivivo.melanoma.phase2.production_transaction.v1" and
            receipt.get("analysis_kind") == "production_88_run" and receipt.get("run_index") == index and
            receipt.get("inventory_row") == row and receipt.get("previous_receipt_sha256") == preceding and
            started.get("state") == "started_production" and started.get("run_index") == index and
            complete.get("state") == "complete_production" and complete.get("run_index") == index and
            complete.get("receipt_sha256") == digest and
            isinstance(complete.get("free_after_bytes"), int) and complete["free_after_bytes"] >= 0,
            f"run {index} not a completed exact inventory transaction")
    if expected_commit is not None:
        require(receipt.get("source_git_commit") == expected_commit, "production source commit changed across 88 runs")
    commit = receipt.get("source_git_commit")
    require(isinstance(commit, str) and len(commit) == 40 and
            all(c in "0123456789abcdef" for c in commit), "invalid source commit")
    hashes = receipt.get("input_hashes", {})
    for key, lock_key in (("inventory_sha256", "inventory_sha256"),
                          ("catalog_fasta_sha256", "catalog_fasta_sha256"),
                          ("catalog_manifest_sha256", "catalog_manifest_sha256"),
                          ("protocol_sha256", "protocol_sha256"),
                          ("comet_template_sha256", "comet_template_sha256"),
                          ("sage_template_sha256", "sage_template_sha256")):
        require(hashes.get(key) == lock[lock_key], f"run {index} frozen {key} mismatch")
    require(hashes.get("production_lock_sha256") == core.sha256(production.LOCK_PATH) and
            hashes.get("production_runner_sha256") and hashes.get("pilot_helpers_sha256"),
            f"run {index} production source hashes absent")
    for item in ("comet", "thermorawfileparser", "sage"):
        require(receipt.get("tool_signatures", {}).get(item, {}).get("sha256") == lock["tools"][item]["sha256"] and
                receipt["tool_signatures"][item].get("cdhash") == lock["tools"][item]["cdhash"],
                f"run {index} {item} binary identity mismatch")
    for item in ("conversion", "comet", "sage"):
        require(receipt.get("commands", {}).get(item, {}).get("exit_code") == 0,
                f"run {index} {item} not terminally successful")
    require(0 < receipt["commands"]["sage"].get("maximum_resident_set_size_bytes", 0)
            <= lock["sage_maximum_rss_bytes"], f"run {index} Sage exceeded RSS limit")
    preflight_hash = receipt.get("sage_resource_preflight_receipt_sha256")
    if expected_preflight is not None:
        require(preflight_hash == expected_preflight, "Sage resource preflight changed across runs")
    preflight_copy = run_dir / "evidence/sage-resource-preflight-receipt.json"
    regular(preflight_copy)
    require(core.sha256(preflight_copy) == preflight_hash,
            f"run {index} resource preflight receipt copy changed")
    preflight = json.loads(preflight_copy.read_text())
    source = lock.get("real_noncohort_sage_resource_input")
    require(isinstance(source, dict) and source.get("source_type") == "real_hla_i_noncohort",
            "real noncohort Sage resource input remains unqualified")
    require(preflight.get("schema") == "numivivo.melanoma.phase2.sage_resource_preflight.v1" and
            preflight.get("status") == "pass_resource_only_no_candidate_inspection" and
            preflight.get("source_git_commit") == commit and
            preflight.get("source_dataset") == source.get("dataset_accession") and
            preflight.get("source_record_url") == source.get("source_record_url") and
            preflight.get("resource_input_mzml_sha256") == source.get("sha256") and
            preflight.get("resource_input_mzml_bytes") == source.get("bytes") and
            preflight.get("mzml_validation", {}).get("ms2_markers", 0) >=
            source.get("minimum_ms2_markers", 10 ** 12) and
            preflight.get("catalog_fasta_sha256") == lock["catalog_fasta_sha256"] and
            preflight.get("sage_template_sha256") == lock["sage_template_sha256"] and
            preflight.get("sage_binary_sha256") == lock["tools"]["sage"]["sha256"] and
            0 < preflight.get("maximum_resident_set_size_bytes", 0) <= lock["sage_maximum_rss_bytes"],
            f"run {index} invalid Sage resource preflight")
    retained = receipt.get("retained_files")
    require(isinstance(retained, list), f"run {index} retained file manifest missing")
    actual = set()
    for dirname in ("evidence", "output"):
        for path in (run_dir / dirname).rglob("*"):
            require(not path.is_symlink(), f"run {index} retained symlink")
            if path.is_file():
                actual.add(path.relative_to(run_dir).as_posix())
    listed = set()
    for item in retained:
        relative = item.get("path")
        require(isinstance(relative, str) and not relative.startswith("/") and
                ".." not in Path(relative).parts and relative.startswith(("evidence/", "output/")) and
                relative not in listed, f"run {index} invalid retained path")
        listed.add(relative)
        path = run_dir / relative
        regular(path)
        require(path.stat().st_size == item.get("bytes") and core.sha256(path) == item.get("sha256"),
                f"run {index} retained file hash mismatch: {relative}")
    require(listed == actual, f"run {index} retained file set mismatch")
    require(hashes.get("rendered_comet_config_sha256") == core.sha256(run_dir / "evidence/comet-txt-only.params") and
            hashes.get("rendered_sage_config_sha256") == core.sha256(run_dir / "evidence/sage-production.json"),
            f"run {index} rendered search configuration changed")
    outputs = receipt.get("compressed_outputs", {})
    require(set(outputs) == EXPECTED_OUTPUTS and
            {path.name for path in (run_dir / "output").iterdir()} ==
            {label + ".gz" for label in EXPECTED_OUTPUTS},
            f"run {index} compressed dual-engine output set incomplete")
    for label in sorted(EXPECTED_OUTPUTS):
        path = run_dir / "output" / (label + ".gz")
        meta = outputs[label]
        require(path.stat().st_size == meta.get("bytes") and core.sha256(path) == meta.get("sha256"),
                f"run {index} {label} compressed hash mismatch")
        uncompressed_hash, uncompressed_bytes = gzip_content_digest(path)
        require(uncompressed_hash == meta.get("decompressed_sha256") and uncompressed_bytes > 0,
                f"run {index} {label} gzip round-trip mismatch")
    require(outputs["comet_txt"]["decompressed_sha256"] == receipt.get("comet_txt", {}).get("sha256") and
            outputs["sage_metadata"]["decompressed_sha256"] ==
            receipt.get("sage_results", {}).get("metadata_sha256"),
            f"run {index} compressed search content unbound")
    comet_meta, comet_header = gzip_header_lines(run_dir / "output/comet_txt.gz", 2)
    comet_fields = comet_meta.split("\t")
    stem = Path(row["raw_file"]).stem
    require(len(comet_fields) == 4 and
            comet_fields[0] == "CometVersion 2026.02 rev. 2 (6edec91)" and
            comet_fields[1] == stem and
            tuple(comet_header.split("\t")) == core.TXT_COLUMNS,
            f"run {index} Comet compressed header changed")
    sage_header = gzip_header_lines(run_dir / "output/sage_psm.gz", 1)[0]
    sage_columns = sage_header.split("\t")
    fragment_header = gzip_header_lines(run_dir / "output/sage_fragments.gz", 1)[0]
    sage_validation = receipt.get("sage_results", {})
    require(len(sage_columns) == len(set(sage_columns)) and
            production.SAGE_REQUIRED_COLUMNS <= set(sage_columns) and
            hashlib.sha256(sage_header.encode()).hexdigest() == sage_validation.get("psm_header_sha256") and
            hashlib.sha256(fragment_header.encode()).hexdigest() == sage_validation.get("fragment_header_sha256") and
            sage_validation.get("psm_header_sha256") ==
            preflight.get("sage_validation", {}).get("psm_header_sha256") and
            sage_validation.get("fragment_header_sha256") ==
            preflight.get("sage_validation", {}).get("fragment_header_sha256"),
            f"run {index} Sage compressed schema differs from resource preflight")
    require(not (run_dir / "_working").exists() and not (run_dir / "STATUS.failed.json").exists(),
            f"run {index} retains failed or unsealed working files")
    return digest, commit, preflight_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--runs-root", type=Path, action="append", required=True,
                        help="repeat for local and independently verified archive roots")
    args = parser.parse_args()
    try:
        lock, _, rows, _, _ = production.static_preflight(args.input_root, 1)
        roots = []
        for root in args.runs_root:
            require(root.is_dir() and not root.is_symlink(), "runs-root is not a real directory")
            roots.append(root.resolve(strict=True))
        require(len(set(roots)) == len(roots), "duplicate runs-root")
        found = {}
        for root in roots:
            for path in root.iterdir():
                if path.name.startswith("production-"):
                    require(path.is_dir() and not path.is_symlink(), "production path is not a real directory")
                    require(path.name not in found, f"duplicate production run: {path.name}")
                    found[path.name] = path
        expected = {production.run_dir_name(i, row) for i, row in enumerate(rows, 1)}
        require(set(found) == expected and len(found) == 88,
                f"incomplete 88-run production set: {len(found)} observed, 88 required")
        preceding = None
        commit = None
        preflight = None
        for index, row in enumerate(rows, 1):
            path = found[production.run_dir_name(index, row)]
            preceding, commit, preflight = verify_one(path, index, row, lock, preceding, commit, preflight)
        print(json.dumps({"status": "complete_88_production_receipts_no_peptide_adjudication",
                          "production_runs": 88, "discovery_donors": 25,
                          "catalog_fasta_sha256": lock["catalog_fasta_sha256"],
                          "source_git_commit": commit,
                          "sage_resource_preflight_receipt_sha256": preflight,
                          "last_receipt_sha256": preceding}, sort_keys=True))
        return 0
    except (OSError, ValueError, TypeError, KeyError, core.PilotError) as error:
        print(f"completion failed closed: {error}", file=os.sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
