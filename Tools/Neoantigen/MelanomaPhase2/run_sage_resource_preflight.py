#!/usr/bin/env python3
"""Bounded, resource-only Sage preflight for the published phase-2 catalog.

Use a known non-PXD004894, non-MSV000084787 mzML. The result is never a
discovery or held-out identification. Default mode only verifies local inputs;
--execute searches after the Mini workload owner releases the machine.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

import run_first_raw_pilot as core
import run_production_raw as production


def require(condition, message):
    if not condition:
        raise core.PilotError(message)


def execute(paths, lock, mzml, mzml_hash, dataset, work_root):
    require(sys.platform == "darwin" and platform.machine() == "arm64", "Sage preflight requires macOS arm64")
    source = lock.get("real_noncohort_sage_resource_input")
    require(isinstance(source, dict) and source.get("source_type") == "real_hla_i_noncohort" and
            source.get("dataset_accession") == dataset and source.get("sha256") == mzml_hash and
            source.get("mzml_file") == mzml.name and
            isinstance(source.get("source_record_url"), str) and
            source["source_record_url"].startswith("https://") and
            all(cohort not in source["source_record_url"] for cohort in
                ("PXD004894", "MSV000084787")) and
            isinstance(source.get("bytes"), int) and source["bytes"] >= 1024 ** 3 and
            isinstance(source.get("minimum_ms2_markers"), int) and source["minimum_ms2_markers"] >= 20000,
            "real comparable-scale noncohort Sage resource input is not pinned")
    source_commit = core.verify_source_checkout()
    lock_hash = core.sha256(production.LOCK_PATH)
    work_root, free_before, persistent_before = core.check_work_root(work_root, lock)
    sage = Path(lock["tools"]["sage"]["path"])
    tool = core.verify_tool(sage, lock["tools"]["sage"], "Sage")
    version = core.tool_version([str(sage), "--version"], lock["tools"]["sage"]["version"])
    require(production.physical_memory_bytes() >= 2 * lock["sage_maximum_rss_bytes"],
            "physical memory below resource preflight gate")
    mzml_validation = core.validate_mzml(mzml)
    require(mzml_validation["sha256"] == mzml_hash and
            mzml_validation["bytes"] == source["bytes"] and
            mzml_validation["ms2_markers"] >= source["minimum_ms2_markers"],
            "real resource input byte identity or comparable MS2 scale changed")
    name = "sage-resource-preflight-" + mzml_hash[:16]
    run_dir = work_root / name
    run_dir.mkdir(mode=0o750, exist_ok=False)
    working = run_dir / "_working"
    evidence = run_dir / "evidence"
    output = run_dir / "output"
    sage_results = working / "sage-results"
    for directory in (working, evidence, output, sage_results):
        directory.mkdir(mode=0o750)
    core.write_json(run_dir / "STATUS.started.json", {"state": "started_resource_preflight", "utc": core.utc_now()})
    try:
        config = evidence / "sage-production.json"
        core.write_durable(config, production.render_sage(paths["catalog"], mzml, sage_results))
        command = [str(sage), "--disable-telemetry-i-dont-want-to-improve-sage",
                   "--annotate-matches", "--batch-size", "1", str(config)]
        run = production.run_sage_bounded(command, sage_results, evidence / "sage.stdout",
                                          evidence / "sage.stderr-time", work_root, lock)
        result = production.validate_sage(sage_results, mzml)
        require(core.verify_source_checkout() == source_commit, "source changed during resource preflight")
        core.verify_sha(paths["catalog"], lock["catalog_fasta_sha256"], "catalog after resource preflight")
        core.verify_sha(sage, lock["tools"]["sage"]["sha256"], "Sage after resource preflight")
        compressed = {}
        for label, filename in (("sage_psm", "results.sage.tsv"),
                                ("sage_fragments", "matched_fragments.sage.tsv"),
                                ("sage_metadata", "results.json")):
            compressed[label] = production.compress_verified_bounded(sage_results / filename,
                                                                      output / (label + ".gz"), work_root, lock)
        require(production.retained_output_bytes(work_root) <= lock["maximum_persistent_output_bytes"],
                "retained phase-2 outputs exceeded 4 GiB before resource receipt")
        require(core.verify_source_checkout() == source_commit and
                core.sha256(production.LOCK_PATH) == lock_hash and
                core.sha256(production.SAGE_TEMPLATE) == lock["sage_template_sha256"],
                "frozen source/config changed before resource receipt")
        core.fsync_directory(evidence)
        core.fsync_directory(output)
        receipt = {
            "schema": "numivivo.melanoma.phase2.sage_resource_preflight.v1",
            "status": "pass_resource_only_no_candidate_inspection",
            "sealed_utc": core.utc_now(), "source_git_commit": source_commit,
            "source_dataset": dataset, "source_record_url": source["source_record_url"],
            "resource_input_mzml_sha256": mzml_hash,
            "resource_input_mzml_bytes": mzml_validation["bytes"],
            "catalog_fasta_sha256": core.sha256(paths["catalog"]),
            "sage_template_sha256": core.sha256(production.SAGE_TEMPLATE),
            "sage_binary_sha256": core.sha256(sage),
            "rendered_sage_config_sha256": core.sha256(config),
            "production_lock_sha256": core.sha256(production.LOCK_PATH),
            "runner_sha256": core.sha256(Path(__file__)),
            "memory_limit_bytes": lock["sage_maximum_rss_bytes"],
            "maximum_resident_set_size_bytes": run["maximum_resident_set_size_bytes"],
            "hardware": {"machine": platform.node(),
                         "arm64_memory_bytes": production.physical_memory_bytes()},
            "disk": {"free_before_bytes": free_before,
                     "persistent_before_bytes": persistent_before,
                     "free_before_cleanup_bytes": shutil.disk_usage(work_root).free},
            "tool_signature": tool, "tool_version": version,
            "mzml_validation": mzml_validation, "sage_validation": result,
            "sage_command": run, "compressed_outputs": compressed,
            "retained_files": core.retained_hashes(run_dir),
        }
        receipt_path = run_dir / "RECEIPT.json"
        core.write_json(receipt_path, receipt)
        receipt_hash = core.sha256(receipt_path)
        core.write_durable(run_dir / "RECEIPT.sha256", (receipt_hash + "  RECEIPT.json\n").encode())
        for item in receipt["retained_files"]:
            path = run_dir / item["path"]
            require(path.stat().st_size == item["bytes"] and core.sha256(path) == item["sha256"],
                    "resource preflight evidence changed")
        require({path.name for path in sage_results.iterdir()} == production.SAGE_RESULTS,
                "resource preflight Sage output set changed")
        for label, filename in (("sage_psm", "results.sage.tsv"),
                                ("sage_fragments", "matched_fragments.sage.tsv"),
                                ("sage_metadata", "results.json")):
            source = sage_results / filename
            require(core.sha256(source) == compressed[label]["decompressed_sha256"],
                    "resource preflight plain output changed")
            source.unlink()
        core.fsync_directory(sage_results)
        sage_results.rmdir()
        core.fsync_directory(working)
        working.rmdir()
        core.write_json(run_dir / "STATUS.complete.json", {"state": "complete_resource_preflight",
                                                          "utc": core.utc_now(),
                                                          "receipt_sha256": receipt_hash})
        for directory in (evidence, output):
            for file in directory.rglob("*"):
                if file.is_file():
                    file.chmod(0o444)
            directory.chmod(0o555)
        for file in (receipt_path, run_dir / "RECEIPT.sha256", run_dir / "STATUS.started.json",
                     run_dir / "STATUS.complete.json"):
            file.chmod(0o444)
        run_dir.chmod(0o555)
        print(json.dumps({"status": "complete_resource_preflight", "receipt_sha256": receipt_hash,
                          "run_dir": str(run_dir), "peak_rss_bytes": run["maximum_resident_set_size_bytes"]},
                         sort_keys=True))
    except Exception as error:
        failure = {"state": "failed_preserved", "utc": core.utc_now(), "error": str(error)}
        if isinstance(error, production.SageBoundedFailure):
            failure["sage_monitor"] = error.observation
        core.write_json(run_dir / "STATUS.failed.json", failure)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--mzml", type=Path, required=True, help="known noncohort tooling mzML")
    parser.add_argument("--mzml-sha256", required=True)
    parser.add_argument("--source-dataset", required=True, help="noncohort dataset provenance label")
    parser.add_argument("--work-root", type=Path, help="external resource evidence directory")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    try:
        require(re.fullmatch(r"[0-9a-f]{64}", args.mzml_sha256) is not None, "invalid pinned mzML SHA-256")
        require(re.fullmatch(r"[A-Za-z0-9._-]+", args.source_dataset) is not None and
                args.source_dataset not in ("PXD004894", "MSV000084787"),
                "resource input must be outside discovery and held-out cohorts")
        lock, paths, rows, _, _ = production.static_preflight(args.input_root, 1)
        if args.execute:
            require(isinstance(lock.get("real_noncohort_sage_resource_input"), dict),
                    "real noncohort Sage resource input is not pinned; execution blocked")
        core.exact_regular(args.mzml, "resource input mzML")
        mzml = args.mzml.resolve(strict=True)
        with paths["msv_inventory"].open(newline="") as handle:
            heldout_names = {row["mzml_file"] for row in csv.DictReader(handle, delimiter="\t")}
        require(mzml.name not in {Path(row["raw_file"]).stem + ".mzML" for row in rows} and
                mzml.name not in heldout_names,
                "resource input filename belongs to discovery or held-out cohort")
        core.verify_sha(mzml, args.mzml_sha256, "resource input mzML")
        if not args.execute:
            print(json.dumps({"status": "sage_resource_static_preflight_pass_no_search",
                              "resource_gate": "unqualified_real_noncohort_input_not_pinned" if
                              lock.get("real_noncohort_sage_resource_input") is None else "pinned_not_executed",
                              "catalog_sha256": lock["catalog_fasta_sha256"],
                              "sage_template_sha256": lock["sage_template_sha256"],
                              "resource_input_mzml_sha256": args.mzml_sha256}, sort_keys=True))
            return 0
        require(args.work_root is not None, "--execute requires --work-root")
        execute(paths, lock, mzml, args.mzml_sha256, args.source_dataset, args.work_root)
        return 0
    except (OSError, ValueError, TypeError, KeyError, core.PilotError, subprocess.TimeoutExpired) as error:
        print(f"Sage resource preflight failed closed: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
