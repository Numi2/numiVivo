"""Read-only terminal count review. Requires the published preparation manifest."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()

def read(path):
    return json.loads(path.read_text())

def inside(root, name):
    path = (root / name).resolve()
    require(path.is_relative_to(root), "path outside study: " + name)
    return path

def review(root, manifest):
    root = root.resolve()
    require(not sys.flags.optimize, "Run without -O: retained bundle checker uses assertions")
    # Reject incomplete live runs before importing their dependencies.
    terminal = read(root / "complete.json")
    require(terminal == dict(status="passed-all-native-ingestions-and-replays",
                            cells=72446, records=124909573,
                            predictionFitted=False, predictionScored=False),
            "invalid terminal declaration")
    published = read(manifest)
    for name, digest in published["members"].items():
        require(sha(inside(root, name)) == digest, "preparation drift: " + name)
    frozen_path = root / "donor-plans/manifest.json"
    frozen = read(frozen_path)
    require((frozen["cells"], frozen["records"], frozen["features"]) ==
            (72446, 124909573, 40352), "cohort dimensions")
    parts = frozen["parts"]
    require(len(parts) == 12 and {p["donor"] for p in parts} ==
            {"Donor" + str(i) for i in range(1, 13)}, "donor roster")
    for name, digest in frozen["boundSourceFiles"].items():
        require(sha(inside(root, name)) == digest, "dependency drift: " + name)
    runtime = read(root / "runtime-freeze.json")
    require(runtime["manifestSHA256"] == sha(frozen_path) and
            runtime["binarySHA256"] == frozen["binarySHA256"] ==
            sha(root / "build-pooled/numivivo-omics"), "runtime identity")
    for name, digest in runtime["sourceFiles"].items():
        require(sha(inside(root, name)) == digest, "driver drift: " + name)
    os.environ["NUMIVIVO_PARSE_STUDY"] = str(root)
    sys.path.insert(0, str(root))
    from run_donors import check_bundle
    all_records = {}
    qc_totals = {}
    for phase in ("ingest", "verify"):
        phase_records = []
        sums = dict(selectedCells=0, selectedNonzeros=0, selectedTotalCounts=0,
                    sourceGeneCountMismatch=0, sourceTscpCountMismatch=0,
                    sourceMinusMatrixTotal=0, sourceMinusMatrixGeneCount=0)
        for part in parts:
            donor = part["donor"]
            prefix = root / "donor-plans" / donor
            for name, key in (("plan.json", "planSHA256"), ("runs.json", "runsSHA256"),
                              ("parent-selected-rows.npy", "rowMapSHA256")):
                require(sha(prefix / name) == part[key], "partition identity")
            rec = read(root / "execution" / phase / donor / "published.json")
            for key, expected in dict(status="passed", phase=phase, donor=donor,
                                      cells=part["cells"], records=part["records"],
                                      planSHA256=part["planSHA256"],
                                      binarySHA256=frozen["binarySHA256"]).items():
                require(rec[key] == expected, "donor receipt: " + key)
            attempt = inside(root, rec["attempt"])
            require(attempt.parent == root / "execution" / phase / donor,
                    "attempt belongs to wrong phase/donor")
            require(read(attempt / "execution.json") == rec, "attempt/pointer mismatch")
            bundle = inside(root, rec["bundle"])
            counts = inside(root, rec["independentCounts"])
            require(counts.parent == attempt, "independent counts ownership")
            require(sha(counts) == rec["independentCountsSHA256"] and
                    sha(bundle / "receipt.json") == rec["receiptSHA256"],
                    "output hash mismatch")
            receipt = check_bundle(bundle, prefix / "plan.json", counts)
            require(receipt["streamBytes"] == rec["streamBytes"] == part["records"] * 16
                    and bytes(receipt["stream"]["bytes"]).hex() == rec["streamSHA256"],
                    "stream identity")
            runs = read(prefix / "runs.json")
            require(len(runs) == part["runs"], "source run count")
            donor_sums = {key: 0 for key in sums}
            for run in runs:
                name = "run-%04d.json" % run["run"]
                qc = read(attempt / name)
                require(qc["selectedCells"] == run["selectedHi"] - run["selectedLo"],
                        "selected row coverage")
                require(bool(qc["ranges"]), "missing source ranges")
                if phase == "verify":
                    old = read(inside(root, all_records[donor]["attempt"]) / name)
                    # Timing/retry details may differ; source and all count QC must not.
                    for key in donor_sums:
                        require(qc[key] == old[key], "replay QC mismatch: " + key)
                    coords = lambda q: [(x["start"], x["stop"], x["SHA256"]) for x in q["ranges"]]
                    require(coords(qc) == coords(old), "source range mismatch")
                for key in donor_sums:
                    donor_sums[key] += qc[key]
            require(donor_sums["selectedCells"] == part["cells"] and
                    donor_sums["selectedNonzeros"] == part["records"] and
                    donor_sums["selectedTotalCounts"] == rec["totalCounts"],
                    "per-run totals mismatch")
            if phase == "verify":
                for key in ("bundle", "streamSHA256", "streamBytes", "totalCounts", "receiptSHA256"):
                    require(rec[key] == all_records[donor][key], "replay mismatch: " + key)
            else:
                require(bundle == attempt / "native", "ingest bundle ownership")
                all_records[donor] = rec
            for key in sums:
                sums[key] += donor_sums[key]
            phase_records.append(rec)
        require(read(root / (phase + "-complete.json")) ==
                dict(status="passed", parts=phase_records, cells=72446,
                     records=124909573, manifestSHA256=sha(frozen_path)),
                "phase completion mismatch")
        qc_totals[phase] = sums
    return dict(status="passed-terminal-count-review", donors=12, cells=72446,
                records=124909573, manifestSHA256=sha(frozen_path),
                qc=qc_totals, predictionFitted=False, predictionScored=False,
                scope="Offline receipt and independent aggregate revalidation; no new source replay")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--preparation-manifest", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = review(args.study, args.preparation_manifest)
    except (OSError, ValueError, KeyError, AssertionError, ImportError) as error:
        print(json.dumps(dict(status="not-verified", errorType=type(error).__name__,
                              error=str(error))))
        return 1
    print(json.dumps(result, indent=2))
    return 0

if __name__ == "__main__":
    sys.exit(main())

