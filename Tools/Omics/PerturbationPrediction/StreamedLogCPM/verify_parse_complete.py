"""Read-only terminal review of Parse normalized means; no fitting or scoring."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def require(condition, message):
    if not condition:
        raise ValueError(message)

def review(study, source, binary, expected_driver):
    protocol = json.loads((study / "protocol.json").read_text())
    require(sha(study / "run.py") == expected_driver == protocol["driverSHA256"],
            "driver differs from trusted published source")
    require(sha(binary) == protocol["binarySHA256"], "native runtime differs")
    manifest_path = source / "donor-plans/manifest.json"
    require(sha(manifest_path) == protocol["manifestSHA256"], "source manifest differs")
    manifest = json.loads(manifest_path.read_text())
    for name, digest in manifest["boundSourceFiles"].items():
        require(sha(source / name) == digest, "source dependency differs: " + name)
    for name, digest in protocol["adapterSHA256"].items():
        require(sha(source / name) == digest, "adapter differs: " + name)
    require(protocol["maximumAbsoluteErrorTolerance"] == 1e-10, "tolerance differs")
    require(protocol["predictionFitted"] is False and protocol["predictionScored"] is False,
            "scope differs")
    complete = json.loads((study / "complete.json").read_text())
    require(complete["status"] == "PASS-all-donors", "not complete")
    require(complete["predictionFitted"] is False and complete["predictionScored"] is False,
            "completion scope differs")
    names = [p["donor"] for p in manifest["parts"]]
    require(len(names) == 12 and len(set(names)) == 12, "source donor coverage")
    require([p["donor"] for p in complete["donors"]] == names, "completed donor coverage")
    rows = []
    for part, receipt in zip(manifest["parts"], complete["donors"]):
        donor = part["donor"]
        require(Path(donor).name == donor and donor not in (".", ".."), "unsafe donor path")
        d = study / donor
        require(not (d / "failure.json").exists(), "retained failed donor")
        require(json.loads((d / "verification.json").read_text()) == receipt,
                "completion differs from donor receipt")
        require(receipt["status"] == "PASS", "donor did not pass")
        old = json.loads((source / "execution/ingest" / donor / "published.json").read_text())
        native_source = source / old["bundle"]
        require(sha(native_source / "receipt.json") == old["receiptSHA256"], "original receipt differs")
        planpath = source / "donor-plans" / donor / "plan.json"
        require(sha(planpath) == part["planSHA256"], "original plan differs")
        plan = json.loads(planpath.read_text())
        independent = source / old["independentCounts"]
        require(sha(independent) == old["independentCountsSHA256"], "original totals differ")
        with np.load(independent) as values:
            totals = values["cellTotals"]
        samples = {v["id"]: v["condition"] for v in plan["metadata"]["samples"]}
        groups = sorted(set(samples.values()))
        require(groups == ["IFN-beta", "PBS"], "condition identity differs")
        assignments = np.array([groups.index(samples[c["sampleID"]]) for c in plan["metadata"]["cells"]])
        expected = dict(featureIDs=[v["id"] for v in plan["metadata"]["features"]],
                        groupIDs=groups, rowGroups=assignments.tolist(), rowTotals=totals.tolist())
        p = json.loads((d / "plan.json").read_text())
        require(p == expected and sha(d / "plan.json") == receipt["planSHA256"], "normalized plan differs")
        require(len(p["featureIDs"]) == 40352, "incomplete feature axis")
        require(receipt["streamSHA256"] == old["streamSHA256"], "canonical stream differs")
        runs_path = source / "donor-plans" / donor / "runs.json"
        require(sha(runs_path) == part["runsSHA256"], "run plan differs")
        runs = json.loads(runs_path.read_text())
        for run in runs:
            name = "run-%04d.json" % run["run"]
            new_qc = json.loads((d / name).read_text())
            old_qc = json.loads((source / old["attempt"] / name).read_text())
            signature = lambda qc: [(v["start"], v["stop"], v["SHA256"]) for v in qc["ranges"]]
            require(signature(new_qc) == signature(old_qc), "source range differs")
            require(new_qc["selectedCells"] == run["selectedHi"] - run["selectedLo"], "run cell count differs")
        require(sha(d / "native.json") == receipt["nativeSHA256"], "native result differs")
        require(sha(d / "reference.npy") == receipt["referenceSHA256"], "reference differs")
        native = json.loads((d / "native.json").read_text())
        require(native["method"] == "mean-per-cell-log1p-cpm-full-axis-v1", "endpoint differs")
        require(native["featureIDs"] == p["featureIDs"] and native["groupIDs"] == groups, "result axes differ")
        require(native["cellCounts"] == np.bincount(assignments, minlength=2).tolist(), "group denominator differs")
        require(native["zeroCellCounts"] == np.bincount(assignments[totals == 0], minlength=2).tolist(), "zero cell handling differs")
        reference = np.load(d / "reference.npy", allow_pickle=False)
        actual = np.array(native["means"])
        require(actual.shape == reference.shape == (2, 40352), "mean shape differs")
        require(np.isfinite(actual).all() and np.isfinite(reference).all(), "nonfinite means")
        error = float(np.max(abs(actual - reference)))
        require(error <= 1e-10 and error == receipt["maximumAbsoluteError"], "numerical comparison differs")
        require(len(totals) == receipt["cells"] == part["cells"], "donor count differs")
        rows.append(dict(donor=donor, cells=len(totals), maximumAbsoluteError=error))
    require(sum(v["cells"] for v in rows) == complete["cells"] == 72446, "complete cohort count differs")
    return dict(status="PASS-terminal-review", donors=rows, cells=72446,
                predictionFitted=False, predictionScored=False,
                limitation="Rechecks retained source range receipts and independent reference arrays; does not refetch raw counts.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("study", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("binary", type=Path)
    parser.add_argument("--driver-sha256", required=True)
    args = parser.parse_args()
    try:
        print(json.dumps(review(args.study, args.source, args.binary, args.driver_sha256), indent=2))
    except (ValueError, KeyError, OSError, TypeError) as error:
        print(json.dumps(dict(status="NOT-VERIFIED", error=str(error))))
        raise SystemExit(1)
