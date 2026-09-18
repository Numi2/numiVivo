#!/usr/bin/env python3
"""Synthetic structural CLI campaign; geometry and both fits checked independently.
These deliberately artificial coordinates are software fixtures, not physical proteins.
"""
import argparse
import copy
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
from unittest.mock import patch

from check_native import close
from fetch_source import FEATURES, TARGETS
from run_published import execute, verify_report


def fixture(length: int, separation: float) -> dict:
    atoms, residues, coordinates = [], [], []
    for chain, count in enumerate((length, 1)):
        for residue in range(count):
            rid = len(residues)
            indices = []
            for j, (name, symbol, number) in enumerate((("N", "N", 7), ("CA", "C", 6), ("C", "C", 6), ("O", "O", 8))):
                index = len(atoms); indices.append(index)
                atoms.append({"index": index, "name": name, "element": {"symbol": symbol, "atomicNumber": number},
                              "formalCharge": 0, "residueIndex": rid, "isHetero": False})
                coordinates.append({"x": residue * 0.5 + j * 0.1, "y": separation if chain else 0, "z": 0})
            residues.append({"index": rid, "name": "GLY", "chainIndex": chain, "sequenceNumber": residue + 1, "atomIndices": indices})
    return {"schemaVersion": 1, "identifier": f"synthetic-G{length}", "atoms": atoms, "bonds": [],
            "residues": residues, "chains": [{"index": 0, "identifier": "B", "residueIndices": list(range(length))},
                                               {"index": 1, "identifier": "T", "residueIndices": [length]}],
            "conformers": [{"identifier": "fixture", "positionsNM": coordinates}], "metadata": {"evidence": "synthetic-only"}}


def geometry_reference(observation: dict, result: dict) -> None:
    structure = observation["structure"]
    points = structure["conformers"][0]["positionsNM"]
    left, right = result["binderHeavyAtomIndices"], result["targetHeavyAtomIndices"]
    pairs = [(i, j, sum((points[i][a] - points[j][a]) ** 2 for a in ("x", "y", "z"))) for i in left for j in right]
    contacts = [(i, j, d) for i, j, d in pairs if d <= 0.45 ** 2]
    assert result["evaluatedAtomPairs"] == len(pairs)
    assert result["contactAtomPairs"] == len(contacts)
    assert result["shortDistanceAtomPairs"] == sum(d < 0.2 ** 2 for _, _, d in pairs)
    assert result["binderContactAtomIndices"] == sorted({i for i, _, _ in contacts})
    assert result["targetContactAtomIndices"] == sorted({j for _, j, _ in contacts})
    close(result["minimumDistanceNM"], math.sqrt(min(d for _, _, d in pairs)))
    center = [sum(points[i][a] for i in left) / len(left) for a in ("x", "y", "z")]
    target = [sum(points[i][a] for i in right) / len(right) for a in ("x", "y", "z")]
    close(result["geometricCentroidDistanceNM"], math.dist(center, target))
    close(result["binderGeometricRadiusOfGyrationNM"], math.sqrt(sum(math.dist(center, list(points[i].values())) ** 2 for i in left) / len(left)))
    expected = {}
    for i, j, d in contacts:
        key = (structure["atoms"][i]["residueIndex"], structure["atoms"][j]["residueIndex"])
        expected[key] = min(expected.get(key, math.inf), d)
    assert len(result["residueContacts"]) == len(expected)
    for r in result["residueContacts"]:
        close(r["minimumDistanceNM"], math.sqrt(expected[(r["binderResidueIndex"], r["targetResidueIndex"])]))


def check(binary: Path, root: Path) -> None:
    root.mkdir()
    def dump(name: str, value: object) -> Path:
        path = root / name
        path.write_text(json.dumps(value, sort_keys=True, allow_nan=False))
        return path
    def run(*args: object, expected: int = 0) -> None:
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True, timeout=120)
        assert result.returncode == expected, (args, result.returncode, result.stderr)
    source = root / "source.csv"
    with source.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["uuid", "target", "sequence", "adaptyv_binding", "twist_binding", *FEATURES])
        for i in range(18):
            writer.writerow([f"s{i:02}", TARGETS[i // 6], "G" * (i + 1),
                             "binder" if i % 2 else "non_binder", "binder" if i % 3 else "non_binder",
                             (i % 5) / 5, (i % 7) / 7, (i % 3) / 3])
    config = dump("import.json", {"schemaVersion": 1, "assay": "adaptyv", "targets": TARGETS, "features": FEATURES})
    original = root / "input"
    run("binder-import", source, config, original)
    plan = {"schemaVersion": 1, "sourceSHA256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "trainingTargets": TARGETS[:2], "testTargets": [TARGETS[2]], "baselineFeature": FEATURES[0],
            "modelFeatures": FEATURES + ["numi.interface.contactAtomPairs", "numi.interface.minimumDistanceNM"],
            "topK": 3, "ridgePenalty": 0.1}
    observations = [{"candidateID": f"s{i:02}", "target": TARGETS[i // 6], "sourceLabel": "synthetic-only",
                     "structure": fixture(i + 1, (0.1, 0.4, 1.0)[i % 3]),
                     "interfacePlan": {"schemaVersion": 1, "binderChains": ["B"], "targetChains": ["T"],
                                       "conformerID": "fixture", "contactDistanceNM": 0.45,
                                       "shortDistanceNM": 0.2, "maximumPairEvaluations": 25_000_000}}
                    for i in range(17)]
    archive = dump("structures.json", {"schemaVersion": 1, "observations": observations})
    request = dump("plan.json", plan); output = root / "result"
    run("binder-evaluate-structures", original, request, archive, output)
    run("binder-verify", output)
    run("binder-evaluate-structures", expected=64)
    run("binder-evaluate-structures", original, request, archive, output, expected=65)
    geometry = json.loads((output / "geometry.json").read_text())
    enhanced = json.loads((output / "report.json").read_text())
    control = json.loads((output / "score-only-report.json").read_text())
    assert geometry["unavailableCandidateIDs"] == ["s17"]
    for observation in observations:
        geometry_reference(observation, geometry["interfaces"][observation["candidateID"]])
    verify_report(geometry, enhanced)
    matched = copy.deepcopy(geometry)
    for row in matched["dataset"]["records"]:
        if not set(plan["modelFeatures"]) <= row["features"].keys():
            row["features"] = {}
    verify_report(matched, control)
    assert enhanced["trainingIDs"] == control["trainingIDs"]
    assert enhanced["targets"][0]["candidateIDs"] == control["targets"][0]["candidateIDs"]
    receipt = json.loads((output / "receipt.json").read_text())
    assert receipt["implementationSHA256"] == hashlib.sha256(binary.read_bytes()).hexdigest()
    for name, digest in receipt["files"].items():
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == digest
    # Alter coordinates AND update their recorded hash: reconstruction still rejects.
    tampered = json.loads((output / "structures.json").read_text())
    tampered["observations"][0]["structure"]["conformers"][0]["positionsNM"][0]["z"] = 0.8
    (output / "structures.json").write_text(json.dumps(tampered))
    receipt["files"]["structures.json"] = hashlib.sha256((output / "structures.json").read_bytes()).hexdigest()
    (output / "receipt.json").write_text(json.dumps(receipt))
    run("binder-verify", output, expected=65)
    bad = copy.deepcopy(observations)
    bad[0]["candidateID"] = "unknown"
    invalid = dump("bad.json", {"schemaVersion": 1, "observations": bad})
    run("binder-evaluate-structures", original, request, invalid, root / "failed", expected=65)
    assert not (root / "failed").exists()
    # Exercise the fixed published-campaign runner on MOCKED source identity only.
    raw = source.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    with patch("run_published.BLOB", blob):
        summary = execute(root, binary, root / "mock-campaign")
        assert len(summary["folds"]) == 6 and summary["sourceRows"] == 18
        assert summary["independentNumericalVerification"] == "passed"
    try:
        execute(root, binary, root / "bad-pin")
    except ValueError:
        pass
    else:
        raise AssertionError("unpublished source was accepted without mocked identity")
    assert not (root / "bad-pin").exists()
    print("PASS: 17 synthetic structures, matched 12-training/5-test candidates, independent geometry/fit checks, tamper rejection; six mocked-source campaign folds")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="binder-geometry-") as tmp:
        check(args.binary.absolute(), Path(tmp) / "check")


if __name__ == "__main__":
    main()
