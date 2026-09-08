#!/usr/bin/env python3
"""Check imported masses and distance constraints against serialized OpenMM inputs.

This audits model parameters only. It performs no force evaluation or dynamics,
and does not substitute input agreement for numerical or ensemble qualification.
"""
import argparse
from collections import Counter
from importlib import metadata
import hashlib
import json
import math
from pathlib import Path
import platform

from run_campaign import digest, write


def compare_parameters(native, reference):
    from openmm import unit
    particles = native["particles"]
    n = reference.getNumParticles()
    if len(particles) != n or n == 0:
        raise ValueError("particle counts differ or are empty")
    expected_masses = [float(reference.getParticleMass(i).value_in_unit(unit.dalton)) for i in range(n)]
    masses = []
    for i, particle in enumerate(particles):
        if type(particle["index"]) is not int or particle["index"] != i:
            raise ValueError("native particle indices are not canonical")
        mass = particle["massDa"]
        if type(mass) not in (int, float) or not math.isfinite(mass) or mass < 0:
            raise ValueError("invalid native mass")
        masses.append(mass)
    if any(not math.isfinite(m) or m < 0 for m in expected_masses):
        raise ValueError("invalid reference mass")
    mismatches = [i for i, (a, b) in enumerate(zip(masses, expected_masses)) if a != b]

    def edge(a, b, distance):
        if type(a) is not int or type(b) is not int or not 0 <= a < n or not 0 <= b < n or a == b:
            raise ValueError("invalid constraint particle indices")
        if type(distance) not in (int, float) or not math.isfinite(distance) or distance <= 0:
            raise ValueError("invalid constraint distance")
        return min(a, b), max(a, b), float(distance)

    actual = Counter(edge(c["a"], c["b"], c["distanceNM"]) for c in native["constraints"])
    expected = Counter()
    for i in range(reference.getNumConstraints()):
        a, b, distance = reference.getConstraintParameters(i)
        expected[edge(a, b, float(distance.value_in_unit(unit.nanometer)))] += 1
    missing, extra = expected - actual, actual - expected
    examples = lambda values: [dict(a=a, b=b, distanceNM=d, count=count) for (a,b,d),count in sorted(values.items())[:5]]
    return dict(particleCount=n, referenceConstraintCount=reference.getNumConstraints(),
        nativeConstraintCount=len(native["constraints"]), particleMassesEqual=not mismatches,
        constraintTargetsEqual=actual == expected, mismatchedMassCount=len(mismatches),
        firstMismatchedMassIndices=mismatches[:10],
        maximumMassDifferenceDa=max(abs(a-b) for a,b in zip(masses,expected_masses)),
        missingReferenceConstraintCount=sum(missing.values()), extraNativeConstraintCount=sum(extra.values()),
        missingReferenceConstraintExamples=examples(missing), extraNativeConstraintExamples=examples(extra),
        passed=not mismatches and actual == expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    import openmm as mm
    root = args.references.resolve(strict=True)
    manifest_path = root / "manifest.json"
    manifest_bytes = manifest_path.read_bytes()
    manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = json.loads(manifest_bytes)
    if manifest["schema"] != "numivivo.org/md-reference-campaign/v1":
        raise ValueError("reference manifest schema")
    names = [c["identifier"] for c in manifest["cases"]]
    if len(set(names)) != len(names) or any(not isinstance(n,str) or not n or not all(c.isalnum() or c in "-_" for c in n) for n in names):
        raise ValueError("unsafe or duplicated case identifier")
    tools = {n:digest(Path(__file__).with_name(n)) for n in ("audit_mass_constraints.py", "run_campaign.py")}
    rows, inputs = [], []
    for case in manifest["cases"]:
        name = case["identifier"]
        row = dict(identifier=name, outcome="failed")
        if case["status"] != "prepared":
            row.update(outcome="preparation-failed", error=case.get("error"))
            rows.append(row)
            continue
        try:
            request_path, xml_path = root/name/"request.json", root/name/"reference-system.xml"
            request_bytes, xml_bytes = request_path.read_bytes(), xml_path.read_bytes()
            request_hash, xml_hash = (hashlib.sha256(value).hexdigest() for value in (request_bytes, xml_bytes))
            if request_hash != case["requestSHA256"] or xml_hash != case["serializedSystemSHA256"]:
                raise ValueError("changed native request or reference model")
            inputs.extend(((request_path, request_hash), (xml_path, xml_hash)))
            request = json.loads(request_bytes)
            if request["identifier"] != name or request["referenceProvenance"]["serializedSystemSHA256"] != case["serializedSystemSHA256"]:
                raise ValueError("native provenance disagrees with the reference manifest")
            result = compare_parameters(request["system"], mm.XmlSerializer.deserialize(xml_bytes.decode("utf-8")))
            row.update(outcome="passed" if result["passed"] else "failed", comparison=result,
                requestSHA256=request_hash, serializedSystemSHA256=xml_hash)
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
        rows.append(row)
        print(json.dumps(row, allow_nan=False), flush=True)
    if digest(manifest_path) != manifest_hash or any(digest(path) != sha for path,sha in inputs) or any(digest(Path(__file__).with_name(n)) != sha for n,sha in tools.items()):
        raise ValueError("input manifest or auditor changed during execution")
    prepared = [r for r in rows if r["outcome"] != "preparation-failed"]
    result = dict(schema="numivivo.org/md-mass-constraint-audit/v1", tools=tools,
        referenceManifestSHA256=manifest_hash,
        environment=dict(python=platform.python_version(), openmm=metadata.version("openmm")),
        cases=rows, preparedPassed=bool(prepared) and all(r["outcome"] == "passed" for r in prepared),
        passed=bool(rows) and all(r["outcome"] == "passed" for r in rows),
        scope="Exact imported particle masses in dalton and unordered distance-constraint endpoints/targets in nanometers, including multiplicity, against the serialized FP64 reference parameters. No GPU arithmetic, general force-field, dynamics or ensemble claim.")
    write(args.out, result)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
