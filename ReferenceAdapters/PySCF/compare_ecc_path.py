#!/usr/bin/env python3
"""Independent PySCF 2.8.0 H2 path references and mandatory native comparisons.
Python is only an oracle dependency. No native energy is used to fit a reference.
"""
import argparse
import hashlib
import json
from pathlib import Path

DISTANCES = [1.2, 1.4, 1.6, 1.8, 2.0]

def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")

def export(directory):
    import numpy as np
    import pyscf
    from pyscf import fci, gto, scf
    if pyscf.__version__ != "2.8.0":
        raise RuntimeError("This fixture protocol requires PySCF 2.8.0")
    records, previous = [], None
    for i, distance in enumerate(DISTANCES):
        molecule = gto.M(atom=f"H 0 0 {-distance/2}; H 0 0 {distance/2}", unit="Bohr", basis="sto-3g", verbose=0)
        reference = scf.RHF(molecule).run(conv_tol=1e-12)
        solver = fci.FCI(reference)
        solver.conv_tol = 1e-12
        energy, _ = solver.kernel()
        if not reference.converged or not solver.converged:
            raise RuntimeError("Unconverged oracle cannot produce path fixtures")
        transport_minimum = 1.0
        if previous is not None:
            prev_mol, prev_coeff = previous
            cross = gto.intor_cross("int1e_ovlp", prev_mol, molecule)
            overlap = prev_coeff.T @ cross @ reference.mo_coeff
            # Two singleton occupied/virtual transport groups, as in the native
            # fixture. Signs of canonical orbitals need not agree across engines.
            transport_minimum = float(min(abs(overlap[0, 0]), abs(overlap[1, 1])))
        records.append({"identifier": f"h2-{i}", "coordinateBohr": distance, "rhf": float(reference.e_tot),
                        "fci": float(energy), "transportMinimumSingularValue": transport_minimum})
        previous = molecule, reference.mo_coeff.copy()
    result = {"schema": "numivivo.org/ecc-path-oracle/v1", "oracle": "PySCF", "version": pyscf.__version__,
              "basis": "STO-3G", "unit": "Bohr", "distances": DISTANCES, "points": records}
    write(directory / "reference.json", result)
    write(directory / "manifest.json", {"sha256": hashlib.sha256((directory / "reference.json").read_bytes()).hexdigest(),
                                        "schema": result["schema"], "count": len(records)})

def compare(native, directory, output):
    import math
    raw = (directory / "reference.json").read_bytes()
    manifest = json.loads((directory / "manifest.json").read_text())
    if hashlib.sha256(raw).hexdigest() != manifest["sha256"]:
        raise RuntimeError("Oracle fixture hash mismatch")
    oracle = json.loads(raw)
    result = json.loads((native / "result.json").read_text())
    request = json.loads((native / "request.json").read_text())
    if (oracle["schema"] != "numivivo.org/ecc-path-oracle/v1" or oracle["version"] != "2.8.0"
        or oracle["distances"] != DISTANCES or manifest["count"] != 5 or len(oracle["points"]) != 5
        or len(result["pointResults"]) != 5 or len(result["transportMinimumSingularValues"]) != 5
        or len(result["relativeElectronicEnergiesHartree"]) != 5 or not result["converged"]
        or request["coordinateUnit"] != "Bohr" or [p["coordinate"] for p in request["snapshots"]] != DISTANCES
        or result["pointIdentifiers"] != [p["identifier"] for p in oracle["points"]]):
        raise RuntimeError("Missing, changed, mismatched or unconverged path fixture")
    rows = []
    def check(name, actual, expected, tolerance):
        error = abs(actual - expected)
        rows.append({"observable": name, "actual": actual, "expected": expected, "absoluteError": error,
                     "tolerance": tolerance, "passed": math.isfinite(actual) and math.isfinite(expected) and error <= tolerance})
    for i, point in enumerate(oracle["points"]):
        check(f"ECC energy {i}", result["pointResults"][i]["frame"]["energyHartree"], point["fci"], 1e-8)
        check(f"physical reference energy {i}", result["pointResults"][i]["frame"]["referencePhysicalEnergyHartree"], point["fci"], 1e-8)
        check(f"relative energy {i}", result["relativeElectronicEnergiesHartree"][i], point["fci"]-oracle["points"][0]["fci"], 1e-8)
        check(f"physical transport {i}", result["transportMinimumSingularValues"][i], point["transportMinimumSingularValue"], 1e-8)
    report = {"schema": "numivivo.org/ecc-path-comparison/v1", "oracleSHA256": manifest["sha256"],
              "nativeResultSHA256": hashlib.sha256((native / "result.json").read_bytes()).hexdigest(),
              "checks": rows, "passed": all(row["passed"] for row in rows),
              "scope": "H2 finite-basis full-bath electronic stretch; not a reduced-fragment accuracy claim, Gibbs barrier or paper reaction"}
    write(output, report)
    if not report["passed"]:
        raise RuntimeError("Native path differs from independent PySCF values; see comparison report")
    print(f"PASS {len(rows)} independent ECC path comparisons")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    make = commands.add_parser("export"); make.add_argument("directory", type=Path)
    run = commands.add_parser("compare"); run.add_argument("native", type=Path); run.add_argument("reference", type=Path); run.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.command == "export": export(args.directory)
    else: compare(args.native, args.reference, args.output)
