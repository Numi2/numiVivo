#!/usr/bin/env python3
"""Freeze complete external reference fits before any qPCR data are loaded."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

SITES = ("AGR", "BGI", "CNL", "COH", "MAY", "NVS")
METHODS = ("edgeR", "limma", "DESeq2")

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", type=Path, required=True)
p.add_argument("--reference-dir", default="reference-closed")
a = p.parse_args()
r = a.root
reference = r / a.reference_dir
h = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
native_freeze = json.loads((r / "native-freeze.json").read_text())
assert native_freeze["scored"] is False
assert {x["site"] for x in native_freeze["records"]} == set(SITES)
records = []
for site in SITES:
    directory = reference / site
    receipt = json.loads((directory / "receipt.json").read_text())
    assert receipt["site"] == site and receipt["features"] == 25794
    files = {}
    for method in METHODS:
        path = directory / (method + ".tsv.gz")
        raw = gzip.decompress(path.read_bytes())
        assert len(raw.splitlines()) == 25795
        files[method] = {"sha256": h(path), "bytes": path.stat().st_size,
                         "uncompressedBytes": len(raw)}
    records.append({"site": site, "directory": str(directory.relative_to(r)),
                    "receipt": receipt, "files": files})

attempts = []
for name in ("reference", "reference-recovery"):
    directory = r / name
    if not directory.exists():
        continue
    files = {}
    for path in sorted(directory.glob("*/*")):
        if path.suffix != ".gz":
            continue
        item = {"bytes": path.stat().st_size}
        try:
            raw = gzip.decompress(path.read_bytes())
            item.update(ok=True, uncompressedBytes=len(raw), rows=len(raw.splitlines()))
        except Exception as error:
            item.update(ok=False, error=type(error).__name__ + ": " + str(error))
        files[str(path.relative_to(r))] = item
    attempts.append({"directory": name, "files": files, "log": name + ".log"})

out = {"scored": False, "qPCRLoaded": False,
       "protocolSHA256": native_freeze["protocolSHA256"],
       "nativeFreezeSHA256": h(r / "native-freeze.json"),
       "referenceDirectory": str(reference.relative_to(r)),
       "records": records, "retainedEarlierAttempts": attempts,
       "scoringBoundary": "No qPCR values are loaded by this freeze."}
(r / "reference-freeze.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps({"referenceDirectory": a.reference_dir, "sites": len(records),
                  "tables": len(records) * len(METHODS),
                  "freeze": str(r / "reference-freeze.json")}))
