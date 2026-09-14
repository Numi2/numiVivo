#!/usr/bin/env python3
"""Verify the retained SEQC benchmark manifests and scoring outputs."""
import argparse
import gzip
import hashlib
import json

import pandas as pd

SITES = ("AGR", "BGI", "CNL", "COH", "MAY", "NVS")
METHODS = ("edgeR", "limma", "DESeq2")
SCORES = ("nativeMLE", "nativeShrinkage", "edgeR", "limma", "DESeq2")

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", type=str, required=True)
a = p.parse_args()
r = __import__("pathlib").Path(a.root)
h = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
native = json.loads((r / "native-checks.json").read_text())
assert len(native) == len(SITES)
assert all(x["features"] == 25794 and x["fits"]["finalFit"] == x["fits"]["effectShrinkageFit"] for x in native)
freeze = json.loads((r / "reference-freeze.json").read_text())
assert freeze["scored"] is False and freeze["qPCRLoaded"] is False
assert len(freeze["records"]) == 6 and {x["site"] for x in freeze["records"]} == set(SITES)
for record in freeze["records"]:
    for method in METHODS:
        path = r / record["directory"] / (method + ".tsv.gz")
        raw = gzip.decompress(path.read_bytes())
        assert len(raw.splitlines()) == 25795
        assert record["files"][method]["sha256"] == h(path)
q = pd.read_csv(r / "taqman-eligibility.tsv", sep="\t", dtype=str, keep_default_na=False)
assert q.shape[0] == 1044 and int((q["eligible"] == "True").sum()) == 785
summary = json.loads((r / "score-summary.json").read_text())
assert summary["qPCRRows"] == 1044 and len(summary["sites"]) == 6
for site in SITES:
    data = summary["sites"][site]
    scores = pd.read_csv(r / "scores" / (site + ".tsv.gz"), sep="\t", compression="gzip")
    assert len(scores) == data["eligibleReferenceFeatures"]
    assert scores["EntrezID"].is_unique
    for method in SCORES:
        item = data["originalFamily"][method]
        assert item["availableFeatures"] == item["n"]
        assert item["n"] <= item["eligibleReferenceFeatures"]
        assert item["descriptiveGate"] is True
    assert data["commonFeatureCount"] == data["commonFamily"]["nativeMLE"]["n"]
print(json.dumps({"nativeSites": len(native), "referenceTables": 18,
                  "qPCREligibleRows": 785, "scoredSites": len(summary["sites"]),
                  "status": "PASS"}))
