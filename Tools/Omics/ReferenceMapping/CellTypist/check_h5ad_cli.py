"""Real retained Parse counts -> AnnData CSR/CSC -> native frozen annotation."""
import argparse
import gzip
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

p = argparse.ArgumentParser()
p.add_argument("cli", type=Path)
p.add_argument("evidence", type=Path)
p.add_argument("output", type=Path)
a = p.parse_args()
a.output.mkdir(exist_ok=False)
out = a.output.resolve()
evidence = a.evidence.resolve()
plan = json.loads(gzip.decompress(Path(__file__).with_name("stream-plan.json.gz").read_bytes()))
metadata = plan["metadata"]
data = json.loads((evidence / "real-input.json").read_text())
reference = json.loads((evidence / "real-output.json").read_text())
indptr = np.concatenate(([0], np.cumsum([len(r["counts"]) for r in data["rows"]])))
matrix = sparse.csr_matrix((np.array([v for r in data["rows"] for v in r["counts"]], dtype=np.int64),
                            np.array([j for r in data["rows"] for j in r["indices"]]), indptr),
                           shape=(len(data["rows"]), len(data["features"])))
assert data["features"] == [f["id"] for f in metadata["features"]]
mapping = {k: metadata[k] for k in ["id", "evidence", "sourceDescription", "countUnit", "samples"]}
mapping.update(schemaVersion=1, matrixPath="X", sampleColumn="sample", mitochondrialFeatureIDs=[])
(out / "mapping.json").write_text(json.dumps(mapping))
obs = pd.DataFrame({"sample": [c["sampleID"] for c in metadata["cells"]]},
                   index=pd.Index([c["barcode"] for c in metadata["cells"]], dtype=str))
var = pd.DataFrame(index=pd.Index(data["features"], dtype=str))
checks = []
for layout in ["csr", "csc"]:
    source = out / (layout + ".h5ad")
    obj = ad.AnnData(X=matrix if layout == "csr" else matrix.tocsc(), obs=obs.copy(), var=var.copy())
    obj.uns["provenance"] = "Retained real Parse run 1; source rows 42:309, unchanged counts and identities"
    obj.write_h5ad(source, compression="gzip")
    bundle = out / (layout + "-bundle")
    cmd = [str(a.cli.resolve()), "singlecell-h5ad-celltypist", str(source), "--plan", str(out / "mapping.json"),
           "--model", str(evidence / "model.json"), "--output", str(bundle)]
    run = subprocess.run(cmd, capture_output=True)
    (out / (layout + ".log")).write_bytes(run.stdout + run.stderr)
    run.check_returncode()
    receipt = json.loads(run.stdout)
    actual = [json.loads(line) for line in (bundle / "results.jsonl").read_text().splitlines()]
    assert len(actual) == len(reference) == 267
    error = 0.0
    for i, (x, y) in enumerate(zip(actual, reference)):
        assert x["sourceRow"] == i and x["label"] == y["label"]
        error = max(error, max(abs(u-v) for u, v in zip(x["decisions"], y["decisions"])))
        np.testing.assert_allclose(x["probabilities"], y["probabilities"], rtol=0, atol=1e-12)
    assert error < 1e-10
    assert (bundle / "original.h5ad").read_bytes() == source.read_bytes()
    retained = json.loads((bundle / "metadata.json").read_text())
    assert retained["cells"] == metadata["cells"]
    assert [f["id"] for f in retained["features"]] == data["features"]
    assert hashlib.sha256((bundle / "results.jsonl").read_bytes()).hexdigest() == receipt["resultsSHA256"]
    for filename, field in [("original.h5ad", "source"), ("mapping.json", "mapping"), ("metadata.json", "metadata")]:
        assert list(hashlib.sha256((bundle / filename).read_bytes()).digest()) == receipt[field]["bytes"]
    assert hashlib.sha256((bundle / "model.json").read_bytes()).hexdigest() == receipt["modelSHA256"]
    assert receipt["records"] == matrix.nnz
    checks.append(dict(layout=layout, cells=len(actual), records=matrix.nnz, maximumDecisionError=error,
                       sourceSHA256=hashlib.sha256(source.read_bytes()).hexdigest(), receipt=receipt))
failures = []
for case in ["negative-count", "missing-feature", "existing-output"]:
    bad = obj.copy()
    destination = out / case
    if case == "negative-count":
        bad.X.data[0] = -1
    elif case == "missing-feature":
        model = json.loads((evidence / "model.json").read_text())
        bad = bad[:, bad.var_names != model["features"][0]].copy()
    else:
        destination = out / "csc-bundle"
    source = out / (case + ".h5ad")
    bad.write_h5ad(source, compression="gzip")
    command = cmd.copy()
    command[2] = str(source)
    command[-1] = str(destination)
    failed = subprocess.run(command, capture_output=True)
    assert failed.returncode == 65, failed.stdout
    assert case == "existing-output" or not destination.exists()
    assert not list(out.glob(".numivivo-celltypist-*"))
    failures.append(dict(case=case, error=failed.stderr.decode()))
assert hashlib.sha256((out / "csc-bundle/results.jsonl").read_bytes()).hexdigest() == checks[-1]["receipt"]["resultsSHA256"]
report = dict(status="PASS", checks=checks, failures=failures, biologicalAccuracyEstablished=False,
              versions={name: importlib.metadata.version(name) for name in ["anndata", "numpy", "scipy", "h5py"]})
(out / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
