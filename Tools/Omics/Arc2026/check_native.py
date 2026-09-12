#!/usr/bin/env python3
"""Synthetic H5AD/CLI conformance on a built Apple native executable; no Arc outcomes.
Requires anndata/numpy/scipy in the separate reference environment.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    root = a.output.absolute(); root.mkdir(parents=False, exist_ok=False)
    native = a.native.absolute()
    genes = [f"synthetic_G{i}" for i in range(18533)]
    targets = genes[:300]
    columns = dict(perturbationColumn="perturbation")
    checked = []

    def save(path, labels, *, extra_count=False, fractional=False, swapped=False):
        rows = len(labels)
        data = np.full(rows, 10.5 if fractional else 10, dtype=np.float64)
        if extra_count: data[0] = 1000001
        matrix = sparse.csr_matrix((data, (np.arange(rows), np.zeros(rows, dtype=np.int64))), shape=(rows, len(genes)))
        var = genes.copy()
        if swapped: var[0], var[1] = var[1], var[0]
        ad.AnnData(matrix, obs=pd.DataFrame({"perturbation": labels}, index=[f"c{i}" for i in range(rows)]),
                   var=pd.DataFrame(index=var)).write_h5ad(path)

    def call(label, *args, success=True):
        with (root / (label + ".log")).open("xb") as out:
            r = subprocess.run([str(native), *map(str, args)], stdout=out, stderr=subprocess.STDOUT, check=False)
        assert (r.returncode == 0) == success, (label, r.returncode)
        checked.append({"case": label, "exit": r.returncode, "expectedSuccess": success})

    prep = dict(schemaVersion=1, phase="validation", contexts=[])
    pack = dict(schemaVersion=1, modelID="synthetic-conformance-only", modelArtifact=str(root/"model.txt"),
                trainingDataManifest=str(root/"training.json"), excludedPerturbedContexts=["A", "B", "C"], contexts=[])
    for cid in ["A", "B", "C"]:
        controls = root / f"{cid}-controls.h5ad"; target_file = root / f"{cid}-targets.json"
        save(controls, ["non-targeting"] * 2); target_file.write_text(json.dumps(targets))
        pred = root / f"{cid}-prediction.h5ad"
        # Counts differ across target groups and are deliberately not 400 cells.
        labels = ["non-targeting"] + targets + targets[:2]
        save(pred, labels)
        prep["contexts"].append(dict(id=cid, controlsH5AD=str(controls), targetsJSON=str(target_file), **columns))
        pack["contexts"].append(dict(id=cid, predictionH5AD=str(pred)))
    (root / "model.txt").write_text("Not a trained model. Synthetic raw-count validation fixture.\n")
    (root / "training.json").write_text('{"evidenceClass":"synthetic-test-only","datasets":[]}')
    (root / "prepare.json").write_text(json.dumps(prep)); (root / "pack.json").write_text(json.dumps(pack))
    call("prepare", "arc2026-prepare", root/"prepare.json", "--output", root/"query")
    call("verify-query", "arc2026-verify-query", root/"query")
    call("pack", "arc2026-pack", root/"query", "--plan", root/"pack.json", "--output", root/"submission")
    call("verify", "arc2026-verify", root/"submission")
    call("existing-output", "arc2026-pack", root/"query", "--plan", root/"pack.json", "--output", root/"submission", success=False)
    for label, opts in [("fractional",dict(fractional=True)),("over-depth",dict(extra_count=True)),("axis-order",dict(swapped=True))]:
        save(root/"A-prediction.h5ad", ["non-targeting"]+targets, **opts)
        call(label, "arc2026-pack", root/"query", "--plan", root/"pack.json", "--output", root/label, success=False)
        assert not (root/label).exists(), "failed package became public"
    save(root/"A-prediction.h5ad", ["non-targeting"]+targets[:-1])
    call("missing-target", "arc2026-pack", root/"query", "--plan", root/"pack.json", "--output", root/"missing-target", success=False)
    altered = root/"submission/contexts/A/prediction.h5ad"
    with altered.open("ab") as f: f.write(b"mutation")
    call("mutation", "arc2026-verify", root/"submission", success=False)
    result = dict(scope="synthetic-native-input-conformance-not-metric-parity-or-biological-validation",
                  nativeSHA256=hashlib.sha256(native.read_bytes()).hexdigest(), cases=checked)
    (root/"result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
