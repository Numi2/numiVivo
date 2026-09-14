#!/usr/bin/env python3
"""Score frozen native and R reference effects against independent TaqMan values."""
import argparse
import collections
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

SITES = ("AGR", "BGI", "CNL", "COH", "MAY", "NVS")
METHODS = ("nativeMLE", "nativeShrinkage", "edgeR", "limma", "DESeq2")
R_METHODS = {"edgeR": ("edgeR.tsv.gz", "logFC"),
             "limma": ("limma.tsv.gz", "logFC"),
             "DESeq2": ("DESeq2.tsv.gz", "log2FoldChange")}

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", type=Path, required=True)
p.add_argument("--reference-dir", default="reference-closed")
p.add_argument("--qPCR", type=Path)
a = p.parse_args()
r = a.root
qpath = a.qPCR or (r / "taqman.tsv")
h = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
freeze = json.loads((r / "reference-freeze.json").read_text())
assert freeze["scored"] is False and freeze["qPCRLoaded"] is False
assert freeze["referenceDirectory"] == a.reference_dir
q = pd.read_csv(qpath, sep="\t", dtype=str, keep_default_na=False)
value_columns = [f"{arm}{i}_value" for arm in "ABCD" for i in range(1, 5)]
detection_columns = [f"{arm}{i}_detection" for arm in "ABCD" for i in range(1, 5)]
for column in value_columns:
    q[column] = pd.to_numeric(q[column], errors="coerce")

def missing_literal(value):
    return value in ("", "NA", "NaN", "<NA>")

q_ids = collections.defaultdict(list)
for i, value in enumerate(q["EntrezID"]):
    if not missing_literal(value):
        q_ids[value].append(i)
q_reasons = []
q_effects = []
for _, row in q.iterrows():
    reasons = []
    eid = row["EntrezID"]
    if missing_literal(eid):
        reasons.append("missingEntrezID")
    elif len(q_ids[eid]) != 1:
        reasons.append("qPCREntrezIDNotUnique")
    if any(row[c] != "P" for c in detection_columns[:8]):
        reasons.append("AorBDetectionNotAllP")
    if any(not np.isfinite(row[c]) or row[c] <= 0 for c in value_columns[:8]):
        reasons.append("AorBValueNotPositiveFinite")
    q_reasons.append(";".join(reasons) if reasons else "")
    q_effects.append(float(np.log2(row[[f"B{i}_value" for i in range(1, 5)]].to_numpy(dtype=float)).mean() -
                          np.log2(row[[f"A{i}_value" for i in range(1, 5)]].to_numpy(dtype=float)).mean())
                    if not reasons else np.nan)
q["eligible"] = np.asarray([not x for x in q_reasons])
q["referenceLog2FoldChange"] = q_effects
q["exclusionReason"] = q_reasons
q.to_csv(r / "taqman-eligibility.tsv", sep="\t", index=False, na_rep="NA")
eligible_q = {q.loc[i, "EntrezID"]: float(q.loc[i, "referenceLog2FoldChange"])
              for i in range(len(q)) if bool(q.loc[i, "eligible"])}

def finite(value):
    return value is not None and np.isfinite(float(value))

def metrics(reference, estimate, denominator):
    ref = np.asarray(reference, dtype=float)
    est = np.asarray(estimate, dtype=float)
    error = est - ref
    result = {"n": int(len(ref)), "denominator": int(denominator),
              "coverage": float(len(ref) / denominator) if denominator else None}
    if not len(ref):
        result.update(spearman=None, pearson=None, rmse=None, medianAbsoluteError=None,
                      medianSignedError=None, directionalN=0, directionalCorrect=0,
                      directionAgreement=None, descriptiveGate=False)
        return result
    result.update(spearman=float(pd.Series(ref).corr(pd.Series(est), method="spearman")) if len(ref) > 1 else None,
                 pearson=float(pd.Series(ref).corr(pd.Series(est), method="pearson")) if len(ref) > 1 else None,
                 rmse=float(np.sqrt(np.mean(error ** 2))),
                 medianAbsoluteError=float(np.median(np.abs(error))),
                 medianSignedError=float(np.median(error)))
    selected = np.abs(ref) >= 1
    result["directionalN"] = int(selected.sum())
    result["directionalCorrect"] = int(np.sum(np.sign(ref[selected]) == np.sign(est[selected])))
    result["directionAgreement"] = (float(result["directionalCorrect"] / result["directionalN"])
                                     if result["directionalN"] else None)
    result["descriptiveGate"] = bool(result["coverage"] >= .90 and
                                      result["spearman"] is not None and result["spearman"] >= .90 and
                                      result["directionAgreement"] is not None and result["directionAgreement"] >= .90)
    return result

summary = {"protocolSHA256": freeze["protocolSHA256"],
           "referenceFreezeSHA256": h(r / "reference-freeze.json"),
           "qPCRSHA256": h(qpath), "qPCRRows": len(q),
           "methods": list(METHODS), "sites": {}}
exclusion_counts = collections.Counter(q_reasons)
(r / "scores").mkdir(exist_ok=True)
for site in SITES:
    count_dir = r / "counts" / site
    feature = pd.read_csv(count_dir / "features.tsv", sep="\t", dtype=str, keep_default_na=False)
    source_ids = collections.defaultdict(list)
    for i, eid in enumerate(feature["EntrezID"]):
        if not missing_literal(eid):
            source_ids[eid].append(i)
    source_unique = {eid: indices[0] for eid, indices in source_ids.items() if len(indices) == 1}
    eligible_ids = {eid for eid in eligible_q if eid in source_unique}
    native = json.loads(gzip.decompress((r / "native" / f"{site}.json.gz").read_bytes()))
    result = native["result"]
    native_features = result["features"]
    diagnostics = result["negativeBinomial"]["features"]
    estimates = {method: {} for method in METHODS}
    for i, item in enumerate(native_features):
        eid = feature.loc[i, "EntrezID"]
        if eid not in eligible_ids:
            continue
        if item.get("status") == "tested" and finite(item.get("log2FoldChange")):
            estimates["nativeMLE"][eid] = float(item["log2FoldChange"])
        fit = diagnostics[i].get("effectShrinkageFit")
        if fit and fit.get("converged") and finite(fit.get("effect")):
            estimates["nativeShrinkage"][eid] = float(fit["effect"]) / math.log(2)
    for method, (filename, column) in R_METHODS.items():
        table = pd.read_csv(r / a.reference_dir / site / filename, sep="\t", compression="gzip",
                            dtype={"featureID": str})
        assert len(table) == 25794 and table["featureID"].is_unique
        for row in table.itertuples(index=False):
            transport = getattr(row, "featureID")
            try:
                index = int(transport.rsplit(":", 1)[1])
            except (ValueError, AttributeError):
                continue
            eid = feature.loc[index, "EntrezID"]
            value = getattr(row, column)
            if eid in eligible_ids and finite(value):
                estimates[method][eid] = float(value)
    per_method = {}
    for method in METHODS:
        ids = sorted(estimates[method])
        item = metrics([eligible_q[eid] for eid in ids],
                       [estimates[method][eid] for eid in ids], len(eligible_ids))
        item["eligibleReferenceFeatures"] = len(eligible_ids)
        item["availableFeatures"] = len(ids)
        item["status"] = "available" if ids else "unavailable"
        per_method[method] = item
    common = set.intersection(*(set(estimates[method]) for method in METHODS))
    common_metrics = {}
    for method in METHODS:
        ids = sorted(common)
        item = metrics([eligible_q[eid] for eid in ids],
                       [estimates[method][eid] for eid in ids], len(common))
        item["family"] = "commonAllMethods"
        common_metrics[method] = item
    scores = pd.DataFrame({"EntrezID": sorted(eligible_ids)})
    scores["Symbol"] = scores["EntrezID"].map(dict(zip(feature["EntrezID"], feature["Symbol"])))
    scores["referenceLog2FoldChange"] = scores["EntrezID"].map(eligible_q)
    for method in METHODS:
        scores[method] = scores["EntrezID"].map(estimates[method])
    scores.to_csv(r / "scores" / f"{site}.tsv.gz", sep="\t", index=False, compression="gzip")
    status_inventory = collections.Counter(item.get("status", "missing") for item in native_features)
    summary["sites"][site] = {"sourceFeatures": len(feature), "sourceUniqueEntrez": len(source_unique),
                              "eligibleReferenceFeatures": len(eligible_ids), "originalFamily": per_method,
                              "commonFamily": common_metrics, "commonFeatureCount": len(common),
                              "nativeStatusInventory": dict(status_inventory),
                              "sourceDuplicateEntrezRows": int(sum(len(v) for v in source_ids.values() if len(v) > 1))}
summary["qPCRExclusionCounts"] = dict(exclusion_counts)
summary["qualification"] = {"coverage": ">=0.90", "spearman": ">=0.90",
                             "directionAgreementAbsReferenceLFCAtLeast1": ">=0.90",
                             "interpretation": "descriptive technical-reference gate; no biological, FDR or clinical claim"}
(r / "score-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps({"sites": len(summary["sites"]), "qPCREligibleRows": len(eligible_q),
                  "scoreSummary": str(r / "score-summary.json")}))
