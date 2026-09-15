#!/usr/bin/env python3
"""Build a failure-targeting ledger from paired annotation-retention results.

This is a diagnostic organizer, not an integration method. In particular,
"supported" means that the frozen annotation-recall diagnostic had sufficient
cells and a sensitive erasure control. It does not establish that cell states
are biologically comparable across donors or that an alignment is warranted.
"""
import argparse
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


def read_json(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def comparison_key(record):
    try:
        return int(record["stratum"]), int(record["label"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("comparison is missing an integer stratum or label") from error


def comparison_index(document, name):
    comparisons = document.get("comparisons")
    if not isinstance(comparisons, list) or not comparisons:
        raise ValueError(name + " does not contain a non-empty comparisons array")
    indexed = {}
    for record in comparisons:
        key = comparison_key(record)
        if key in indexed:
            raise ValueError(name + " contains duplicate comparison " + repr(key))
        indexed[key] = record
    return indexed


def require_equal(first, second, field, key):
    if first.get(field) != second.get(field):
        raise ValueError(
            "paired results disagree on " + field + " for comparison " + repr(key)
        )


def support_status(record):
    sufficient = bool(record.get("sufficientSupport"))
    sensitive = bool(record.get("controlSensitive"))
    if sufficient and sensitive:
        return "supported_control_sensitive"
    reasons = []
    if not sufficient:
        reasons.append("insufficient_support")
    if not sensitive:
        reasons.append("insufficient_control_sensitivity")
    return "unsupported_" + "_and_".join(reasons)


def outcome(native, candidate):
    if support_status(native) != "supported_control_sensitive":
        return "unsupported_no_promotion_decision"
    native_failure = bool(native.get("marginFailure"))
    candidate_failure = bool(candidate.get("marginFailure"))
    if native_failure and candidate_failure:
        return "shared_margin_failure"
    if native_failure:
        return "resolved_by_candidate"
    if candidate_failure:
        return "new_candidate_margin_failure"
    return "neither_margin_failure"


def metric_view(record):
    return {
        "baselineMeanRecall": record["baselineMeanRecall"],
        "resultMeanRecall": record["candidateMeanRecall"],
        "meanRecallLoss": record["meanRecallLoss"],
        "maximumFoldRecallLoss": record["maximumFoldRecallLoss"],
        "marginFailure": record["marginFailure"],
    }


def count_outcomes(records):
    counts = Counter(record["outcome"] for record in records)
    return {
        "sharedMarginFailures": counts["shared_margin_failure"],
        "newCandidateMarginFailures": counts["new_candidate_margin_failure"],
        "resolvedByCandidate": counts["resolved_by_candidate"],
        "neitherMarginFailure": counts["neither_margin_failure"],
        "unsupportedNoPromotionDecision": counts["unsupported_no_promotion_decision"],
    }


def freeze_extension(base_freeze, candidate_freeze, candidate_matrix, candidate_sha):
    if not isinstance(base_freeze, dict) or not isinstance(candidate_freeze, dict):
        raise ValueError("freeze documents must be JSON objects")
    base_without_matrices = {
        key: value for key, value in base_freeze.items() if key != "matrices"
    }
    candidate_without_matrices = {
        key: value for key, value in candidate_freeze.items() if key != "matrices"
    }
    if base_without_matrices != candidate_without_matrices:
        raise ValueError("candidate freeze changes fields outside matrices")
    base_matrices = base_freeze.get("matrices")
    candidate_matrices = candidate_freeze.get("matrices")
    if not isinstance(base_matrices, dict) or not isinstance(candidate_matrices, dict):
        raise ValueError("freeze documents must contain matrices objects")
    for matrix, value in base_matrices.items():
        if candidate_matrices.get(matrix) != value:
            raise ValueError("candidate freeze changes existing matrix " + matrix)
    additions = sorted(set(candidate_matrices) - set(base_matrices))
    if additions != [candidate_matrix]:
        raise ValueError(
            "candidate freeze additions must be exactly " + repr([candidate_matrix])
        )
    candidate_entry = candidate_matrices[candidate_matrix]
    if not isinstance(candidate_entry, dict) or candidate_entry.get("SHA256") != candidate_sha:
        raise ValueError("candidate freeze does not bind the candidate score matrix")
    return {
        "relationship": "candidate_freeze_extends_native_freeze",
        "addedMatrix": candidate_matrix,
        "addedMatrixSHA256": candidate_sha,
    }


def build_summary(
    ledger,
    native_document,
    candidate_document,
    input_hashes,
    freeze_provenance=None,
):
    native_index = comparison_index(native_document, "native result")
    candidate_index = comparison_index(candidate_document, "candidate result")
    for field in ("cells", "foldCount"):
        if native_document.get(field) != candidate_document.get(field):
            raise ValueError("paired results disagree on " + field)
    if native_index.keys() != candidate_index.keys():
        missing = sorted(native_index.keys() - candidate_index.keys())
        extra = sorted(candidate_index.keys() - native_index.keys())
        raise ValueError(
            "paired comparison keys differ; missing=" + repr(missing) +
            " extra=" + repr(extra)
        )

    labels = ledger.get("labels")
    strata = ledger.get("strata")
    if not isinstance(labels, list) or not isinstance(strata, list):
        raise ValueError("ledger must contain labels and strata arrays")

    records = []
    for key in sorted(native_index):
        native = native_index[key]
        candidate = candidate_index[key]
        for field in (
            "labelName",
            "rare",
            "sufficientSupport",
            "controlSensitive",
            "baselineMeanRecall",
            "erasureMeanRecall",
        ):
            require_equal(native, candidate, field, key)
        stratum, label = key
        if not 0 <= stratum < len(strata) or not 0 <= label < len(labels):
            raise ValueError("ledger lacks metadata for comparison " + repr(key))
        if labels[label] != native["labelName"]:
            raise ValueError("ledger label disagrees with comparison " + repr(key))
        metadata = strata[stratum]
        if not isinstance(metadata, dict):
            raise ValueError("ledger stratum is not an object for comparison " + repr(key))
        records.append({
            "stratumIndex": stratum,
            "preparation": metadata["preparation"],
            "treatment": metadata["treatment"],
            "labelIndex": label,
            "label": native["labelName"],
            "rare": bool(native["rare"]),
            "diagnosticSupport": {
                "status": support_status(native),
                "sufficientSupport": bool(native["sufficientSupport"]),
                "controlSensitive": bool(native["controlSensitive"]),
            },
            "outcome": outcome(native, candidate),
            "native": metric_view(native),
            "candidate": metric_view(candidate),
        })

    by_stratum = defaultdict(list)
    for record in records:
        by_stratum[record["stratumIndex"]].append(record)
    stratum_summary = []
    for stratum_index in sorted(by_stratum):
        group = by_stratum[stratum_index]
        row = {
            "stratumIndex": stratum_index,
            "preparation": group[0]["preparation"],
            "treatment": group[0]["treatment"],
            "comparisons": len(group),
            "supportedControlSensitive": sum(
                record["diagnosticSupport"]["status"] ==
                "supported_control_sensitive"
                for record in group
            ),
        }
        row.update(count_outcomes(group))
        stratum_summary.append(row)

    supported = [
        record for record in records
        if record["diagnosticSupport"]["status"] == "supported_control_sensitive"
    ]
    counts = {
        "allComparisons": len(records),
        "supportedControlSensitive": len(supported),
        "unsupported": len(records) - len(supported),
        "nativeMarginFailures": sum(
            record["native"]["marginFailure"] for record in supported
        ),
        "candidateMarginFailures": sum(
            record["candidate"]["marginFailure"] for record in supported
        ),
    }
    counts.update(count_outcomes(records))

    return {
        "schemaVersion": 1,
        "scope": (
            "Failure-targeting ledger for the frozen transductive author-annotation "
            "retention diagnostic. It is not a biological comparability, phenotype, "
            "or prospective integration qualification."
        ),
        "inputSHA256": input_hashes,
        "freezeProvenance": freeze_provenance,
        "matrices": {
            "native": {
                "id": native_document.get("matrix"),
                "scoreSHA256": native_document.get("matrixSHA256"),
            },
            "candidate": {
                "id": candidate_document.get("matrix"),
                "scoreSHA256": candidate_document.get("matrixSHA256"),
            },
        },
        "counts": counts,
        "strata": stratum_summary,
        "records": records,
        "decision": {
            "candidatePromotion": "reject",
            "nextStep": (
                "Investigate shared and newly introduced supported failures before "
                "another correction sweep. Keep unsupported comparisons explicit; "
                "their absence from a promotion decision is not evidence that the "
                "underlying states are comparable or safe to align."
            ),
        },
        "limitations": [
            (
                "Diagnostic support is a cell-count and erasure-control criterion, "
                "not evidence of biological comparability across donors."
            ),
            (
                "Author labels are recoverability controls, not biological ground truth."
            ),
            (
                "The ledger does not apply, tune, or validate a correction method."
            ),
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", required=True, type=Path)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--native-freeze", required=True, type=Path)
    parser.add_argument("--candidate-freeze", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    native_document = read_json(args.native)
    candidate_document = read_json(args.candidate)
    native_freeze_hash = sha256(args.native_freeze)
    candidate_freeze_hash = sha256(args.candidate_freeze)
    if native_document.get("freezeSHA256") != native_freeze_hash:
        raise ValueError("native result does not match the supplied native freeze")
    if candidate_document.get("freezeSHA256") != candidate_freeze_hash:
        raise ValueError("candidate result does not match the supplied candidate freeze")
    native_freeze = read_json(args.native_freeze)
    candidate_freeze = read_json(args.candidate_freeze)
    native_matrix = native_document.get("matrix")
    native_matrix_entry = native_freeze.get("matrices", {}).get(native_matrix)
    if not isinstance(native_matrix_entry, dict) or (
        native_matrix_entry.get("SHA256") != native_document.get("matrixSHA256")
    ):
        raise ValueError("native freeze does not bind the native score matrix")
    freeze_provenance = freeze_extension(
        native_freeze,
        candidate_freeze,
        candidate_document.get("matrix"),
        candidate_document.get("matrixSHA256"),
    )
    freeze_provenance["nativeFreezeSHA256"] = native_freeze_hash
    freeze_provenance["candidateFreezeSHA256"] = candidate_freeze_hash
    summary = build_summary(
        read_json(args.ledger),
        native_document,
        candidate_document,
        {
            "ledger": sha256(args.ledger),
            "native": sha256(args.native),
            "candidate": sha256(args.candidate),
            "nativeFreeze": native_freeze_hash,
            "candidateFreeze": candidate_freeze_hash,
        },
        freeze_provenance,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary["counts"], sort_keys=True))


if __name__ == "__main__":
    main()
