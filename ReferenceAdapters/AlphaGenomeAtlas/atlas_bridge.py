#!/usr/bin/env python3
"""Controlled AlphaGenome Atlas capture bridge for NumiVivo public-reference research.

The native NumiVivo executable prepares atlas-request.json and never performs a
network call. This bridge checks the exact request, local GRCh38 FASTA identity,
usage declaration, pinned official AlphaGenome client, and writes a bounded
atlas-capture.json plus a completion record for native import.

Predictions are molecular-impact evidence only. They do not establish antigen
presentation, immune recognition, safety, clinical benefit, or treatment choice.
"""
from __future__ import annotations

import argparse
import bisect
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from typing import Any, Callable

SDK_COMMIT = "aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d"
ATLAS_GIT_BLOB = "d413f29ebfa8f28a5d770e1dbc87cdc0e15315bb"
TERMS_URI = "https://deepmind.google.com/science/alphagenome/terms"
MAX_DOCUMENT = 32 * 1024 * 1024
MAX_REQUEST = 1024 * 1024
MAX_CELLS = 250_000
MAX_VARIANTS = 128


def encoded(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      allow_nan=False).encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Duplicate JSON key")
        result[key] = value
    return result


def loads(data: bytes) -> Any:
    if len(data) > MAX_DOCUMENT:
        raise ValueError("JSON exceeds 32 MiB")
    return json.loads(data, object_pairs_hook=no_duplicates,
                      parse_constant=lambda _: (_ for _ in ()).throw(ValueError("Nonfinite JSON number")))


def read(path: Path, limit: int = MAX_DOCUMENT) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise ValueError("A regular, nonsymlink input file is required")
    with path.open("rb") as handle:
        data = handle.read(limit + 1)
    if len(data) > limit:
        raise ValueError("Input exceeds its size bound")
    return data


def write_new(path: Path, data: bytes) -> None:
    if len(data) > MAX_DOCUMENT:
        raise ValueError("Output exceeds 32 MiB")
    with path.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def is_digest(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[a-f0-9]{64}", value) is not None


def validate_request(request: dict[str, Any]) -> None:
    if set(request) != {"schema", "binding", "requestedScorers", "ontologyTerms", "variants", "excludedCandidates"}:
        raise ValueError("Unexpected Atlas request fields")
    binding = request["binding"]
    required_binding = {"parentReportSHA256", "caseFingerprint", "dataClass", "assembly", "referenceSHA256", "hlaAlleles", "tumorRNASampleID"}
    if (request["schema"] != "numivivo.org/atlas-request/v1" or set(binding) != required_binding
            or binding["dataClass"] not in {"synthetic", "publicReference"}
            or binding["assembly"] not in {"GRCh37", "GRCh38"}
            or not all(is_digest(binding[k]) for k in ["parentReportSHA256", "caseFingerprint", "referenceSHA256"])):
        raise ValueError("Unsupported Atlas request binding")
    scorers = request["requestedScorers"]
    terms = request["ontologyTerms"]
    if (not isinstance(scorers, list) or not 1 <= len(scorers) <= 16 or len(set(scorers)) != len(scorers)
            or any(not isinstance(x, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}", x) is None for x in scorers)
            or not isinstance(terms, list) or len(terms) > 32 or len(set(terms)) != len(terms)
            or any(re.fullmatch(r"(?:UBERON|CL|CLO|EFO):[0-9]+", x) is None for x in terms)):
        raise ValueError("Invalid Atlas scorer or ontology selection")
    variants = request["variants"]
    if not isinstance(variants, list) or len(variants) > MAX_VARIANTS:
        raise ValueError("Atlas request exceeds the 128-variant bound")
    seen_keys: set[str] = set()
    seen_candidates: set[str] = set()
    for item in variants:
        required = {"key", "chromosome", "position1", "reference", "alternate", "candidateIDs"}
        if set(item) != required or not isinstance(item["candidateIDs"], list) or not item["candidateIDs"]:
            raise ValueError("Malformed Atlas query")
        key = f"GRCh38|{item['chromosome']}|{item['position1']}|{item['reference']}|{item['alternate']}"
        if (binding["assembly"] != "GRCh38" or item["key"] != key or item["key"] in seen_keys
                or item["chromosome"] not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
                or type(item["position1"]) is not int or not 1 <= item["position1"] <= 300_000_000
                or item["reference"] not in "ACGT" or item["alternate"] not in "ACGT"
                or item["reference"] == item["alternate"]):
            raise ValueError("Unsupported or inconsistent Atlas variant identity")
        seen_keys.add(item["key"])
        for candidate in item["candidateIDs"]:
            if not is_digest(candidate) or candidate in seen_candidates:
                raise ValueError("Duplicate or invalid candidate mapping")
            seen_candidates.add(candidate)
    excluded = request["excludedCandidates"]
    if not isinstance(excluded, dict) or len(excluded) > 10_000:
        raise ValueError("Invalid excluded-candidate map")
    for candidate, reason in excluded.items():
        if not is_digest(candidate) or candidate in seen_candidates or not isinstance(reason, str) or not reason or len(reason.encode()) > 512:
            raise ValueError("Invalid excluded-candidate record")
        seen_candidates.add(candidate)


def validate_usage(usage: dict[str, Any], terms: bytes) -> None:
    required = {"purpose", "termsURI", "termsSHA256", "acknowledgedBy", "externalSharingApproved", "snapshotStorageApproved", "modelTrainingAllowed"}
    if (set(usage) != required or usage["purpose"] != "noncommercialResearch"
            or usage["termsURI"] != TERMS_URI or usage["termsSHA256"] != sha(terms)
            or not isinstance(usage["acknowledgedBy"], str) or not 1 <= len(usage["acknowledgedBy"]) <= 128
            or usage["externalSharingApproved"] is not True or usage["snapshotStorageApproved"] is not True
            or usage["modelTrainingAllowed"] is not False):
        raise ValueError("Explicit noncommercial research usage and matching reviewed terms bytes are required")


def verify_reference(fasta: Path, request: dict[str, Any]) -> list[bool]:
    if fasta.is_symlink() or not fasta.is_file():
        raise ValueError("Reference must be a regular uncompressed FASTA")
    positions: dict[str, list[int]] = {}
    for query in request["variants"]:
        positions.setdefault(query["chromosome"], []).append(query["position1"])
    positions = {key: sorted(set(value)) for key, value in positions.items()}
    bases: dict[tuple[str, int], str] = {}
    visited: set[str] = set()
    digest = hashlib.sha256()
    contig = ""
    offset = 0
    with fasta.open("rb") as handle:
        while True:
            raw = handle.readline(16 * 1024 * 1024 + 1)
            if not raw:
                break
            if len(raw) > 16 * 1024 * 1024:
                raise ValueError("FASTA line exceeds 16 MiB")
            digest.update(raw)
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                name = line[1:].split(maxsplit=1)[0].decode("ascii")
                contig = name if name.startswith("chr") else "chr" + name
                offset = 0
                if contig in positions:
                    if contig in visited:
                        raise ValueError("Ambiguous reference contig aliases")
                    visited.add(contig)
                continue
            if not contig or re.fullmatch(rb"[ACGTNacgtn]*", line) is None:
                raise ValueError("Invalid uncompressed reference FASTA")
            if contig in positions:
                targets = positions[contig]
                lo = bisect.bisect_right(targets, offset)
                hi = bisect.bisect_right(targets, offset + len(line))
                for position in targets[lo:hi]:
                    bases[(contig, position)] = chr(line[position - offset - 1]).upper()
            offset += len(line)
    if digest.hexdigest() != request["binding"]["referenceSHA256"]:
        raise ValueError("FASTA digest differs from the case reference; refusing all queries")
    return [bases.get((query["chromosome"], query["position1"])) == query["reference"] for query in request["variants"]]


def normalize_json(value: Any, depth: int = 0) -> Any:
    if depth > 8:
        raise ValueError("Provider metadata nesting exceeds bounds")
    if value is None:
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if not math.isfinite(value):
            raise ValueError("Nonfinite provider metadata")
        return value
    if isinstance(value, str):
        if "\0" in value or len(value.encode()) > 4096:
            raise ValueError("Oversized provider metadata")
        return value
    if isinstance(value, (list, tuple)):
        if len(value) > 128:
            raise ValueError("Provider metadata array exceeds bounds")
        return [normalize_json(item, depth + 1) for item in value]
    if isinstance(value, dict):
        if len(value) > 128:
            raise ValueError("Provider metadata object exceeds bounds")
        result = {}
        for key, item in value.items():
            key = str(key)
            if len(key.encode()) > 256:
                raise ValueError("Provider metadata key exceeds bounds")
            result[key] = normalize_json(item, depth + 1)
        return result
    return normalize_json(str(value), depth + 1)


def finite_score(value: Any) -> float | None:
    number = float(value)
    if math.isnan(number):
        return None
    if not math.isfinite(number):
        raise ValueError("Infinite provider score")
    return number


def sdk_identity() -> tuple[str, dict[str, Any]]:
    from alphagenome.atlas import atlas
    distribution = importlib.metadata.distribution("alphagenome")
    direct = json.loads(distribution.read_text("direct_url.json") or "{}")
    if (direct.get("vcs_info", {}).get("commit_id") != SDK_COMMIT
            or direct.get("url", "").removesuffix(".git") != "https://github.com/google-deepmind/alphagenome"):
        raise ValueError("Install the exact official AlphaGenome SDK commit from requirements.txt")
    module = Path(atlas.__file__).read_bytes()
    git_blob = hashlib.sha1(b"blob " + str(len(module)).encode() + b"\0" + module).hexdigest()
    if git_blob != ATLAS_GIT_BLOB:
        raise ValueError("Installed Atlas client source differs from the inspected client source")
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}", distribution.version) is None:
        raise ValueError("Installed SDK version is not representable by the native provenance token")
    return distribution.version, atlas.scorer_metadata if hasattr(atlas, "scorer_metadata") else {}


def metadata_catalog(client: Any) -> dict[str, Any]:
    metadata = client.scorer_metadata()
    result: dict[str, Any] = {}
    for name, descriptor in sorted(metadata.items()):
        tracks = []
        for _, row in descriptor.track_metadata.iterrows():
            tracks.append({str(k): normalize_json(v) for k, v in dict(row).items()})
        result[name] = {"isSigned": bool(descriptor.is_signed), "tracks": tracks}
    return result


def matrix_to_table(name: str, matrix: Any, query: dict[str, Any], descriptor: Any) -> dict[str, Any] | None:
    if matrix.shape != (len(matrix.obs), len(matrix.var)):
        raise ValueError("Provider score matrix dimensions disagree")
    if len(matrix.obs) == 0 or len(matrix.var) == 0:
        return None
    observations: list[dict[str, Any]] = []
    for i in range(len(matrix.obs)):
        values = dict(matrix.obs.iloc[i])
        returned = values.pop("variant", None)
        if returned is None or (returned.chromosome, returned.position, returned.reference_bases, returned.alternate_bases) != (
                query["chromosome"], query["position1"], query["reference"], query["alternate"]):
            raise ValueError("Provider observation variant differs from the requested variant")
        observations.append({"variantKey": query["key"], "metadata": {str(k): normalize_json(v) for k, v in values.items()}})
    tracks = [{str(k): normalize_json(v) for k, v in dict(matrix.var.iloc[j]).items()} for j in range(len(matrix.var))]
    raw: list[list[float | None]] = []
    quantiles: list[list[float | None]] | None = [] if "quantiles" in matrix.layers else None
    for i in range(len(matrix.obs)):
        raw_row: list[float | None] = []
        quantile_row: list[float | None] = []
        for j in range(len(matrix.var)):
            r = finite_score(matrix.X[i, j])
            if r is not None and not descriptor.is_signed and r < 0:
                raise ValueError("Negative raw score from an unsigned scorer")
            raw_row.append(r)
            if quantiles is not None:
                q = finite_score(matrix.layers["quantiles"][i, j])
                if q is not None and not (-1 if descriptor.is_signed else 0) <= q <= 1:
                    raise ValueError("Calibrated quantile is outside its signed/unsigned range")
                if q is not None and r is None:
                    raise ValueError("Calibrated score cannot conceal an absent raw score")
                quantile_row.append(q)
        raw.append(raw_row)
        if quantiles is not None:
            quantiles.append(quantile_row)
    return {"scorer": name, "isSigned": bool(descriptor.is_signed), "observations": observations,
            "tracks": tracks, "rawScores": raw, "quantiles": quantiles}


def worker(payload: dict[str, Any]) -> dict[str, Any]:
    from alphagenome.atlas import atlas
    from alphagenome.data import genome
    version, _ = sdk_identity()
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        raise PermissionError("Missing API key")
    client = atlas.create(key, timeout=20)
    metadata = client.scorer_metadata()
    if payload["operation"] == "catalog":
        return {"clientVersion": version, "clientCommit": SDK_COMMIT, "scorers": metadata_catalog(client)}
    query, requested, terms = payload["query"], payload["requestedScorers"], payload["ontologyTerms"]
    if any(name not in metadata for name in requested):
        raise ValueError("Requested scorer is not present in the live catalog")
    variant = genome.Variant(chromosome=query["chromosome"], position=query["position1"],
                             reference_bases=query["reference"], alternate_bases=query["alternate"])
    scores = client.query_variant(variant, requested_scorers=requested, ontology_terms=terms or None)
    tables = []
    for name, matrix in sorted(scores.items()):
        if name not in requested or name not in metadata:
            raise ValueError("Provider returned an unrequested scorer")
        table = matrix_to_table(name, matrix, query, metadata[name])
        if table is not None:
            tables.append(table)
    return {"clientVersion": version, "clientCommit": SDK_COMMIT, "tables": tables}


def invoke(payload: dict[str, Any], timeout: int) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="numivivo-atlas-") as directory:
        request_path = Path(directory) / "request.json"
        write_new(request_path, encoded(payload))
        result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "_worker", str(request_path)],
                                capture_output=True, check=False, timeout=timeout)
        if result.returncode != 0:
            raise RuntimeError("Atlas worker failed; check SDK installation, permission, quota and selected scorers")
        return loads(result.stdout)


def fetch(request_data: bytes, fasta: Path, usage: dict[str, Any], terms: bytes, timeout: int,
          call: Callable[[dict[str, Any], int], dict[str, Any]] = invoke) -> dict[str, Any]:
    if len(request_data) > MAX_REQUEST:
        raise ValueError("Atlas request exceeds 1 MiB")
    request = loads(request_data)
    validate_request(request)
    if request["binding"]["dataClass"] != "publicReference":
        raise ValueError("Live Atlas retrieval is restricted to public-reference research cases")
    validate_usage(usage, terms)
    reference_ok = verify_reference(fasta, request)
    catalog = call({"operation": "catalog"}, timeout)
    if catalog.get("clientCommit") != SDK_COMMIT or any(name not in catalog.get("scorers", {}) for name in request["requestedScorers"]):
        raise ValueError("Pinned Atlas client/catalog does not contain every requested scorer")
    outcomes: list[dict[str, Any]] = []
    tables: list[dict[str, Any]] = []
    cells = 0
    cancelled = False
    for query, reference_matches in zip(request["variants"], reference_ok, strict=True):
        if cancelled:
            outcomes.append({"variantKey": query["key"], "status": "failed", "diagnostic": "cancelled"})
            continue
        if not reference_matches:
            outcomes.append({"variantKey": query["key"], "status": "failed", "diagnostic": "referenceMismatch"})
            continue
        try:
            answer = call({"operation": "query", "query": query, "requestedScorers": request["requestedScorers"],
                           "ontologyTerms": request["ontologyTerms"]}, timeout)
            if answer.get("clientCommit") != SDK_COMMIT or answer.get("clientVersion") != catalog.get("clientVersion"):
                raise ValueError("Provider/client identity changed during retrieval")
            new_tables = answer.get("tables", [])
            if not isinstance(new_tables, list):
                raise ValueError("Malformed provider response")
            added = 0
            for table in new_tables:
                if table["scorer"] not in request["requestedScorers"]:
                    raise ValueError("Unrequested scorer table")
                product = len(table["observations"]) * len(table["tracks"])
                added += product
            if cells + added > MAX_CELLS:
                raise ValueError("Atlas capture exceeds the native 250,000-cell budget")
            cells += added
            tables.extend(new_tables)
            outcomes.append({"variantKey": query["key"], "status": "available" if added else "noData", "diagnostic": None})
        except KeyboardInterrupt:
            cancelled = True
            outcomes.append({"variantKey": query["key"], "status": "failed", "diagnostic": "cancelled"})
        except subprocess.TimeoutExpired:
            outcomes.append({"variantKey": query["key"], "status": "failed", "diagnostic": "deadlineExceeded"})
        except Exception:
            outcomes.append({"variantKey": query["key"], "status": "failed", "diagnostic": "providerOrContractFailure"})
    return {"schema": "numivivo.org/atlas-capture/v1", "requestSHA256": sha(request_data), "mode": "live",
            "provider": "google-deepmind/alphagenome-atlas", "clientCommit": SDK_COMMIT,
            "clientVersion": catalog["clientVersion"], "adapterSHA256": sha(Path(__file__).read_bytes()),
            "retrievedAt": now(), "serviceVersion": None, "serviceVersionStatus": "not-exposed-by-pinned-api",
            "termsURI": TERMS_URI, "permittedUse": "noncommercialResearch", "trainingPermitted": False,
            "referenceCheck": "local-fasta-sha256-and-ref-bases", "referenceSHA256": request["binding"]["referenceSHA256"],
            "outcomes": outcomes, "tables": tables}


def completion(files: dict[str, bytes]) -> bytes:
    return encoded({"schema": "numivivo.org/external-capture-complete/v1",
                    "files": {name: sha(data) for name, data in sorted(files.items())}})


def main() -> int:
    if len(sys.argv) == 3 and sys.argv[1] == "_worker":
        try:
            sys.stdout.buffer.write(encoded(worker(loads(read(Path(sys.argv[2]))))))
            return 0
        except Exception:
            print("Atlas SDK call failed; details withheld to protect credentials.", file=sys.stderr)
            return 1
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=["catalog", "fetch"])
    parser.add_argument("--request", type=Path)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--usage", type=Path, required=True)
    parser.add_argument("--terms", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=120)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink() or not 1 <= args.timeout_seconds <= 600:
        parser.error("Output must be new and timeout must be 1–600 seconds")
    usage_data, terms = read(args.usage, 128 * 1024), read(args.terms, 1024 * 1024)
    usage = loads(usage_data)
    validate_usage(usage, terms)
    if not os.environ.get("ALPHAGENOME_API_KEY"):
        parser.error("ALPHAGENOME_API_KEY is required; never store it in a manifest or command argument")
    if args.operation == "catalog":
        result = invoke({"operation": "catalog"}, args.timeout_seconds)
        output = encoded(result)
        args.output.mkdir(mode=0o700)
        write_new(args.output / "catalog.json", output)
        write_new(args.output / "usage.json", usage_data)
        write_new(args.output / "terms.txt", terms)
        print(json.dumps({"status": "catalog-written", "file": str(args.output / "catalog.json"), "sha256": sha(output)}))
        return 0
    if args.request is None or args.fasta is None:
        parser.error("fetch requires --request and --fasta")
    request_data = read(args.request, MAX_REQUEST)
    capture = fetch(request_data, args.fasta, usage, terms, args.timeout_seconds)
    capture_data = encoded(capture)
    evidence = {"atlas-request.json": request_data, "atlas-capture.json": capture_data}
    args.output.mkdir(mode=0o700)
    write_new(args.output / "atlas-request.json", request_data)
    write_new(args.output / "atlas-capture.json", capture_data)
    write_new(args.output / "usage.json", usage_data)
    write_new(args.output / "terms.txt", terms)
    write_new(args.output / "complete.json", completion(evidence))  # written last
    print(json.dumps({"status": "atlas-capture-written", "captureSHA256": sha(capture_data),
                      "clinicalQualification": "none", "serviceVersionPinned": False}))
    return 2 if any(item["status"] == "failed" for item in capture["outcomes"]) else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError, KeyError, TypeError) as error:
        print(f"atlas-bridge: {error}", file=sys.stderr)
        raise SystemExit(1)
