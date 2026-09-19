#!/usr/bin/env python3
"""Fetch an outcome-blind, bounded public structural development panel.

Six hash-selected designs per existing BBF-14/MBP/EGFR target, one published
seed-best Boltz-2 model per design. This is source acquisition, not inference,
new-target validation, or improved candidate selection. No assay label or score
is used to choose the panel. Original manifest, construct FASTA and bytes remain.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import tempfile
import urllib.parse
import urllib.request

from fetch_source import BLOB, REVISION, TARGETS, URL

BASE = f"https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/{REVISION}/"
PER_TARGET = 6
SELECTION = "numivivo-binder-source-panel-v1"
MAX_FILE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024
METADATA = {
    "structures_manifest.csv": "data/manifests/structures_manifest.csv",
    "target_constructs.fasta": "data/docs/insilico_target_constructs.fasta",
    "INSILICO.md": "data/docs/INSILICO.md",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        raw = response.read(MAX_FILE_BYTES + 1)
    if not raw or len(raw) > MAX_FILE_BYTES:
        raise ValueError("empty source or per-file acquisition limit exceeded")
    return raw


def rows(raw: bytes) -> list[dict[str, str]]:
    text = raw.decode("utf-8", errors="strict")
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError("missing/duplicate source columns")
    result = list(reader)
    if len(result) > 100_000 or any(None in row or any(v is None for v in row.values()) for row in result):
        raise ValueError("ragged or oversized source table")
    return result


def select(source: list[dict[str, str]]) -> list[dict[str, str]]:
    # Project out all score and outcome fields before ordering.
    identity = [{key: row[key] for key in ("uuid", "target", "sequence")}
                for row in source if row.get("target") in TARGETS]
    if len({r["uuid"] for r in identity}) != len(identity) or any(not r["uuid"] or not r["sequence"] for r in identity):
        raise ValueError("missing/duplicate candidate identity or sequence")
    selected = []
    for target in TARGETS:
        pool = [r for r in identity if r["target"] == target]
        if len(pool) < PER_TARGET:
            raise ValueError(f"insufficient source candidates for fixed panel: {target}")
        pool.sort(key=lambda r: (sha256(f"{SELECTION}\0{target}\0{r['uuid']}".encode()), r["uuid"]))
        selected.extend(pool[:PER_TARGET])
    return selected


def source_path(text: str) -> str:
    if not isinstance(text, str) or "\\" in text:
        raise ValueError("invalid published relative path")
    pieces = text.split("/")
    if any(p in ("", ".", "..") for p in pieces) or PurePosixPath(text).is_absolute():
        raise ValueError("published path is absolute or traverses directories")
    if pieces[0] == "designs":
        pieces.insert(0, "data")
    if pieces[:2] != ["data", "designs"] or pieces[-1] != "predicted_boltz2_1to1.cif":
        raise ValueError("model is not the declared published Boltz-2 1to1 source")
    return "/".join(pieces)


def declaration(source: list[dict[str, str]], manifest: list[dict[str, str]]) -> dict:
    selected = select(source)
    requested = {r["uuid"]: r for r in selected}
    models: dict[str, dict[str, str]] = {}
    for row in manifest:
        if row.get("uuid") not in requested or row.get("kind") != "predicted" or row.get("predictor") != "boltz2" or row.get("stoich") != "1to1":
            continue
        key = row["uuid"]
        if key in models or row.get("target") != requested[key]["target"]:
            raise ValueError("duplicate model or target identity mismatch")
        path = source_path(row["rel_path"])
        digest = row["sha256"]
        if len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("invalid source model digest")
        count = int(row["bytes"])
        if not 0 < count <= MAX_FILE_BYTES:
            raise ValueError("source model exceeds byte limit")
        models[key] = {**row, "resolvedPath": path}
    total = sum(int(r["bytes"]) for r in models.values())
    if total > MAX_TOTAL_BYTES:
        raise ValueError("fixed panel exceeds total source budget; never silently downsample")
    return {"schemaVersion": 1, "selection": SELECTION, "perTarget": PER_TARGET,
            "targets": TARGETS, "predictor": "boltz2", "stoichiometry": "1to1",
            "requested": selected,
            "models": [{"candidateID": r["uuid"], "target": r["target"],
                        "sourcePath": f"structures/{models[r['uuid']]['sha256']}.cif",
                        "published": models[r["uuid"]]}
                       for r in selected if r["uuid"] in models],
            "unavailable": [{"candidateID": r["uuid"], "reason": "no declared published model"}
                            for r in selected if r["uuid"] not in models],
            "limitations": ["Outcome-blind fixed subset of already-inspected development targets, not the full 270-design panel.",
                "Published seed-best Boltz-2 predictions are inputs; no model weights are executed.",
                "Acquisition does not imply successful native parsing, geometry qualification or binding prediction."]}


def acquire(output: Path) -> dict:
    output = output.absolute()
    if output.exists() or not output.parent.is_dir():
        raise ValueError("output must be new and its parent must exist")
    raw = get(URL)
    actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if actual != BLOB:
        raise ValueError("published table does not match pinned Git blob")
    metadata = {name: get(BASE + path) for name, path in METADATA.items()}
    panel = declaration(rows(raw), rows(metadata["structures_manifest.csv"]))
    if not panel["models"]:
        raise ValueError("no published models in requested panel")
    with tempfile.TemporaryDirectory(prefix=".binder-structure-panel-", dir=output.parent) as tmp:
        stage = Path(tmp) / "panel"
        stage.mkdir(); (stage / "structures").mkdir()
        (stage / "source.csv").write_bytes(raw)
        for name, data in metadata.items():
            (stage / name).write_bytes(data)
        # Fixed source/model selection is written before model download or analysis.
        (stage / "panel.json").write_text(json.dumps(panel, sort_keys=True, indent=2) + "\n")
        total = 0
        for model in panel["models"]:
            record = model["published"]
            path = urllib.parse.quote(record["resolvedPath"], safe="/")
            data = get(BASE + path)
            total += len(data)
            if len(data) != int(record["bytes"]) or sha256(data) != record["sha256"] or total > MAX_TOTAL_BYTES:
                raise ValueError(f"published model byte/hash mismatch: {model['candidateID']}")
            file = stage / model["sourcePath"]
            if file.exists() and file.read_bytes() != data:
                raise ValueError("inconsistent content-addressed model")
            file.write_bytes(data)
        provenance = {"schemaVersion": 1, "revision": REVISION, "tableURL": URL,
            "tableGitBlobSHA1": BLOB, "metadataURLs": {n: BASE + p for n, p in METADATA.items()},
            "attribution": "Anthropic, claude-protein-binder-design dataset, CC BY 4.0",
            "files": {p.relative_to(stage).as_posix(): sha256(p.read_bytes()) for p in sorted(stage.rglob("*")) if p.is_file()},
            "bytes": total, "requestedCount": len(panel["requested"]), "availableCount": len(panel["models"]),
            "evidence": "source acquisition only; no numerical or biological qualification"}
        (stage / "SOURCE.json").write_text(json.dumps(provenance, sort_keys=True, indent=2) + "\n")
        if output.exists():
            raise FileExistsError(output)
        os.rename(stage, output)
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        print(json.dumps(acquire(args.output), sort_keys=True))
    except (OSError, ValueError, KeyError) as error:
        parser.exit(2, f"published structure panel rejected: {error}\n")


if __name__ == "__main__":
    main()
