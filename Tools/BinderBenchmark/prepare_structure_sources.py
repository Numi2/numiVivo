#!/usr/bin/env python3
"""Retain original PDB/mmCIF text for native binder structure evaluation.

This is transport only: Swift's existing PDB/mmCIF owners parse coordinates and
check units/topology. Expected target sequences and source hashes must come from
an explicit input manifest; no candidate or chain is inferred from a filename.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile

MAX_SOURCE_BYTES = 16 * 1024 * 1024
MAX_TOTAL_BYTES = 64 * 1024 * 1024
MAX_INPUT_BYTES = 128 * 1024 * 1024


def read_regular(path: Path, limit: int) -> bytes:
    """Bound each read and reject symlinks and non-regular files before reading."""
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    fd = os.open(path, flags)
    with os.fdopen(fd, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or not 0 <= before.st_size <= limit:
            raise ValueError(f"unsupported file type or size: {path.name}")
        data = stream.read(limit + 1)
        after = os.fstat(stream.fileno())
        if len(data) > limit or len(data) != before.st_size or (
            before.st_size, before.st_mtime_ns, before.st_ino
        ) != (after.st_size, after.st_mtime_ns, after.st_ino):
            raise ValueError(f"source changed or exceeded limit: {path.name}")
        return data


def local_source(root: Path, name: str) -> Path:
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError("sourcePath must be a relative POSIX path")
    parts = name.split("/")
    if any(part in ("", ".", "..") for part in parts):
        raise ValueError("sourcePath cannot be absolute or traverse directories")
    path = root
    for part in parts:
        path = path / part
        if path.is_symlink():
            raise ValueError("sourcePath cannot contain symbolic links")
    return path


def prepare(manifest: dict, root: Path) -> dict:
    if not isinstance(manifest, dict) or manifest.get("schemaVersion") != 1:
        raise ValueError("source manifest schemaVersion must be 1")
    rows = manifest.get("sources")
    if not isinstance(rows, list) or not 1 <= len(rows) <= 10_000:
        raise ValueError("source count must be between 1 and 10000")
    output, identities, total = [], set(), 0
    fields = {"candidateID", "target", "sourceLabel", "format", "sha256",
              "targetChainSequences", "interfacePlan", "sourcePath"}
    for row in rows:
        if not isinstance(row, dict) or set(row) != fields:
            raise ValueError("source manifest has missing or unknown fields")
        for field in ("candidateID", "target", "sourceLabel"):
            if not isinstance(row[field], str) or not row[field] or len(row[field].encode()) > 1024:
                raise ValueError(f"invalid {field}")
        if row["candidateID"] in identities:
            raise ValueError("duplicate candidateID")
        identities.add(row["candidateID"])
        if row["format"] not in ("pdb", "mmcif"):
            raise ValueError("format must be pdb or mmcif")
        digest = row["sha256"]
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError("source SHA-256 must be declared as lowercase hex")
        plan, sequences = row["interfacePlan"], row["targetChainSequences"]
        chains = plan.get("targetChains") if isinstance(plan, dict) else None
        if (not isinstance(sequences, dict) or not sequences or not isinstance(chains, list)
                or not all(isinstance(c, str) and c for c in chains)
                or len(set(chains)) != len(chains) or set(chains) != set(sequences)
                or not all(isinstance(s, str) and 0 < len(s) <= 100_000 for s in sequences.values())):
            raise ValueError("exact expected sequence required for every target chain")
        data = read_regular(local_source(root, row["sourcePath"]), MAX_SOURCE_BYTES)
        total += len(data)
        if not data or total > MAX_TOTAL_BYTES or b"\x00" in data:
            raise ValueError("source text is empty, contains NUL or exceeds archive capacity")
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError(f"source SHA-256 mismatch: {row['candidateID']}")
        text = data.decode("utf-8", errors="strict")  # Preserve CRLF and all original bytes.
        item = {k: v for k, v in row.items() if k != "sourcePath"}
        item["contents"] = text
        output.append(item)
    return {"schemaVersion": 1, "sources": output}


def write_new(path: Path, value: dict) -> None:
    data = (json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")) + "\n").encode("utf-8")
    if len(data) > MAX_INPUT_BYTES:
        raise ValueError("encoded structural source input exceeds 128 MiB")
    # Publish a complete file without replacing an existing destination. Hard
    # linking a same-directory temporary is atomic and fails if the name exists.
    fd, temporary = tempfile.mkstemp(prefix=".binder-sources-", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        os.unlink(temporary)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        manifest = json.loads(read_regular(args.manifest, MAX_SOURCE_BYTES))
        write_new(args.output, prepare(manifest, args.source_root.resolve(strict=True)))
    except (OSError, ValueError, UnicodeError) as exc:
        parser.exit(2, f"structure source preparation rejected: {exc}\n")


if __name__ == "__main__":
    main()
