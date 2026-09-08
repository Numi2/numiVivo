#!/usr/bin/env python3
"""Verify listed evidence against an independently pinned seal, without replay.

The evidence tree must be quiescent. Exclusions and unlisted files are outside
this check. A byte-identical failed result remains failed; integrity is not a
scientific verdict or authentication of whoever created the pinned manifest.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat


SCHEMA = "numivivo.org/audit/raw-evidence-manifest/v1"
RESULT_SCHEMA = "numivivo.org/audit/evidence-integrity-verification/v1"
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
SCOPE = ("Integrity of listed entries only. Exclusions and unlisted files are not "
         "verified. Existing scientific outcomes are unchanged; no trajectory is replayed.")


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key: " + key)
        result[key] = value
    return result


def validate_manifest(data):
    manifest = json.loads(data, object_pairs_hook=unique_object)
    if not isinstance(manifest, dict) or manifest.get("schema") != SCHEMA:
        raise ValueError("unsupported evidence manifest schema")
    entries = manifest.get("files")
    if not isinstance(entries, list) or not entries:
        raise ValueError("a nonempty evidence inventory is required")
    paths = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("invalid evidence entry")
        name = entry.get("path")
        if (not isinstance(name, str) or not name or "\\" in name or "\0" in name
                or name.startswith("/") or any(p in ("", ".", "..") for p in name.split("/"))):
            raise ValueError("unsafe or noncanonical evidence path")
        if name in paths:
            raise ValueError("duplicate evidence path: " + name)
        paths.add(name)
        kind = entry.get("type")
        if kind == "file":
            if (set(entry) != {"path", "type", "bytes", "sha256"}
                    or type(entry["bytes"]) is not int or entry["bytes"] < 0
                    or not isinstance(entry["sha256"], str)
                    or not SHA256.fullmatch(entry["sha256"])):
                raise ValueError("invalid regular-file identity: " + name)
        elif kind == "symlink":
            if (set(entry) != {"path", "type", "target"}
                    or not isinstance(entry["target"], str) or not entry["target"]
                    or "\0" in entry["target"]):
                raise ValueError("invalid symlink identity: " + name)
        else:
            raise ValueError("unsupported evidence entry type: " + name)
    # An inventoried file/link cannot also be a directory leading to another row.
    for name in paths:
        if any(str(parent) in paths for parent in PurePosixPath(name).parents):
            raise ValueError("overlapping evidence paths: " + name)
    return manifest


def metadata_identity(value):
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)


def verify_entry(root_fd, entry):
    """Use directory handles so a substituted parent link is never followed."""
    parts = entry["path"].split("/")
    parent_fd = os.dup(root_fd)
    try:
        for part in parts[:-1]:
            next_fd = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                              dir_fd=parent_fd)
            os.close(parent_fd)
            parent_fd = next_fd
        leaf = parts[-1]
        before = os.stat(leaf, dir_fd=parent_fd, follow_symlinks=False)
        if entry["type"] == "symlink":
            if not stat.S_ISLNK(before.st_mode):
                raise ValueError("expected a symbolic link")
            if os.readlink(leaf, dir_fd=parent_fd) != entry["target"]:
                raise ValueError("symbolic-link target differs")
        else:
            if not stat.S_ISREG(before.st_mode):
                raise ValueError("expected a regular file")
            value = hashlib.sha256()
            # NONBLOCK also avoids waiting on a file replaced by a FIFO.
            fd = os.open(leaf, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK,
                         dir_fd=parent_fd)
            with os.fdopen(fd, "rb") as stream:
                if metadata_identity(os.fstat(stream.fileno())) != metadata_identity(before):
                    raise ValueError("file changed before reading")
                if before.st_size != entry["bytes"]:
                    raise ValueError("file size differs")
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    value.update(block)
                if metadata_identity(os.fstat(stream.fileno())) != metadata_identity(before):
                    raise ValueError("file changed while reading")
            if value.hexdigest() != entry["sha256"]:
                raise ValueError("file SHA256 differs")
        after = os.stat(leaf, dir_fd=parent_fd, follow_symlinks=False)
        if metadata_identity(after) != metadata_identity(before):
            raise ValueError("entry changed during verification")
    finally:
        os.close(parent_fd)


def verify(root, manifest_path, expected_sha256):
    if not isinstance(expected_sha256, str) or not SHA256.fullmatch(expected_sha256):
        raise ValueError("an independently recorded lowercase SHA256 is required")
    data = manifest_path.read_bytes()
    observed = hashlib.sha256(data).hexdigest()
    if observed != expected_sha256:
        raise ValueError("manifest SHA256 differs from the pinned identity")
    manifest = validate_manifest(data)
    root = root.resolve(strict=True)
    root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    failures = []
    try:
        for entry in manifest["files"]:
            try:
                verify_entry(root_fd, entry)
            except (OSError, ValueError) as error:
                failures.append(dict(path=entry["path"], error=f"{type(error).__name__}: {error}"))
    finally:
        os.close(root_fd)
    if manifest_path.read_bytes() != data:
        failures.append(dict(path=str(manifest_path), error="manifest changed during verification"))
    return dict(schema=RESULT_SCHEMA, manifestSHA256=observed, root=str(root),
        entries=len(manifest["files"]),
        regularFiles=sum(e["type"] == "file" for e in manifest["files"]),
        symlinks=sum(e["type"] == "symlink" for e in manifest["files"]),
        declaredFileBytes=sum(e.get("bytes", 0) for e in manifest["files"]),
        failures=failures, passed=not failures, exclusions=manifest.get("exclusions", []), scope=SCOPE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--out", type=Path, help="New JSON path; otherwise write to stdout")
    args = parser.parse_args()
    try:
        result = verify(args.root, args.manifest, args.expected_sha256)
    except (OSError, ValueError) as error:
        result = dict(schema=RESULT_SCHEMA, passed=False, error=f"{type(error).__name__}: {error}", scope=SCOPE)
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.out:
        with args.out.open("x") as stream:
            stream.write(encoded)
    else:
        print(encoded, end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
