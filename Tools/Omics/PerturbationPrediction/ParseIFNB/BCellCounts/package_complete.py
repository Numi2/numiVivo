"""Archive a terminal study under its controller lock; preserve incomplete attempts."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import tarfile
from verify_complete import review, sha, require

def inventory(root):
    entries = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if "__pycache__" in relative.parts or relative.as_posix() == "controller.lock":
            continue
        require(not path.is_symlink(), "symlink in evidence: " + str(relative))
        if path.is_file():
            entries[relative.as_posix()] = dict(bytes=path.stat().st_size, SHA256=sha(path))
        else:
            require(path.is_dir(), "nonregular evidence: " + str(relative))
    require(bool(entries), "empty evidence inventory")
    return entries

def verify_contents(archive, expected):
    seen = set()
    with tarfile.open(archive, "r|gz") as stream:
        for member in stream:
            require(member.isfile() and member.name in expected and member.name not in seen,
                    "unexpected, duplicate or nonregular archive member: " + member.name)
            spec = expected[member.name]
            require(member.size == spec["bytes"], "archive member size mismatch")
            digest = hashlib.sha256()
            with stream.extractfile(member) as content:
                for block in iter(lambda: content.read(1048576), b""):
                    digest.update(block)
            require(digest.hexdigest() == spec["SHA256"], "archive member hash mismatch")
            seen.add(member.name)
    require(seen == set(expected), "missing archive members")

def package(root, manifest, output):
    root = root.resolve()
    output = output.resolve()
    require(not output.is_relative_to(root), "archive must be outside the frozen study")
    require(output.parent.is_dir(), "output directory does not exist")
    sidecar = output.with_name(output.name + ".manifest.json")
    partial = output.with_name(output.name + ".partial")
    require(not any(p.exists() for p in (output, sidecar, partial)), "output already exists")
    # Opening read-only does not truncate or replace the live driver's lock file.
    with (root / "controller.lock").open("rb") as lock:
        fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
        result = review(root, manifest)
        entries = inventory(root)
        # No uncompressed intermediate and no extraction into the study.
        with partial.open("xb") as destination:
            with tarfile.open(fileobj=destination, mode="w|gz") as archive:
                for name in entries:
                    archive.add(root / name, arcname=name, recursive=False)
        require(inventory(root) == entries, "evidence changed during packaging")
        verify_contents(partial, entries)
        receipt = dict(status="verified-complete-count-archive", review=result,
                       preparationManifestSHA256=sha(manifest),
                       archiveSHA256=sha(partial), archiveBytes=partial.stat().st_size,
                       members=entries, excluded=["controller.lock", "**/__pycache__/**"],
                       allStudyDependenciesEmbedded=True, predictionScored=False)
        # Rename only after content verification. A failed attempt retains .partial.
        partial.rename(output)
        with sidecar.open("x") as stream:
            json.dump(receipt, stream, indent=2)
            stream.write("\n")
        return receipt

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("--preparation-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    receipt = package(args.study, args.preparation_manifest, args.output)
    print(json.dumps({k: v for k, v in receipt.items() if k != "members"}, indent=2))

if __name__ == "__main__":
    main()

