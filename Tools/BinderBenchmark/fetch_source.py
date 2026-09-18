#!/usr/bin/env python3
"""Fetch a revision-pinned public table; verify its Git blob identity before publication."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import urllib.request

REVISION = "9e1b81696da46835e9e9cde9a3da976e0abc92ab"
BLOB = "d1573ba03e8322c70ccb3a40e46418e86b40e2dc"
URL = f"https://huggingface.co/datasets/Anthropic/claude-protein-binder-design/resolve/{REVISION}/data/tables/design_summary.csv"
FEATURES = ["ipsae_min_boltz2", "ipsae_min_ptxv2", "ipsae_min_ef2full"]
TARGETS = ["BBF-14", "MBP", "EGFR"]

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new source directory; parent must exist")
    args = parser.parse_args()
    out = args.output.absolute()
    if out.exists() or not out.parent.is_dir():
        parser.error("output must be new and its parent must exist")
    with urllib.request.urlopen(URL, timeout=120) as response:
        raw = response.read(64 * 1024 * 1024 + 1)
    if len(raw) > 64 * 1024 * 1024:
        raise ValueError("source exceeds the 64 MiB contract")
    git_blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if git_blob != BLOB:
        raise ValueError(f"pinned source mismatch: expected {BLOB}, received {git_blob}")
    with tempfile.TemporaryDirectory(prefix=".binder-source-", dir=out.parent) as tmp:
        staging = Path(tmp) / "source"
        staging.mkdir()
        (staging / "source.csv").write_bytes(raw)
        for assay in ("adaptyv", "twist"):
            config = {"schemaVersion": 1, "assay": assay, "targets": TARGETS, "features": FEATURES}
            (staging / f"import-{assay}.json").write_text(json.dumps(config, sort_keys=True, indent=2) + "\n")
        provenance = {"url": URL, "revision": REVISION, "gitBlobSHA1": git_blob,
                      "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw),
                      "attribution": "Anthropic, claude-protein-binder-design dataset, CC BY 4.0",
                      "scope": "BBF-14, MBP, EGFR; source preserved, other targets excluded by import config"}
        (staging / "SOURCE.json").write_text(json.dumps(provenance, sort_keys=True, indent=2) + "\n")
        if out.exists():
            raise FileExistsError(out)
        os.rename(staging, out)
    print(out)

if __name__ == "__main__":
    main()
