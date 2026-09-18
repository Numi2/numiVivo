#!/usr/bin/env python3
"""Make an explicit whole-target holdout plan without inspecting outcome values."""
import argparse
import json
from pathlib import Path

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("bundle", type=Path)
    p.add_argument("--test-target", action="append", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--top-k", type=int, default=10)
    args = p.parse_args()
    config = json.loads((args.bundle / "import.json").read_text())
    imported = json.loads((args.bundle / "imported.json").read_text())
    tests = set(args.test_target)
    available = set(config["targets"])
    if not tests < available or len(tests) != len(args.test_target) or args.top_k < 1:
        p.error("test targets must be a unique, nonempty proper subset of imported targets; K must be positive")
    plan = {"schemaVersion": 1, "sourceSHA256": imported["dataset"]["sourceSHA256"],
            "trainingTargets": sorted(available - tests), "testTargets": sorted(tests),
            "baselineFeature": "ipsae_min_boltz2", "modelFeatures": config["features"],
            "topK": args.top_k, "ridgePenalty": 0.1}
    with args.output.open("x") as handle:
        json.dump(plan, handle, sort_keys=True, indent=2)
        handle.write("\n")
    print("Plan written; freeze it before evaluating or examining test outcomes.")

if __name__ == "__main__":
    main()
