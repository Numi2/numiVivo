"""Replay the retained real window through the actual CLI; no network required."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import struct
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("cli", type=Path)
parser.add_argument("evidence", type=Path, help="unpacked original evidence.tar.gz")
parser.add_argument("output", type=Path, help="new directory")
args = parser.parse_args()
args.output.mkdir(exist_ok=False)
out = args.output.resolve()
evidence = args.evidence.resolve()
plan = json.loads(gzip.decompress(Path(__file__).with_name("stream-plan.json.gz").read_bytes()))
rows = json.loads((evidence / "real-input.json").read_text())["rows"]
reference = json.loads((evidence / "real-output.json").read_text())
payload = b"".join(struct.pack("<QQ", i | (j << 32), c)
                   for i, row in enumerate(rows)
                   for j, c in zip(row["indices"], row["counts"]))
digest = hashlib.sha256(payload).hexdigest()
(out / "plan.json").write_text(json.dumps(plan))
base = [str(args.cli.resolve()), "singlecell-celltypist-stream", "--model",
        str(evidence / "model.json"), "--plan", str(out / "plan.json"),
        "--stream-sha256", digest, "--output", str(out / "bundle")]
result = subprocess.run(base, input=payload, capture_output=True, check=True)
receipt = json.loads(result.stdout)
actual = [json.loads(line) for line in (out / "bundle/results.jsonl").read_text().splitlines()]
assert len(actual) == len(reference) == len(rows) == 267
for i, (a, b) in enumerate(zip(actual, reference)):
    assert a.pop("sourceRow") == i and a == b
for name, field in [("model.json", "modelSHA256"), ("plan.json", "planSHA256"),
                    ("results.jsonl", "resultsSHA256")]:
    assert hashlib.sha256((out / "bundle" / name).read_bytes()).hexdigest() == receipt[field]
assert receipt["streamSHA256"] == digest and receipt["records"] == len(payload) // 16
failures = []
for name in ["wrong-digest", "truncated", "wrong-total", "wrong-cardinality", "existing"]:
    command = base.copy()
    command[-1] = str(out / name) if name != "existing" else str(out / "bundle")
    data = payload
    expected = {"wrong-digest": "differs from frozen input", "truncated": "invalidCounts",
                "wrong-total": "stream row total", "wrong-cardinality": "stream row cardinality",
                "existing": "output path"}[name]
    if name == "wrong-digest":
        command[7] = "0" * 64
    elif name == "truncated":
        data = payload[:-1]
    elif name in ("wrong-total", "wrong-cardinality"):
        changed = json.loads(json.dumps(plan))
        field = "rowTotals" if name == "wrong-total" else "rowNonzeros"
        changed[field][0] += 1 if name == "wrong-total" else -1
        path = out / (name + ".json")
        path.write_text(json.dumps(changed))
        command[5] = str(path)
    failed = subprocess.run(command, input=data, capture_output=True)
    assert failed.returncode == 65 and expected in failed.stderr.decode(), failed.stderr
    assert name == "existing" or not Path(command[-1]).exists()
    assert not list(out.glob(".numivivo-celltypist-*"))
    failures.append({"case": name, "error": failed.stderr.decode()})
assert hashlib.sha256((out / "bundle/results.jsonl").read_bytes()).hexdigest() == receipt["resultsSHA256"]
report = dict(status="PASS", cells=len(rows), records=len(payload)//16,
              receipt=receipt, exactRetainedResults=True, failures=failures,
              biologicalAccuracyEstablished=False)
(out / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
