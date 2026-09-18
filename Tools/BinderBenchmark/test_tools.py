#!/usr/bin/env python3
"""Offline regression for source retrieval and plan construction, not a live download."""
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FETCH = load("fetch_source")
PLAN = load("make_plan")


class SourceTests(unittest.TestCase):
    raw = b"uuid,target\na,MBP\n"

    def invoke(self, output: Path, blob: str | None = None, error: Exception | None = None):
        digest = hashlib.sha1(b"blob " + str(len(self.raw)).encode() + b"\0" + self.raw).hexdigest()
        with patch("sys.argv", ["fetch_source.py", str(output)]), \
             patch.object(FETCH, "BLOB", blob or digest), \
             patch.object(FETCH.urllib.request, "urlopen", side_effect=error,
                          return_value=io.BytesIO(self.raw)) as request, \
             redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            FETCH.main()
            request.assert_called_once_with(FETCH.URL, timeout=120)

    def test_verified_bytes_published_with_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "source"
            self.invoke(output)
            self.assertEqual((output / "source.csv").read_bytes(), self.raw)
            provenance = json.loads((output / "SOURCE.json").read_text())
            self.assertEqual(provenance["sha256"], hashlib.sha256(self.raw).hexdigest())
            self.assertEqual(provenance["revision"], FETCH.REVISION)
            self.assertEqual(provenance["bytes"], len(self.raw))
            for assay in ("adaptyv", "twist"):
                config = json.loads((output / f"import-{assay}.json").read_text())
                self.assertEqual(config["assay"], assay)
                self.assertEqual(config["targets"], FETCH.TARGETS)
            self.assertEqual([p.name for p in Path(tmp).iterdir()], ["source"])

    def test_bad_pin_leaves_no_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "source"
            with self.assertRaisesRegex(ValueError, "pinned source mismatch"):
                self.invoke(output, blob="0" * 40)
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_network_failure_leaves_no_artifact(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "source"
            with self.assertRaises(OSError):
                self.invoke(output, error=OSError("synthetic network error"))
            self.assertFalse(output.exists())
            self.assertEqual(list(Path(tmp).iterdir()), [])

    def test_existing_directory_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "source"
            output.mkdir()
            (output / "existing").write_text("keep")
            with self.assertRaises(SystemExit):
                self.invoke(output)
            self.assertEqual((output / "existing").read_text(), "keep")


class PlanTests(unittest.TestCase):
    def workspace(self, root: Path):
        bundle = root / "bundle"
        bundle.mkdir()
        (bundle / "import.json").write_text(json.dumps({
            "targets": ["A", "B", "C"], "features": ["ipsae_min_boltz2", "ipsae_min_ptxv2"]}))
        (bundle / "imported.json").write_text(json.dumps({
            "dataset": {"sourceSHA256": "a" * 64, "records": [{"outcome": "binder"}]}}))
        return bundle

    def invoke(self, bundle: Path, output: Path, tests: list[str]):
        arguments = ["make_plan.py", str(bundle), "--output", str(output)]
        for target in tests:
            arguments += ["--test-target", target]
        with patch("sys.argv", arguments), redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            PLAN.main()

    def test_plan_does_not_depend_on_test_outcomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = self.workspace(root)
            first, second = root / "first.json", root / "second.json"
            self.invoke(bundle, first, ["B"])
            imported = json.loads((bundle / "imported.json").read_text())
            imported["dataset"]["records"][0]["outcome"] = "nonBinder"
            (bundle / "imported.json").write_text(json.dumps(imported))
            self.invoke(bundle, second, ["B"])
            self.assertEqual(first.read_bytes(), second.read_bytes())
            result = json.loads(first.read_text())
            self.assertEqual(result["trainingTargets"], ["A", "C"])
            self.assertEqual(result["testTargets"], ["B"])

    def test_empty_training_and_duplicate_or_unknown_target_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = self.workspace(root)
            output = root / "plan.json"
            for tests in (["A", "B", "C"], ["B", "B"], ["unknown"]):
                with self.assertRaises(SystemExit):
                    self.invoke(bundle, output, tests)
                self.assertFalse(output.exists())

    def test_existing_plan_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bundle = self.workspace(root)
            output = root / "plan.json"
            self.invoke(bundle, output, ["B"])
            original = output.read_bytes()
            with self.assertRaises(FileExistsError):
                self.invoke(bundle, output, ["A"])
            self.assertEqual(output.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
