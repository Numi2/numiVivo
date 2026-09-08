"""Evidence-copy controls; no native runtime or statistical result is changed."""
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import verify_evidence as verifier


class EvidenceVerificationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.directory = self.root / "retained-failure"
        self.directory.mkdir()
        self.report = self.directory / "report.json"
        self.report.write_bytes(b'{"passed":false}\n')
        self.link = self.root / "dangling-link"
        self.link.symlink_to("absent-rejection-output")
        self.manifest = dict(schema=verifier.SCHEMA, exclusions=["environment"], files=[
            dict(path="retained-failure/report.json", type="file", bytes=self.report.stat().st_size,
                 sha256=hashlib.sha256(self.report.read_bytes()).hexdigest()),
            dict(path="dangling-link", type="symlink", target="absent-rejection-output")])
        self.manifest_path = self.root / "seal.json"
        self.save_manifest()

    def save_manifest(self):
        self.manifest_path.write_text(json.dumps(self.manifest))
        self.expected = hashlib.sha256(self.manifest_path.read_bytes()).hexdigest()

    def verify(self):
        return verifier.verify(self.root, self.manifest_path, self.expected)

    def test_verified_copy_preserves_failed_result_and_dangling_link(self):
        result = self.verify()
        self.assertTrue(result["passed"])
        self.assertEqual((result["regularFiles"], result["symlinks"]), (1, 1))
        self.assertFalse(json.loads(self.report.read_text())["passed"])
        self.assertIn("unlisted files are not verified", result["scope"])

    def test_changed_bytes_size_and_missing_file_each_fail(self):
        original = self.report.read_bytes()
        for data in (b"X" * len(original), original + b"\n", None):
            with self.subTest(data=data):
                if data is None:
                    self.report.unlink()
                else:
                    self.report.write_bytes(data)
                result = self.verify()
                self.assertFalse(result["passed"])
                self.assertEqual(len(result["failures"]), 1)

    def test_changed_pin_rejects_self_consistent_replacement_manifest(self):
        self.manifest["files"].pop()
        self.manifest_path.write_text(json.dumps(self.manifest))
        with self.assertRaisesRegex(ValueError, "pinned identity"):
            self.verify()

    def test_link_target_and_type_substitutions_are_rejected(self):
        self.link.unlink()
        self.link.symlink_to("different-absent-target")
        self.assertFalse(self.verify()["passed"])
        self.link.unlink()
        self.link.write_text("absent-rejection-output")
        self.assertFalse(self.verify()["passed"])
        self.report.unlink()
        self.report.symlink_to(self.manifest_path)
        self.assertEqual(len(self.verify()["failures"]), 2)

    def test_substituted_parent_link_is_not_followed(self):
        renamed = self.root / "moved-directory"
        self.directory.rename(renamed)
        self.directory.symlink_to(renamed, target_is_directory=True)
        self.assertFalse(self.verify()["passed"])

    def test_unsafe_duplicate_overlapping_and_malformed_entries_rejected(self):
        for name in ("/outside", "../outside", "a/../outside", "a//b", "a/./b", "a\\b", ""):
            candidate = copy.deepcopy(self.manifest)
            candidate["files"][0]["path"] = name
            with self.subTest(path=name), self.assertRaises(ValueError):
                verifier.validate_manifest(json.dumps(candidate))
        variants = []
        duplicate = copy.deepcopy(self.manifest)
        duplicate["files"].append(duplicate["files"][0])
        variants.append(duplicate)
        overlap = copy.deepcopy(self.manifest)
        overlap["files"][1]["path"] = "retained-failure"
        variants.append(overlap)
        for key, value in (("bytes", True), ("bytes", -1), ("sha256", "bad"), ("type", "directory")):
            candidate = copy.deepcopy(self.manifest)
            candidate["files"][0][key] = value
            variants.append(candidate)
        variants.append(dict(schema=verifier.SCHEMA, files=[]))
        for candidate in variants:
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                verifier.validate_manifest(json.dumps(candidate))
        with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
            verifier.validate_manifest('{"schema":"a","schema":"b"}')

    def test_manifest_changed_during_check_fails(self):
        original = verifier.verify_entry
        def mutate(*args):
            original(*args)
            self.manifest_path.write_text("{}")
        with patch.object(verifier, "verify_entry", side_effect=mutate):
            self.assertFalse(self.verify()["passed"])

    def test_file_changed_during_hash_is_rejected(self):
        hasher = hashlib.sha256()
        before = self.report.stat()
        report = self.report
        class MutatingHash:
            def update(self, data):
                hasher.update(data)
                os.utime(report, ns=(before.st_atime_ns, before.st_mtime_ns + 1_000_000_000))
            def hexdigest(self):
                return hasher.hexdigest()
        fd = os.open(self.root, os.O_RDONLY | os.O_DIRECTORY)
        try:
            with patch.object(verifier.hashlib, "sha256", return_value=MutatingHash()):
                with self.assertRaisesRegex(ValueError, "changed while reading"):
                    verifier.verify_entry(fd, self.manifest["files"][0])
        finally:
            os.close(fd)

    def test_cli_does_not_overwrite_an_existing_verification(self):
        out = self.root / "verification.json"
        command = [sys.executable, str(Path(verifier.__file__)), "--root", str(self.root),
                   "--manifest", str(self.manifest_path), "--expected-sha256", self.expected,
                   "--out", str(out)]
        first = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(first.returncode, 0, first.stderr)
        original = out.read_bytes()
        self.report.write_text("altered")
        second = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(second.returncode, 0)
        self.assertEqual(out.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
