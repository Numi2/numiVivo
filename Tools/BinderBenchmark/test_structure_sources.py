#!/usr/bin/env python3
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import prepare_structure_sources as source


class StructureSourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.raw = b"HEADER test\r\nATOM source\r\nEND\r\n"
        (self.root / "input.pdb").write_bytes(self.raw)
        self.row = dict(candidateID="candidate-1", target="declared-target", sourceLabel="published model",
                        format="pdb", sha256=hashlib.sha256(self.raw).hexdigest(),
                        targetChainSequences={"A": "GG"}, interfacePlan={"targetChains": ["A"]},
                        sourcePath="input.pdb")
        self.manifest = {"schemaVersion": 1, "sources": [self.row]}

    def test_preserves_original_bytes(self):
        result = source.prepare(self.manifest, self.root)
        self.assertEqual(result["sources"][0]["contents"].encode(), self.raw)
        self.assertNotIn("sourcePath", result["sources"][0])

    def test_mismatch_rejected(self):
        self.row["sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "mismatch"):
            source.prepare(self.manifest, self.root)

    def test_duplicate_rejected(self):
        self.manifest["sources"].append(copy.deepcopy(self.row))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            source.prepare(self.manifest, self.root)

    def test_path_traversal_and_absolute_rejected(self):
        for path in ("../input.pdb", "/input.pdb", "./input.pdb", "folder//input.pdb", "folder\\input.pdb"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                source.local_source(self.root, path)

    def test_symlink_rejected(self):
        (self.root / "link.pdb").symlink_to(self.root / "input.pdb")
        self.row["sourcePath"] = "link.pdb"
        with self.assertRaisesRegex(ValueError, "symbolic"):
            source.prepare(self.manifest, self.root)

    def test_target_chain_reference_required(self):
        self.row["targetChainSequences"] = {"B": "GG"}
        with self.assertRaisesRegex(ValueError, "target chain"):
            source.prepare(self.manifest, self.root)

    def test_missing_and_unknown_fields_rejected(self):
        for field, value in (("contents", "invented"), ("outcome", True)):
            manifest = copy.deepcopy(self.manifest)
            manifest["sources"][0][field] = value
            with self.assertRaises(ValueError):
                source.prepare(manifest, self.root)

    def test_invalid_utf8_rejected(self):
        data = b"\xff\xfe"
        (self.root / "input.pdb").write_bytes(data)
        self.row["sha256"] = hashlib.sha256(data).hexdigest()
        with self.assertRaises(UnicodeDecodeError):
            source.prepare(self.manifest, self.root)

    def test_bounded_read(self):
        with self.assertRaises(ValueError):
            source.read_regular(self.root / "input.pdb", 2)

    def test_atomic_publication_and_no_overwrite(self):
        output = self.root / "sources.json"
        result = source.prepare(self.manifest, self.root)
        source.write_new(output, result)
        self.assertEqual(json.loads(output.read_bytes()), result)
        original = output.read_bytes()
        with self.assertRaises(FileExistsError):
            source.write_new(output, {"different": True})
        self.assertEqual(output.read_bytes(), original)
        self.assertEqual(list(self.root.glob(".binder-sources-*")), [])


if __name__ == "__main__":
    unittest.main()
