#!/usr/bin/env python3
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import fetch_structure_panel as panel


def csv_bytes(rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
    writer.writeheader(); writer.writerows(rows)
    return stream.getvalue().encode()


class StructurePanelTests(unittest.TestCase):
    def setUp(self):
        self.source = [{"uuid": f"{target}-{i}", "target": target, "sequence": "G" * (i + 1),
                        "adaptyv_binding": "binder", "ipsae_min_boltz2": "0.8"}
                       for target in panel.TARGETS for i in range(10)]
        self.selected = panel.select(self.source)
        self.model = b"data_synthetic\n# not a parsed protein\n"
        self.manifest = [{"uuid": r["uuid"], "target": r["target"], "kind": "predicted", "predictor": "boltz2",
                          "stoich": "1to1", "rel_path": f"designs/{r['target']}/{r['uuid']}/insilico/predicted_boltz2_1to1.cif",
                          "bytes": str(len(self.model)), "sha256": panel.sha256(self.model)} for r in self.selected]

    def test_selection_is_order_and_outcome_independent(self):
        changed = copy.deepcopy(self.source)[::-1]
        for r in changed:
            r["adaptyv_binding"] = "non_binder"; r["ipsae_min_boltz2"] = "0.01"
        self.assertEqual(panel.select(changed), self.selected)
        self.assertEqual(len(self.selected), 18)
        self.assertTrue(all(set(r) == {"uuid", "target", "sequence"} for r in self.selected))

    def test_missing_model_remains_explicit(self):
        result = panel.declaration(self.source, self.manifest[:-1])
        self.assertEqual(len(result["requested"]), 18)
        self.assertEqual(len(result["models"]), 17)
        self.assertEqual(result["unavailable"][0]["candidateID"], self.manifest[-1]["uuid"])

    def test_duplicate_candidates_rejected(self):
        with self.assertRaises(ValueError):
            panel.select(self.source + self.source[:1])

    def test_duplicate_models_rejected(self):
        with self.assertRaises(ValueError):
            panel.declaration(self.source, self.manifest + self.manifest[:1])

    def test_target_mismatch_rejected(self):
        self.manifest[0]["target"] = "other"
        with self.assertRaises(ValueError):
            panel.declaration(self.source, self.manifest)

    def test_unsafe_paths_rejected(self):
        for path in ("../x", "/designs/x", "designs/../predicted_boltz2_1to1.cif", "https://other/file", "designs/x/model.cif"):
            with self.subTest(path=path), self.assertRaises(ValueError):
                panel.source_path(path)

    def test_hash_and_capacity_rejected(self):
        for field, value in (("sha256", "bad"), ("bytes", "0"), ("bytes", str(panel.MAX_FILE_BYTES + 1))):
            rows = copy.deepcopy(self.manifest); rows[0][field] = value
            with self.assertRaises(ValueError):
                panel.declaration(self.source, rows)

    def test_ragged_and_duplicate_columns_rejected(self):
        for raw in (b"a,a\n1,2\n", b"a,b\n1\n", b"a\n1,2\n"):
            with self.assertRaises(ValueError):
                panel.rows(raw)

    def test_acquisition_preserves_bytes_and_refuses_overwrite(self):
        raw = csv_bytes(self.source)
        blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        metadata = {"structures_manifest.csv": csv_bytes(self.manifest),
                    "target_constructs.fasta": b">synthetic\nG\n", "INSILICO.md": b"synthetic documentation\n"}
        def get(url):
            if url == panel.URL: return raw
            for name, path in panel.METADATA.items():
                if url == panel.BASE + path: return metadata[name]
            self.assertTrue(url.startswith(panel.BASE + "data/designs/"))
            return self.model
        with tempfile.TemporaryDirectory() as temp, patch.object(panel, "BLOB", blob), patch.object(panel, "get", get):
            dest = Path(temp) / "new"
            receipt = panel.acquire(dest)
            self.assertEqual(receipt["availableCount"], 18)
            for path, digest in receipt["files"].items():
                self.assertEqual(panel.sha256((dest/path).read_bytes()), digest)
            with self.assertRaises(ValueError):
                panel.acquire(dest)
            self.assertEqual((dest/"source.csv").read_bytes(), raw)

    def test_failed_acquisition_publishes_nothing(self):
        with tempfile.TemporaryDirectory() as temp, patch.object(panel, "get", return_value=b"wrong source"):
            dest = Path(temp) / "new"
            with self.assertRaises(ValueError):
                panel.acquire(dest)
            self.assertFalse(dest.exists())


if __name__ == "__main__":
    unittest.main()
