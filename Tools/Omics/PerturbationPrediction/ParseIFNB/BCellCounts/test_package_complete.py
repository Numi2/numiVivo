"""Archive integrity tests are software checks, not complete cohort evidence."""
import fcntl
import hashlib
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from package_complete import inventory, package, verify_contents

class ArchiveIntegrity(unittest.TestCase):
    def test_archive_members(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            archive = root / "evidence.tar.gz"
            data = b"retained counts"
            expected = {"counts": dict(bytes=len(data), SHA256=hashlib.sha256(data).hexdigest())}
            def create(names):
                with tarfile.open(archive, "w:gz") as stream:
                    for name in names:
                        member = tarfile.TarInfo(name)
                        member.size = len(data)
                        stream.addfile(member, io.BytesIO(data))
            create(["counts"])
            verify_contents(archive, expected)
            for names in ([], ["unexpected"], ["counts", "counts"]):
                create(names)
                with self.assertRaises(ValueError):
                    verify_contents(archive, expected)
            create(["counts"])
            with self.assertRaises(ValueError):
                verify_contents(archive, {"counts": dict(bytes=len(data), SHA256="0" * 64)})

    def test_live_controller_and_incomplete_state(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "study"
            root.mkdir()
            lockpath = root / "controller.lock"
            lockpath.touch()
            output = parent / "complete.tar.gz"
            with lockpath.open("rb") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                with self.assertRaises(BlockingIOError):
                    package(root, parent / "manifest", output)
            with self.assertRaises(FileNotFoundError):
                package(root, parent / "manifest", output)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_name(output.name + ".partial").exists())

    def test_inventory_rejects_symlinks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "link").symlink_to(root.parent)
            with self.assertRaises(ValueError):
                inventory(root)

if __name__ == "__main__":
    unittest.main()

