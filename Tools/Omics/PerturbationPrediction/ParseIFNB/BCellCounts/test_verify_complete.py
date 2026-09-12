"""Admission regressions only; a completed real cohort must test the success path."""
import json
from pathlib import Path
import tempfile
import unittest
from verify_complete import inside, review


class TerminalAdmission(unittest.TestCase):
    def test_missing_terminal(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                review(Path(d), Path(d) / 'manifest.json')

    def test_forged_terminal_cannot_replace_preparation(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            terminal = dict(status='passed-all-native-ingestions-and-replays',
                            cells=72446, records=124909573,
                            predictionFitted=False, predictionScored=False)
            (root / 'complete.json').write_text(json.dumps(terminal))
            (root / 'payload').write_text('altered')
            (root / 'manifest.json').write_text(json.dumps({'members': {'payload': '0' * 64}}))
            with self.assertRaisesRegex(ValueError, 'preparation drift'):
                review(root, root / 'manifest.json')
            terminal['cells'] = 72445
            (root / 'complete.json').write_text(json.dumps(terminal))
            with self.assertRaisesRegex(ValueError, 'terminal declaration'):
                review(root, root / 'manifest.json')

    def test_paths_cannot_escape(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            with self.assertRaisesRegex(ValueError, 'outside study'):
                inside(root, '../another-study')
            (root / 'escape').symlink_to(root.parent)
            with self.assertRaisesRegex(ValueError, 'outside study'):
                inside(root, 'escape/another-study')


if __name__ == '__main__':
    unittest.main()
