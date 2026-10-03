import copy
import tempfile
import unittest
from pathlib import Path
from wetlab import pairs, score_values, write, inventory, check_seal, FORMAT


class RecordTests(unittest.TestCase):
    def setUp(self):
        self.assay = {'cellGroup': 'B', 'minimumCellsPerArm': 1,
                      'controlCondition': 'ctrl', 'treatmentCondition': 'stim'}
        self.bulk = {'featureIDs': ['a', 'b'], 'groups': []}
        for d in range(3):
            for c in ('ctrl', 'stim'):
                self.bulk['groups'].append({'cellGroup': 'B', 'donorID': str(d),
                    'biologicalReplicateID': str(d), 'condition': c,
                    'sourceCellIndices': [len(self.bulk['groups'])]})

    def test_pairs_and_unsupported_population(self):
        self.assertEqual(pairs(self.bulk, self.assay)['1'], [2, 3])
        self.assay['cellGroup'] = 'T'
        with self.assertRaisesRegex(ValueError, 'three eligible'): pairs(self.bulk, self.assay)

    def test_duplicate_arms_are_not_silently_combined(self):
        self.bulk['groups'].append(copy.deepcopy(self.bulk['groups'][0]))
        with self.assertRaisesRegex(ValueError, 'three eligible'): pairs(self.bulk, self.assay)

    def test_aliased_biological_units_rejected(self):
        for g in self.bulk['groups']: g['biologicalReplicateID'] = 'same-unit'
        with self.assertRaisesRegex(ValueError, 'aliases'): pairs(self.bulk, self.assay)

    def test_shared_cells_rejected(self):
        self.bulk['groups'][1]['sourceCellIndices'] = [0]
        with self.assertRaisesRegex(ValueError, 'overlap'): pairs(self.bulk, self.assay)

    def test_scoring_uses_held_out_counts_and_seal_rejects_tampering(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); (root / 'prediction/source').mkdir(parents=True)
            (root / 'prediction/folds/000').mkdir(parents=True)
            self.bulk['matrix'] = {'rowOffsets': [0, 2, 4], 'featureIndices': [0, 1, 0, 1], 'counts': [1, 1, 1, 1]}
            control = [__import__('math').log1p(500000)] * 2
            write(root / 'registration.json', {'donor': '0', 'heldOutRow': 1, 'controlRow': 0, 'assay': self.assay})
            write(root / 'prediction/source/report.json', {'pseudobulk': self.bulk})
            write(root / 'prediction/folds/000/prediction.json', {'featureIDs': ['a', 'b'], 'predictions': [
                {'group': self.bulk['groups'][0], 'control': control, 'estimates': [
                    {'baseline': 'contextRidge', 'predictedResponse': [3, 4]},
                    {'baseline': 'noChange', 'predictedResponse': [0, 0]},
                    {'baseline': 'meanResponse', 'predictedResponse': [1, 1]}]}]})
            write(root / 'seal.json', {'format': FORMAT, 'files': inventory(root)})
            result = score_values(root)
            self.assertAlmostEqual(result['metrics'][0]['rmse'], (12.5)**0.5)
            self.assertEqual(result['metrics'][0]['mae'], 3.5)
            self.assertEqual(result['verdict'], 'does-not-beat-both-baselines')
            (root / 'prediction/source/report.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'Sealed artifact changed'): check_seal(root)


if __name__ == '__main__': unittest.main()
