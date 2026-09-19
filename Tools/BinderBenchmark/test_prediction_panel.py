#!/usr/bin/env python3
"""Offline contracts for the independent prediction-panel checker."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import check_prediction_panel as checker


class PredictionPanelTests(unittest.TestCase):
    def cif(self) -> str:
        headers = ('group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_comp_id',
                   'label_asym_id', 'label_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z',
                   'occupancy', 'B_iso_or_equiv', 'auth_asym_id', 'label_alt_id',
                   'pdbx_PDB_ins_code', 'pdbx_PDB_model_num')
        text = 'data_fixture\nloop_\n' + ''.join('_atom_site.' + h + '\n' for h in headers)
        serial = 0
        for chain, residue, names, y in (
                ('B', 'GLY', ['N', 'CA', 'C', 'O'], 0),
                ('T', 'ALA', ['N', 'CA', 'C', 'O', 'CB'], 4)):
            for i, name in enumerate(names):
                serial += 1
                element = name[0]
                text += f'ATOM {serial} {element} {name} {residue} {chain} 1 {i} {y} 0 1 0 {chain} . ? 1\n'
        return text

    def model(self) -> dict:
        return {'candidateID': 'test', 'target': 'target', 'published':
                {'target_construct_id': 'construct', 'predictor': 'fixture', 'seed': '1',
                 'resolvedPath': 'synthetic.cif'}}

    def test_fasta_rejects_duplicates_and_unidentified_sequence(self):
        self.assertEqual(checker.fasta('>a label\nAG\nC\n'), {'a': 'AGC'})
        for text in ('>a\nA\n>a\nG\n', 'A\n>a\nG\n'):
            with self.assertRaises(ValueError):
                checker.fasta(text)

    def test_explicit_units_and_radius_profile(self):
        p = checker.surface_plan(960)
        self.assertEqual(p['probeRadiusNM'], .14)
        self.assertEqual(p['radiiNM']['C'], .17)
        self.assertEqual(p['pointsPerAtom'], 960)
        self.assertNotIn('affinity', p)

    def test_unique_sequence_assignment_and_raw_hash(self):
        raw = self.cif()
        request, _ = checker.prepare(self.model(), {'sequence': 'G'},
                                     {'construct|protein_chain_1': 'A'}, raw, 96)
        self.assertEqual(request['source']['interfacePlan']['binderChains'], ['B'])
        self.assertEqual(request['source']['interfacePlan']['targetChains'], ['T'])
        self.assertEqual(request['source']['sha256'], hashlib.sha256(raw.encode()).hexdigest())
        self.assertEqual(request['identity']['sampleID'], 'published-seed-1')
        self.assertEqual(request['source']['contents'], raw)

    def test_ambiguous_chain_assignment_rejected(self):
        with self.assertRaises(ValueError):
            checker.prepare(self.model(), {'sequence': 'G'},
                            {'construct|protein_chain_1': 'G'}, self.cif(), 96)

    def test_unsupported_atom_population_rejected(self):
        raw = self.cif().replace('ATOM 1 N N ', 'ATOM 1 H N ')
        with self.assertRaises(ValueError):
            checker.prepare(self.model(), {'sequence': 'G'},
                            {'construct|protein_chain_1': 'A'}, raw, 96)

    def test_checker_rejects_area_tampering(self):
        request, structure = checker.prepare(self.model(), {'sequence': 'G'},
                                             {'construct|protein_chain_1': 'A'}, self.cif(), 96)
        report = {'surface': {'atoms': [{'atomIndex': i, 'complexAreaNM2': 1000,
                                       'isolatedAreaNM2': 1000} for i in range(9)]}}
        with self.assertRaisesRegex(ValueError, 'atomic SASA mismatch'):
            checker.compare(report, structure, request)

    def test_manifest_traversal_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'SOURCE.json').write_text(json.dumps({'files': {'../bad': '0'*64}}))
            with patch.object(checker.subprocess, 'run') as run:
                with self.assertRaisesRegex(ValueError, 'unsafe source manifest path'):
                    checker.run(root, root/'no-binary', root/'output', 96)
                run.assert_not_called()
            self.assertFalse((root/'output').exists())

    def test_manifest_hash_and_symlink_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root/'data').write_text('actual')
            (root/'SOURCE.json').write_text(json.dumps({'files': {'data': '0'*64}}))
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                checker.run(root, root/'no-binary', root/'output', 96)
            (root/'alias').symlink_to(root/'data')
            (root/'SOURCE.json').write_text(json.dumps({'files': {'alias': checker.digest(b'actual')}}))
            with self.assertRaisesRegex(ValueError, 'unsafe source manifest path'):
                checker.run(root, root/'no-binary', root/'output', 96)

    def test_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FileExistsError):
                checker.run(root, root/'no-binary', root, 96)
            file = root/'result.json'
            checker.write_json(file, {'a': 1})
            with self.assertRaises(FileExistsError):
                checker.write_json(file, {'a': 2})
            self.assertEqual(json.loads(file.read_text()), {'a': 1})


if __name__ == '__main__':
    unittest.main()
