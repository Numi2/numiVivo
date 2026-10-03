import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from adapters import WetLabExperimentAdapter, adapter_for_family, adapter_for_run

class AdapterTests(unittest.TestCase):
    def test_boundary_is_abstract(self):
        with self.assertRaises(TypeError):WetLabExperimentAdapter()
    def test_owners_and_unsupported_family(self):
        self.assertEqual(adapter_for_family('cell-response').adapter_id,'rna-response/v1')
        self.assertEqual(adapter_for_family('spatial-tissue').adapter_id,'spatial-tissue/v1')
        with self.assertRaises(ValueError):adapter_for_family('molecular')
    def test_legacy_record_dispatch_and_version_rejection(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'registration.json';p.write_text(json.dumps({'assay':{'family':'cell-response'}}))
            self.assertEqual(adapter_for_run(d).adapter_id,'rna-response/v1')
            p.write_text(json.dumps({'family':'spatial-tissue','adapterID':'spatial-tissue/v2'}))
            with self.assertRaises(ValueError):adapter_for_run(d)
    def test_rna_delegates_existing_verifier(self):
        with patch('adapters.rna.verify',return_value={'status':'verified'}) as verify:
            result=adapter_for_family('cell-response').verify(Path('/record'),{'binary':Path('/native')})
            verify.assert_called_once_with(Path('/record'),Path('/native'))
            self.assertEqual(result,{'status':'verified'})

if __name__=='__main__':unittest.main()
