"""End-to-end refusal checks for incomplete or altered evidence matrices."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from ensemble_campaign import read, write, digest, analysis_environment

TOOLS=Path(__file__).parent


class EnsembleAnalysisTests(unittest.TestCase):
    def campaign(self,root):
        result=dict(schema="numivivo.org/md-ensemble-campaign/v1",policy=read(TOOLS/"ensemble_policy.json"),runs=[],
            inputFailures=[dict(identifier="original-water",status="preparation-failed",error="retained fixture failure")],
            excludedPreparedCases=[dict(identifier="protein",status="prepared")],identities=dict(analysisEnvironment=analysis_environment(),
            tools={n:digest(TOOLS/n) for n in ("analyze_ensemble.py","ensemble_campaign.py","ensemble_statistics.py","run_campaign.py")}))
        for name,key,value in (("policy.json","policySHA256",result["policy"]),("binary-manifest.json","binaryManifestSHA256",{}),
                               ("parent-manifest.json","parentManifestSHA256",{})):
            write(root/name,value);result["identities"][key]=digest(root/name)
        return result

    def test_missing_matrix_cannot_pass_or_erase_parent_failures(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);campaign=self.campaign(root);write(root/"campaign.json",campaign)
            result=subprocess.run([sys.executable,str(TOOLS/"analyze_ensemble.py"),"--campaign",str(root),"--out",str(root/"analysis.json")],capture_output=True,text=True)
            self.assertEqual(result.returncode,1,result.stderr)
            report=read(root/"analysis.json");self.assertFalse(report["passed"])
            self.assertEqual(report["preparedOutcome"],"failed")
            self.assertEqual(report["inputFailures"],campaign["inputFailures"])
            self.assertEqual(report["excludedPreparedCases"],campaign["excludedPreparedCases"])
            self.assertEqual(len(report["cases"]),2)
            self.assertTrue(all(r["outcome"]=="failed" and "error" in r for r in report["cases"]))

    def test_changed_implementation_or_duplicate_cell_refused(self):
        for corruption in ("implementation","duplicate","policy"):
            with self.subTest(corruption=corruption),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);campaign=self.campaign(root)
                if corruption=="implementation":campaign["identities"]["tools"]["ensemble_statistics.py"]="changed"
                elif corruption=="policy":(root/"policy.json").write_text("{}")
                else:
                    row=dict(timeStepPS=.001,temperatureK=300,seed=1729,directory="run-00")
                    campaign["runs"]=[row,row]
                write(root/"campaign.json",campaign)
                result=subprocess.run([sys.executable,str(TOOLS/"analyze_ensemble.py"),"--campaign",str(root),"--out",str(root/"analysis.json")],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertFalse((root/"analysis.json").exists())


if __name__=="__main__":unittest.main()
