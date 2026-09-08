"""Reference parameter controls; no simulated trajectory is used as an oracle."""
import copy
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import openmm as mm
from audit_mass_constraints import compare_parameters
from run_campaign import digest, write

TOOLS=Path(__file__).parent


def fixture():
    reference=mm.System()
    for mass in (16,1,1):reference.addParticle(mass)
    constraints=[dict(a=0,b=1,distanceNM=.1),dict(a=0,b=2,distanceNM=.1),dict(a=1,b=2,distanceNM=math.sqrt(.02))]
    for c in constraints:reference.addConstraint(c["a"],c["b"],c["distanceNM"])
    native=dict(particles=[dict(index=i,massDa=m) for i,m in enumerate((16,1,1))],constraints=constraints)
    return native,reference


class MassConstraintTests(unittest.TestCase):
    def test_exact_parameters_allow_reordered_symmetric_endpoints(self):
        native,reference=fixture()
        native["constraints"]=[dict(a=c["b"],b=c["a"],distanceNM=c["distanceNM"]) for c in reversed(native["constraints"])]
        result=compare_parameters(native,reference)
        self.assertTrue(result["passed"])
        self.assertEqual(result["particleCount"],3)
        self.assertEqual(result["maximumMassDifferenceDa"],0)

    def test_changed_mass_is_not_hidden_by_matching_constraints(self):
        native,reference=fixture();native["particles"][0]["massDa"]+=1e-12
        result=compare_parameters(native,reference)
        self.assertFalse(result["passed"])
        self.assertTrue(result["constraintTargetsEqual"])
        self.assertEqual(result["firstMismatchedMassIndices"],[0])

    def test_changed_missing_and_duplicate_constraints_refused(self):
        for variant in ("target","missing","duplicate"):
            native,reference=fixture()
            if variant=="target":native["constraints"][0]["distanceNM"]*=10
            elif variant=="missing":native["constraints"].pop()
            else:native["constraints"].append(copy.deepcopy(native["constraints"][0]))
            with self.subTest(variant=variant):
                result=compare_parameters(native,reference)
                self.assertFalse(result["passed"])
                self.assertTrue(result["particleMassesEqual"])

    def test_invalid_counts_indices_and_nonfinite_values_refused(self):
        mutations=[lambda n:n["particles"].pop(),lambda n:n["particles"][1].update(index=True),
            lambda n:n["particles"][0].update(massDa=math.nan),lambda n:n["constraints"][0].update(a=-1),
            lambda n:n["constraints"][0].update(distanceNM=math.inf)]
        for i,mutate in enumerate(mutations):
            native,reference=fixture();mutate(native)
            with self.subTest(mutation=i),self.assertRaises(ValueError):compare_parameters(native,reference)

    def test_cli_retains_input_failure_and_rejects_changed_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);case=root/"fixture";case.mkdir();native,reference=fixture()
            (case/"reference-system.xml").write_text(mm.XmlSerializer.serialize(reference))
            xml_hash=digest(case/"reference-system.xml")
            request=dict(identifier="fixture",system=native,referenceProvenance=dict(serializedSystemSHA256=xml_hash))
            write(case/"request.json",request)
            failure=dict(identifier="original-water",status="preparation-failed",error="retained fixture failure")
            write(root/"manifest.json",dict(schema="numivivo.org/md-reference-campaign/v1",cases=[
                dict(identifier="fixture",status="prepared",requestSHA256=digest(case/"request.json"),serializedSystemSHA256=xml_hash),failure]))
            def run(name):
                path=root/name
                result=subprocess.run([sys.executable,str(TOOLS/"audit_mass_constraints.py"),"--references",str(root),"--out",str(path)],capture_output=True,text=True)
                self.assertEqual(result.returncode,1,result.stderr)
                return json.loads(path.read_text())
            result=run("audit.json")
            self.assertTrue(result["preparedPassed"]);self.assertFalse(result["passed"])
            self.assertEqual(result["cases"][-1]["error"],failure["error"])
            request["system"]["particles"][0]["massDa"]=160
            (case/"request.json").write_text(json.dumps(request))
            result=run("altered-audit.json")
            self.assertFalse(result["preparedPassed"])
            self.assertIn("changed native request",result["cases"][0]["error"])


if __name__=="__main__":unittest.main()
