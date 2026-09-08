import copy
import math
import unittest
from audit_endpoints import kinetic_check


class EndpointEnergyTests(unittest.TestCase):
    def setUp(self):
        self.system=dict(particles=[dict(index=0,role="atom",massDa=2),dict(index=1,role="atom",massDa=3)])
        self.checkpoint=dict(velocitiesNMPerPS=[dict(x=1,y=2,z=3),dict(x=-2,y=0,z=1)])

    def test_independent_mass_and_units_accounting(self):
        result=kinetic_check(self.system,self.checkpoint,dict(kineticEnergyKJPerMol=21.5))
        self.assertTrue(result["passed"]);self.assertEqual(result["referenceKineticEnergyKJPerMol"],21.5)

    def test_plausible_but_wrong_reported_energy_rejected(self):
        self.assertFalse(kinetic_check(self.system,self.checkpoint,dict(kineticEnergyKJPerMol=21.6))["passed"])
        self.assertFalse(kinetic_check(self.system,self.checkpoint,dict(kineticEnergyKJPerMol=43))["passed"])

    def test_invalid_particle_state_refused(self):
        for key,value in (("massDa",0),("role","virtualSite"),("index",1)):
            system=copy.deepcopy(self.system);system["particles"][0][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):kinetic_check(system,self.checkpoint,dict(kineticEnergyKJPerMol=21.5))
        with self.assertRaises(ValueError):kinetic_check(self.system,self.checkpoint,dict(kineticEnergyKJPerMol=math.nan))


if __name__=="__main__":unittest.main()
