import copy
import math
import unittest
from audit_endpoints import kinetic_check, whole_molecule_positions


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

    def test_whole_molecular_images_preserve_boundary_bonds(self):
        system=dict(particles=[{}, {}, {}],constraints=[dict(a=0,b=1),dict(a=1,b=2),dict(a=0,b=2)],nonbondedExceptions=[dict(a=0,b=2)])
        checkpoint=dict(positionsNM=[dict(x=.9375,y=.25,z=0),dict(x=.0625,y=.25,z=0),dict(x=.9375,y=.375,z=0)],
            periodicCell=dict(a=dict(x=1,y=0,z=0),b=dict(x=0,y=1,z=0),c=dict(x=0,y=0,z=1)))
        positions,images=whole_molecule_positions(system,checkpoint)
        self.assertEqual(positions,[[.9375,.25,0],[1.0625,.25,0],[.9375,.375,0]])
        self.assertEqual(images["shiftedParticles"],1);self.assertEqual(images["components"],1)
        # Arbitrary atom-wise periodic copies give the same relative geometry.
        checkpoint["positionsNM"][1]["x"]+=3
        self.assertEqual(whole_molecule_positions(system,checkpoint)[0],positions)

    def test_winding_cycle_and_ambiguous_half_box_refused(self):
        system=dict(particles=[{}, {}, {}],bonds=[dict(a=0,b=1),dict(a=1,b=2),dict(a=0,b=2)])
        checkpoint=dict(positionsNM=[dict(x=x,y=0,z=0) for x in (0,.4,.8)],
            periodicCell=dict(a=dict(x=1,y=0,z=0),b=dict(x=0,y=1,z=0),c=dict(x=0,y=0,z=1)))
        with self.assertRaises(ValueError):whole_molecule_positions(system,checkpoint)
        checkpoint["positionsNM"][1]["x"]=.5
        with self.assertRaises(ValueError):whole_molecule_positions(system,checkpoint)


if __name__=="__main__":unittest.main()
