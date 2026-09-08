import json
from pathlib import Path
import unittest
import numpy as np
from ensemble_statistics import kinetic_assessment, slope_assessment, effective_samples


class EnsembleStatisticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.policy=json.loads(Path(__file__).with_name("ensemble_policy.json").read_text())
        cls.policy["bootstrapReplicates"]=400
        cls.rng=np.random.default_rng(879123)
        cls.dof=5370;cls.temperature=300;cls.k=cls.policy["boltzmannKJPerMolK"]
        cls.kinetic=cls.rng.gamma(cls.dof/2,cls.k*cls.temperature,size=(3,4000))

    def assess(self,values):
        return kinetic_assessment(values,self.dof,self.temperature,self.policy)

    def test_canonical_gamma_control_passes(self):
        result=self.assess(self.kinetic)
        self.assertEqual(result["outcome"],"passed",result)

    def test_wrong_temperature_is_rejected(self):
        result=self.assess(self.kinetic*1.02)
        self.assertEqual(result["outcome"],"failed",result)
        self.assertIn("kinetic mean misses the prescribed temperature",result["failureReasons"])

    def test_suppressed_fluctuations_are_rejected(self):
        expected=self.dof/2*self.k*self.temperature
        result=self.assess(expected+0.25*(self.kinetic-self.kinetic.mean()))
        self.assertEqual(result["outcome"],"failed",result)

    def test_matching_moments_do_not_hide_wrong_distribution_shape(self):
        rng=np.random.default_rng(27488);expected=self.dof/2*self.k*self.temperature
        sigma=np.sqrt(self.dof/2)*self.k*self.temperature
        mixture=rng.choice([-1,1],size=(3,4000))*np.sqrt(0.96)+rng.normal(0,0.2,size=(3,4000))
        result=self.assess(expected+sigma*mixture)
        self.assertEqual(result["outcome"],"failed",result)
        self.assertTrue(any("shape rejected" in reason for reason in result["failureReasons"]))

    def test_constant_and_repeated_samples_cannot_certify_sampling(self):
        constant=np.full((3,4000),self.dof/2*self.k*self.temperature)
        self.assertEqual(self.assess(constant)["outcome"],"inconclusive")
        repeated=np.repeat(self.kinetic[:,:4],1000,axis=1)
        self.assertLess(effective_samples(repeated[0]),40)
        self.assertEqual(self.assess(repeated)["outcome"],"inconclusive")

    def test_drift_is_retained_as_insufficient_stationarity(self):
        values=self.kinetic+np.linspace(-600,600,4000)
        result=self.assess(values)
        self.assertEqual(result["outcome"],"inconclusive",result)
        self.assertTrue(any("drift" in reason for reason in result["inconclusiveReasons"]))

    def test_declared_temperature_ratio_is_recovered(self):
        rng=np.random.default_rng(6271)
        low=rng.gamma(1000,self.k*300,size=(3,2000))
        high=rng.gamma(1000,self.k*305,size=(3,2000))
        result=slope_assessment(low,high,300,305,self.policy)
        self.assertEqual(result["outcome"],"passed",result)

    def test_unchanged_distribution_at_two_declared_temperatures_fails(self):
        rng=np.random.default_rng(7271)
        low=rng.gamma(1000,self.k*300,size=(3,2000))
        high=rng.gamma(1000,self.k*300,size=(3,2000))
        result=slope_assessment(low,high,300,305,self.policy)
        self.assertEqual(result["outcome"],"failed",result)

    def test_no_overlap_is_inconclusive(self):
        rng=np.random.default_rng(23)
        low=rng.normal(0,1,size=(3,1000));high=rng.normal(10000,1,size=(3,1000))
        result=slope_assessment(low,high,300,305,self.policy)
        self.assertEqual(result["outcome"],"inconclusive",result)

    def test_nonfinite_data_and_partial_blocks_are_rejected(self):
        bad=self.kinetic.copy();bad[0,0]=np.nan
        with self.assertRaises(ValueError):self.assess(bad)
        with self.assertRaises(ValueError):self.assess(self.kinetic[:,:3999])


if __name__=="__main__":unittest.main()
