#!/usr/bin/env python3
"""Identification and invariance checks for the protected regression experiment."""
import unittest
import numpy as np
from fit_protected_regression import fit

class Checks(unittest.TestCase):
 def test_single_stratum_matches_categorical_ridge(self):
  mass=np.array([[2.,3.,5.]]);total=np.array([[[4.,2.],[9.,-3.],[0.,10.]]])
  actual,_,_=fit(mass,total);m=mass[0];s=total[0]
  center=(s/(m+1)[:,None]).sum(axis=0)/(m/(m+1)).sum()
  np.testing.assert_allclose(actual,(s-m[:,None]*center)/(m+1)[:,None],rtol=1e-12,atol=1e-12)
 def test_perfect_stratum_donor_confounding_does_not_invent_effect(self):
  mass=np.array([[100.,0.],[0.,200.]])
  total=np.array([[[900.,200.],[0.,0.]],[[0.,0.],[-1400.,400.]]])
  actual,_,_=fit(mass,total);np.testing.assert_allclose(actual,0,atol=1e-12)
 def test_stratum_shifts_do_not_change_donor_effect(self):
  mass=np.array([[8.,4.,2.],[5.,9.,10.],[1.,20.,6.]])
  means=np.arange(18,dtype=float).reshape(3,3,2)/7
  total=means*mass[:,:,None];first,_,a=fit(mass,total)
  offsets=np.array([[20.,-3.],[-11.,4.],[9.,-18.]])
  second,_,b=fit(mass,total+mass[:,:,None]*offsets[:,None,:])
  np.testing.assert_allclose(first,second,rtol=1e-11,atol=1e-11)
  np.testing.assert_allclose(np.array(b['intercepts'])-a['intercepts'],offsets,rtol=1e-11,atol=1e-11)

if __name__=='__main__':unittest.main()
