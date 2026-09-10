#!/usr/bin/env python3
"""Numerical and missing-class controls for the annotation diagnostic."""
import unittest
import numpy as np
from check_annotation_retention import fit,accumulate,measures,training_stats
from check_annotation_retention_oracle import direct_fit,direct_confusion


def sufficient(x,y,c):
 n=np.bincount(y,minlength=c);sx=np.zeros((c,x.shape[1]));xx=np.zeros((c,x.shape[1],x.shape[1]))
 for k in range(c):
  z=x[y==k];sx[k]=z.sum(axis=0);xx[k]=z.T@z
 return n,sx,xx

class Checks(unittest.TestCase):
 def test_imbalanced_fit_matches_direct_svd(self):
  rng=np.random.default_rng(19);x=rng.normal(size=(113,5));y=np.array([0]*100+[2]*10+[4]*3)
  model=fit(*sufficient(x,y,6));classes,w,b=direct_fit(x,y)
  np.testing.assert_allclose(model['weight'],w,rtol=1e-11,atol=1e-11);np.testing.assert_allclose(model['intercept'],b,rtol=1e-11,atol=1e-11)
  q=rng.normal(size=(37,5));truth=rng.integers(0,6,37);actual=np.zeros((6,6));accumulate(actual,q,truth,model)
  np.testing.assert_array_equal(actual,direct_confusion(q,truth,classes,w,b,6))
 def test_erasure_uniform_and_unseen_class_zero_recall(self):
  x=np.zeros((12,3));y=np.array([0]*10+[2]*2);model=fit(*sufficient(x,y,4))
  result=np.zeros((4,4));accumulate(result,np.zeros((3,3)),np.array([0,1,2]),model)
  expected=np.array([[.5,0,.5,0],[.5,0,.5,0],[.5,0,.5,0],[0,0,0,0]])
  np.testing.assert_array_equal(result,expected)
  self.assertEqual(measures(result,np.array([1,1,1,0]))[1]['recall'],0)
 def test_balancing_unchanged_by_repeating_one_class(self):
  x=np.array([[1.,2.],[2.,4.],[4.,3.],[6.,2.]]);y=np.array([0,0,1,1]);a=fit(*sufficient(x,y,2))
  b=fit(*sufficient(np.concatenate([x,x[:2],x[:2]]),np.r_[y,0,0,0,0],2))
  np.testing.assert_allclose(a['weight'],b['weight'],rtol=1e-12,atol=1e-12)
 def test_held_donor_changes_cannot_change_training_inputs(self):
  n=np.array([[[12,20],[5,6],[18,22]]]);sx=np.arange(12,dtype=float).reshape(1,3,2,2);xx=np.ones((1,3,2,2,2))
  original=training_stats(n,sx,xx,0,1)
  n[0,1]=999;sx[0,1]=1e12;xx[0,1]=1e24
  for a,b in zip(original,training_stats(n,sx,xx,0,1)):np.testing.assert_array_equal(a,b)
 def test_no_training_not_a_success(self):
  self.assertIsNone(fit(*sufficient(np.empty((0,3)),np.array([],dtype=int),4)))

if __name__=='__main__':unittest.main()
