#!/usr/bin/env python3
import unittest
import numpy as np
from sklearn.linear_model import Ridge
from test_integration_programs import Programs
from check_integration_programs import evaluate as original_evaluate,compare
from check_program_calibration import fit_within,evaluate


class Calibration(unittest.TestCase):
    def setUp(self):
        self.f=Programs();self.f.setUp()
    def tearDown(self):self.f.tearDown()

    def test_direct_centered_rows_in_every_fold(self):
        f=self.f
        for donor in range(4):
            selected=[i for i in range(8) if i//2!=donor];train=np.isin(f.codes,selected)
            codes=f.codes[train];x=f.x[train].copy();y=f.y[train].copy()
            weights=np.array([1/(len(selected)*sum(codes==c)) for c in codes])
            for k in selected:
                mask=codes==k;x[mask]-=x[mask].mean(axis=0);y[mask]-=y[mask].mean(axis=0)
            scale=np.sqrt(np.average(x*x,axis=0,weights=weights));scale[scale<1e-12]=1
            oracle=Ridge(alpha=1.,fit_intercept=False,solver='svd').fit(x/scale,y,sample_weight=weights)
            w,_,_=fit_within(f.stats,selected,1.)
            np.testing.assert_allclose(w,oracle.coef_.T/scale[:,None],rtol=1e-12,atol=1e-12)

    def test_training_library_offsets_do_not_change_gradient_fit(self):
        f=self.f;selected=list(range(6));w,_,_=fit_within(f.stats,selected,1.)
        changed={k:v.copy() for k,v in f.stats.items()}
        dx=np.arange(8)[:,None]*np.array([[2.,-3.,0.]])
        dy=np.arange(8)[:,None]*np.array([[5.,-7.]])
        changed['xx']+=np.einsum('si,sj->sij',f.stats['x'],dx)+np.einsum('si,sj->sij',dx,f.stats['x'])+np.einsum('si,sj->sij',dx,dx)
        changed['xy']+=np.einsum('si,sj->sij',f.stats['x'],dy)+np.einsum('si,sj->sij',dx,f.stats['y'])+np.einsum('si,sj->sij',dx,dy)
        changed['x']+=dx;changed['y']+=dy
        actual,_,_=fit_within(changed,selected,1.)
        np.testing.assert_allclose(w,actual,rtol=1e-10,atol=1e-10)

    def test_controls_and_held_out_target_independence(self):
        f=self.f;old=original_evaluate(f.stats,f.codes,f.ids,f.design,f.protocol)
        base=evaluate(old,f.ids,f.protocol);identity=evaluate(old,f.ids,f.protocol)
        self.assertEqual(base,identity)
        self.assertTrue(all(all(c['gates'].values()) for c in compare(base,identity,f.protocol)))
        erased=evaluate(old,f.ids,f.protocol,True)
        self.assertTrue(all(m['withinLibraryR2']==0 for fold in erased['folds'] for m in fold['measurements']))
        self.assertTrue(all(not c['gates']['meanPreserved'] for c in compare(base,erased,f.protocol)))
        changed={k:v.copy() for k,v in f.stats.items()};changed['y'][6:]+=100;changed['xy'][6:]+=200
        w,b,_=fit_within(f.stats,list(range(6)),1.);w2,b2,_=fit_within(changed,list(range(6)),1.)
        np.testing.assert_array_equal(w,w2);np.testing.assert_array_equal(b,b2)


if __name__=='__main__':unittest.main(verbosity=2)
