#!/usr/bin/env python3
"""Decision instrumentation qualification with known corruptions and signs."""
import numpy as np
from design_objective import objective,utility,selection_score,calibrated_controls,acquisition_backtest

def check():
    ix=objective(['A','B'],{'genes':['A'],'preserveGenes':['B'],'penalty':2})
    score,benefit,outside=utility(np.array([[-3.,0.],[-1.,0.]]),np.array([[0.,2.],[0.,0.]]),ix)
    assert score.tolist()==[-1.,1.] and benefit.tolist()==[3.,1.]
    try:utility(np.zeros((2,2)),None,ix);raise AssertionError('missing receiver accepted')
    except ValueError:pass
    result=selection_score(['A','B','no-intervention'],[2,1,0],[-1,3,0],[0,1,0])
    assert result['selected']=='A' and result['regret']==4 and result['gainOverTrainingMean']==-4
    assert result['gainOverNoIntervention']==-1 and not result['reliableWinner']
    calibration=calibrated_controls([[-2,1],[1,-3],[.2,.4]])
    assert calibration['status']=='passed'
    assert calibration['controls']['identity-positive']['meanProgramRegret']==0
    assert calibration['controls']['opposite-response']['meanProgramRegret']>0
    assert not acquisition_backtest([[1,0,3],[0,2,2]],[0,1,2],[10,3,2],1)['activeLearningBenefitEstablished']
    print('PASS: objective signs, preservation availability, failed selection regret, corruption calibration, acquisition boundary')

if __name__=='__main__':check()
