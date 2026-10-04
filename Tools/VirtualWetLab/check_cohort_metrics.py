#!/usr/bin/env python3
"""Positive/negative controls for cohort diagnostics; synthetic software checks."""
import numpy as np
from train_cohort_recovery import measure
rows=[{'source':'s','context':'c','target':t} for t in ('A','B')]
a={'context':np.zeros((2,3)), 'observed':np.array([[1.,-1.,1000.],[-1.,1.,-1000.]]),'mask':np.array([[1,1,0],[1,1,0]])}
correct=measure(a['observed'],a,rows);wrong=measure(-a['observed'],a,rows);flat=measure(np.zeros((2,3)),a,rows)
assert correct['equalGroupRMSE']==0 and correct['identityAccuracy']==1 and correct['directionAccuracy']==1
assert wrong['equalGroupRMSE']==2 and wrong['identityAccuracy']==0 and wrong['directionAccuracy']==0
assert flat['equalGroupRMSE']==1
corrupt=a['observed'].copy();corrupt[:,2]=1e9
assert measure(corrupt,a,rows)['equalGroupRMSE']==0
print('Masked feature, direction and target-identity positive/negative controls passed; no biological evidence.')
