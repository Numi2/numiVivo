"""Independent FP64 matrix-exponential numerical reference; never measured data.

Read only geometry and authored plan. No native result/model is consumed.
The result represents the discrete graph model, not continuum-mesh convergence.
"""
import argparse,itertools
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from wetlab import read,write,sha

def reference(plan,geometry):
    p=read(plan);g=read(geometry);nodes=np.array(g['nodesMetres']);tets=np.array(g['tetrahedra']);n=len(tets)
    xyz=nodes[tets];centers=xyz.mean(axis=1);volumes=np.linalg.det(np.moveaxis(xyz[:,1:]-xyz[:,0,None],1,2))/6
    for c in p['cells']:volumes[c['tetrahedron']]-=c['volumeCubicMetres']
    volumes=np.r_[volumes,[c['volumeCubicMetres'] for c in p['cells']]]
    matrix=np.zeros((len(volumes),len(volumes)))
    def edge(i,j,k):
        matrix[i,i]-=k/volumes[i];matrix[i,j]+=k/volumes[i]
        matrix[j,j]-=k/volumes[j];matrix[j,i]+=k/volumes[j]
    for i,j in itertools.combinations(range(n),2):
        face=sorted(set(tets[i])&set(tets[j]))
        if len(face)==3:
            v=nodes[face];area=np.linalg.norm(np.cross(v[1]-v[0],v[2]-v[0]))/2
            edge(i,j,p['diffusionSquareMetresPerSecond']*area/np.linalg.norm(centers[i]-centers[j]))
    for i,c in enumerate(p['cells']):
        radius=(3*c['volumeCubicMetres']/(4*np.pi))**(1/3)
        edge(c['tetrahedron'],n+i,c['membranePermeabilityMetresPerSecond']*4*np.pi*radius**2)
    arms=[]
    for treated in (False,True):
        x=np.array(p['extracellularInitialMolPerM3']+p['intracellularInitialMolPerM3'],dtype=float);time=0;samples=[]
        for event in sorted(set(p['sampleTimesSeconds']+[e['timeSeconds'] for e in p['pulses']])):
            x=expm(matrix*(event-time))@x;time=event
            if treated:
                for pulse in p['pulses']:
                    if pulse['timeSeconds']==event:x[pulse['extracellularIndices']]+=pulse['concentrationIncrementMolPerM3']
            if event in p['sampleTimesSeconds']:samples.append({'timeSeconds':event,'concentrationsMolPerM3':x.tolist()})
        arms.append({'name':'intervention' if treated else 'control','samples':samples})
    return {'format':'wet-lab-spatial-observations/v1','specimenID':p['specimenID'],'geometrySHA256':sha(geometry),
            'units':'mol/m3','compartmentIDs':[f'ecs-{i}' for i in range(n)]+[c['id'] for c in p['cells']],
            'evidenceClass':'numerical-reference','provenance':'Independent SciPy FP64 matrix exponential on the authored graph. Numerical calibration only; no measured biological observations.',
            'arms':arms}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('plan',type=Path);p.add_argument('geometry',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    write(a.output,reference(a.plan,a.geometry))
