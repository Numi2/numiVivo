from pathlib import Path
import json,subprocess,copy,math
s=Path(__file__).parent;out=s/'tests';out.mkdir(exist_ok=False)
a={'featureIDs':['g'],'trainingDonorIDs':['a1','a2','b1','b2','b3'],'trainingStudies':['A','A','B','B','B'],'controls':[[1]]*5,'treated':[[2],[2],[4],[4],[4]],'queryStudy':'C','queryDonorIDs':['c1'],'queryControls':[[5]]}
records=[]
def run(name,a,passed):
 p=out/(name+'.json');dest=out/(name+'-output.json');p.write_text(json.dumps(a));r=subprocess.run([str(s/'study-shrinkage'),str(p),str(dest)],capture_output=True,text=True);assert (r.returncode==0)==passed
 if not passed:assert not dest.exists()
 records.append({'case':name,'expectedSuccess':passed,'returnCode':r.returncode});return json.loads(dest.read_text()) if passed else None
b=run('equal-study-unequal-donor-weight',a,True);assert math.isclose(b['predictedTreated'][0][0],7,abs_tol=1e-12);assert math.isclose(sum(b['model']['donorWeights'][:2]),.5,abs_tol=1e-12);assert math.isclose(sum(b['model']['donorWeights'][2:]),.5,abs_tol=1e-12);assert b['model'].get('selectedPenalty') is None
for name,modify in [('query-study-in-training',lambda z:z.update(queryStudy='A')),('duplicate-training-donor',lambda z:z['trainingDonorIDs'].__setitem__(1,'a1')),('query-donor-in-training',lambda z:z.update(queryDonorIDs=['a1'])),('negative-log-expression',lambda z:z.update(queryControls=[[-1]])),('wrong-query-width',lambda z:z.update(queryControls=[[1,2]])),('one-training-study',lambda z:z.update(trainingStudies=['A']*5))]:
 z=copy.deepcopy(a);modify(z);run(name,z,False)
(s/'tests.json').write_text(json.dumps({'status':'pass','scope':'Software edge cases; biological evidence comes from the real three-study experiment.','cases':records},indent=2)+'\n');print('PASS',len(records))
