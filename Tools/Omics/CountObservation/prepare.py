"""Freeze original control-cell numerical cases, without reading treated values."""
from pathlib import Path
import hashlib,json,sys,time
import numpy as np
import anndata as ad
from scipy import sparse
study=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name,obj):(out/name).write_text(json.dumps(obj,sort_keys=True,indent=2,allow_nan=False)+'\n')
paths=['handoff/B-source-codes-RNA.h5ad','handoff/reference-groups.json','handoff/reference-pseudobulk.npz','prediction-inputs/cohort.json','prediction-inputs/panel.json']
protocol={'createdUnix':time.time(),'purpose':'Conditional numerical qualification of native NB2 count observation posterior; not a perturbation prediction or biological coverage test','source':str(study),'sources':{p:sha(study/p) for p in paths},'selection':'16 frozen panel gene IDs with lexicographically smallest SHA256 digests; all 62 admitted control donors','dispersionSensitivity':[0,1,10],'gammaShape':0.5,'gammaRatePerCPM':0.001,'priorInterpretation':'Explicit sensitivity prior with mean 500 CPM; not a learned prior or qualified default','futureSampling':'20 independent same-condition cells, each planned for 1000 full-RNA UMIs; no observed treated depth','treatedCountValuesUsed':False,'GSE181897KnownDevelopmentData':True,'nativeReferenceTolerance':{'relativeMoments':2e-7,'absoluteLogQuantiles':2e-7}}
write('protocol.json',protocol)
cohort=json.loads((study/'prediction-inputs/cohort.json').read_bytes());donors=cohort['eligibleDonors'];assert len(donors)==62
panel=json.loads((study/'prediction-inputs/panel.json').read_bytes())
features=sorted(panel,key=lambda x:hashlib.sha256(x.encode()).hexdigest())[:16]
groups=json.loads((study/'handoff/reference-groups.json').read_bytes());lookup={(g['expID'],g['conditionCode']):i for i,g in enumerate(groups)}
source=ad.read_h5ad(study/'handoff/B-source-codes-RNA.h5ad',backed='r')
reference=sparse.load_npz(study/'handoff/reference-pseudobulk.npz')
order=[source.var_names.get_loc(g.removeprefix('symbol|')) for g in features]
cases=[];records=[]
for donor in donors:
 gi=lookup[donor,'C'];g=groups[gi];rows=g['sourceCellIndices']
 assert all(source.obs.iloc[rows]['cond'].astype(str)=='C')
 assert all(source.obs.iloc[rows]['exp_id'].astype(str)==donor)
 X=source.X[rows,:].astype(np.int64).tocsr()
 np.testing.assert_array_equal(np.asarray(X.sum(axis=0)).ravel(),reference[gi].toarray().ravel())
 depths=np.asarray(X.sum(axis=1)).ravel();assert np.all(depths>0)
 for j,feature in zip(order,features):
  counts=X[:,j].toarray().ravel()
  for phi in protocol['dispersionSensitivity']:
   cases.append(dict(id=f'{donor}:{feature}:phi={phi}',counts=counts.tolist(),libraryCounts=depths.tolist(),cellDispersion=phi,gammaPriorShape=0.5,gammaPriorRatePerCPM=0.001,plannedLibraryCounts=[1000]*20))
 records.append(dict(donor=donor,cells=len(rows),libraryCounts=int(depths.sum()),sourceCellIndices=rows))
source.file.close()
write('cases.json',cases);write('selection.json',{'features':features,'donors':records,'cases':len(cases),'sourceGroupsChecked':len(records),'sourceCells':sum(r['cells'] for r in records),'zeroCountCases':sum(sum(c['counts'])==0 for c in cases)})
# Distinct numerical boundary cases supplement, never replace, experimental counts.
boundary=[]
for counts,depths in [([0],[1]),([0,0],[1000,10000]),([5,10],[1000,2000]),([1,500,0],[10,10000,100]),([100000000],[1000000000])]:
 for phi in [0,1e-8,0.001,1,100]:
  for a,b in [(0.1,1e-12),(0.5,0.001),(10,10)]:
   boundary.append(dict(id=f'boundary-{len(boundary)}',counts=counts,libraryCounts=depths,cellDispersion=phi,gammaPriorShape=a,gammaPriorRatePerCPM=b,plannedLibraryCounts=[1000,2000]))
write('boundary-cases.json',boundary)
write('freeze.json',{'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()},'prepareSHA256':sha(__file__)})
print(json.dumps(json.loads((out/'selection.json').read_bytes())|{'donors':'retained in selection.json'}))
