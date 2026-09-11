"""Freeze disjoint complete donor partitions for restartable count execution."""
from pathlib import Path
import os,json,hashlib,numpy as np
ROOT=Path(os.environ['NUMIVIVO_PARSE_STUDY']);out=ROOT/'donor-plans';out.mkdir(exist_ok=False)
raw=(ROOT/'prepared/plan.json').read_bytes();parent=json.loads(raw);runs=json.loads((ROOT/'prepared/runs.json').read_text());metadata=parent['metadata'];samples={s['id']:s for s in metadata['samples']};donors=json.loads((ROOT/'sources/roster.json').read_text())['categories']['donor'];all_indices=[];all_runs=[];parts=[]
for donor in donors:
 chosen=[r for r in runs if samples[r['sample']]['donorID']==donor];indices=[];local_runs=[]
 for run in chosen:
  mapped=dict(run,globalLo=run['lo'],globalHi=run['hi'],lo=len(indices),hi=len(indices)+run['hi']-run['lo']);local_runs.append(mapped);indices+=list(range(run['lo'],run['hi']))
 cells=[metadata['cells'][i] for i in indices];selected_samples=[s for s in metadata['samples'] if s['donorID']==donor];assert len(selected_samples)==2 and {s['condition'] for s in selected_samples}=={'PBS','IFN-beta'}
 assert all(samples[c['sampleID']]['donorID']==donor for c in cells)
 declaration=json.loads(parent['sourceDeclaration']);declaration.update(partitionDonor=donor,parentPlanSHA256=hashlib.sha256(raw).hexdigest(),globalSelectedRowIndicesSHA256=hashlib.sha256(np.array(indices,dtype='<u8').tobytes()).hexdigest())
 plan=dict(parent,metadata=dict(metadata,id=metadata['id']+'-'+donor,samples=selected_samples,cells=cells),sourceDeclaration=json.dumps(declaration,sort_keys=True),rowNonzeros=[parent['rowNonzeros'][i] for i in indices])
 if 'rowTotals' in parent:plan['rowTotals']=[parent['rowTotals'][i] for i in indices]
 directory=out/donor;directory.mkdir();data=json.dumps(plan,separators=(',',':')).encode();(directory/'plan.json').write_bytes(data);(directory/'runs.json').write_text(json.dumps(local_runs,indent=2)+'\n');np.save(directory/'global-selected-rows.npy',np.array(indices,dtype='<u8'))
 parts.append(dict(donor=donor,cells=len(indices),records=sum(plan['rowNonzeros']),runs=len(local_runs),planSHA256=hashlib.sha256(data).hexdigest(),rowMapSHA256=declaration['globalSelectedRowIndicesSHA256']));all_indices+=indices;all_runs += [r['run'] for r in chosen]
assert len(parts)==12 and sorted(all_indices)==list(range(725031)) and sorted(all_runs)==list(range(3456)) and sum(p['records'] for p in parts)==1373870697
(out/'manifest.json').write_text(json.dumps(dict(status='frozen-disjoint-complete-donor-partitions',parentPlanSHA256=hashlib.sha256(raw).hexdigest(),cells=725031,features=40352,records=1373870697,parts=parts,scope='Execution partition only; same complete source cohort and feature axis. Each donor stream has local row indices and its own SHA256; no monolithic stream hash is claimed.'),indent=2)+'\n');print(json.dumps(parts))
