"""Bind every frozen donor axis and source range map to immutable preparation."""
import argparse,hashlib,json,tarfile,re
from pathlib import Path
import numpy as np
from retain import sha,write

PREPARATION_SHA256='f08c78ca7a74ab01eab6c3f7748615f45e217cbf95446ba35ef830eaf9c75d02'


def validate_partition(parent,parent_sha,source_runs,donor,part,plan,runs,mapping):
 samples={s['id']:s for s in parent['metadata']['samples']}
 chosen=[r for r in source_runs if samples[r['sample']]['donorID']==donor]
 expected_rows=[];expected_runs=[]
 for run in chosen:
  expected_runs.append(dict(run,globalLo=run['lo'],globalHi=run['hi'],lo=len(expected_rows),hi=len(expected_rows)+run['hi']-run['lo']))
  expected_rows.extend(range(run['lo'],run['hi']))
 assert runs==expected_runs,'Donor source ranges differ from immutable preparation'
 expected_map=np.array(expected_rows,dtype='<u8')
 assert mapping.dtype==np.dtype('<u8') and np.array_equal(mapping,expected_map),'Donor row map differs from immutable source selection'
 map_sha=hashlib.sha256(expected_map.tobytes()).hexdigest()
 assert part['rowMapSHA256']==map_sha and part['donor']==donor and part['cells']==len(expected_rows) and part['runs']==len(chosen)
 declaration=json.loads(parent['sourceDeclaration']);declaration.update(partitionDonor=donor,parentPlanSHA256=parent_sha,globalSelectedRowIndicesSHA256=map_sha)
 metadata=dict(parent['metadata'],id=parent['metadata']['id']+'-'+donor,samples=[s for s in parent['metadata']['samples'] if s['donorID']==donor],cells=[parent['metadata']['cells'][i] for i in expected_rows])
 expected=dict(parent,metadata=metadata,sourceDeclaration=json.dumps(declaration,sort_keys=True),rowNonzeros=[parent['rowNonzeros'][i] for i in expected_rows])
 if 'rowTotals' in parent:expected['rowTotals']=[parent['rowTotals'][i] for i in expected_rows]
 assert plan==expected,'Donor plan differs from exact immutable parent projection'
 assert part['records']==sum(expected['rowNonzeros'])
 return expected_rows,[r['run'] for r in chosen]


def verify(study,archive):
 assert sha(archive)==PREPARATION_SHA256,'Preparation archive identity changed'
 with tarfile.open(archive,'r:gz') as t:members=json.load(t.extractfile('members.json'))
 # These are all source inputs used by the count adapter and donor partitioner.
 inputs={name:v for name,v in members.items() if name.startswith(('prepared/','axes/','chunks/')) or name in ['sources/asset-head.json','sources/roster.json']}
 assert {n for n in inputs if re.fullmatch(r'axes/run-\d+\.npz',n)}=={f'axes/run-{i:04d}.npz' for i in range(3456)}
 assert {n for n in inputs if re.fullmatch(r'chunks/run-\d+\.json',n)}=={f'chunks/run-{i:04d}.json' for i in range(3456)}
 for name,entry in inputs.items():
  path=study/name;assert path.is_file() and not path.is_symlink() and path.stat().st_size==entry['bytes'] and sha(path)==entry['SHA256'],name
 parent=json.loads((study/'prepared/plan.json').read_text());parent_sha=sha(study/'prepared/plan.json');source_runs=json.loads((study/'prepared/runs.json').read_text());frozen=json.loads((study/'donor-plans/manifest.json').read_text());donors=json.loads((study/'sources/roster.json').read_text())['categories']['donor']
 assert frozen['status']=='frozen-disjoint-complete-donor-partitions' and frozen['parentPlanSHA256']==parent_sha and len(frozen['parts'])==len(donors)==12
 assert [p['donor'] for p in frozen['parts']]==donors and frozen['cells']==725031 and frozen['features']==40352 and frozen['records']==1373870697
 rows=[];runs=[];part_files={'donor-plans/manifest.json':sha(study/'donor-plans/manifest.json')}
 for donor,part in zip(donors,frozen['parts']):
  directory=study/'donor-plans'/donor;path=directory/'plan.json';assert sha(path)==part['planSHA256']
  chosen_rows,chosen_runs=validate_partition(parent,parent_sha,source_runs,donor,part,json.loads(path.read_text()),json.loads((directory/'runs.json').read_text()),np.load(directory/'global-selected-rows.npy',allow_pickle=False))
  rows+=chosen_rows;runs+=chosen_runs
  for name in ['plan.json','runs.json','global-selected-rows.npy']:
   path=directory/name;assert path.is_file() and not path.is_symlink();part_files[str(path.relative_to(study))]=sha(path)
 assert sorted(rows)==list(range(725031)) and sorted(runs)==list(range(3456)) and sum(p['records'] for p in frozen['parts'])==1373870697
 return dict(status='passed-input-bindings',preparationArchiveSHA256=PREPARATION_SHA256,preparationInputs=len(inputs),donorFiles=len(part_files),donors=12,cells=725031,records=1373870697,allOriginalPreparationInputsExact=True,allDonorPlansExactParentProjections=True,allSourceRowsCoveredExactlyOnce=True,inputSHA256={**{n:v['SHA256'] for n,v in inputs.items()},**part_files},scope='Source identities, metadata, cardinalities and range maps only. Does not claim count ingestion, native replay or prediction completion.')


def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--preparation',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();result=verify(a.study,a.preparation);assert not a.out.exists();write(a.out,result);print(json.dumps({k:v for k,v in result.items() if k!='inputSHA256'}))
if __name__=='__main__':main()
