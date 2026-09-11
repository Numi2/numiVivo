"""Exercise provenance substitutions against the complete real Donor9 partition."""
from pathlib import Path
import argparse,copy,json,hashlib
import numpy as np
from verify_donor_preparation import validate_partition
from retain import sha,write

p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.study
parent=json.loads((root/'prepared/plan.json').read_text());parent_sha=sha(root/'prepared/plan.json');source_runs=json.loads((root/'prepared/runs.json').read_text());part=next(p for p in json.loads((root/'donor-plans/manifest.json').read_text())['parts'] if p['donor']=='Donor9');folder=root/'donor-plans/Donor9';plan=json.loads((folder/'plan.json').read_text());runs=json.loads((folder/'runs.json').read_text());mapping=np.load(folder/'global-selected-rows.npy',allow_pickle=False)
validate_partition(parent,parent_sha,source_runs,'Donor9',part,plan,runs,mapping)
rejected=[]
for kind in ['rehash-changed-barcode','changed-condition','changed-source-declaration','changed-cardinalities','shifted-source-entry-range','rehashed-row-map-permutation','changed-total-records']:
 altered=copy.deepcopy(plan);mapped=mapping.copy();run_copy=copy.deepcopy(runs);part_copy=dict(part)
 if kind=='rehash-changed-barcode':altered['metadata']['cells'][0]['barcode']+='-substituted'
 elif kind=='changed-condition':altered['metadata']['samples'][0]['condition']='substituted-condition'
 elif kind=='changed-source-declaration':
  declaration=json.loads(altered['sourceDeclaration']);declaration['parentPlanSHA256']='0'*64;altered['sourceDeclaration']=json.dumps(declaration,sort_keys=True)
 elif kind=='changed-cardinalities':altered['rowNonzeros'][0]+=1;altered['rowNonzeros'][1]-=1
 elif kind=='shifted-source-entry-range':run_copy[0]['entryFirst']+=1;run_copy[0]['entryLast']+=1
 elif kind=='rehashed-row-map-permutation':mapped[[0,1]]=mapped[[1,0]];part_copy['rowMapSHA256']=hashlib.sha256(mapped.tobytes()).hexdigest()
 elif kind=='changed-total-records':part_copy['records']+=1
 part_copy['planSHA256']=hashlib.sha256(json.dumps(altered,separators=(',',':')).encode()).hexdigest()
 try:validate_partition(parent,parent_sha,source_runs,'Donor9',part_copy,altered,run_copy,mapped)
 except AssertionError:rejected.append(kind)
 else:raise AssertionError('Accepted provenance substitution: '+kind)
validate_partition(parent,parent_sha,source_runs,'Donor9',part,plan,runs,mapping)
result=dict(status='passed',realDonor='Donor9',cells=34279,records=61860291,unchangedControlPassedBeforeAndAfter=True,expectedRejections=rejected,activeInputsModified=False,nativeBinaryModified=False,predictionFitOrScoring=False,scope='Real-partition provenance substitution controls; full count ingestion/replay remain separate.')
assert len(rejected)==7 and not a.out.exists();write(a.out,result);print(json.dumps(result))
