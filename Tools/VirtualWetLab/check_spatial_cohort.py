#!/usr/bin/env python3
"""Source equality and leakage qualification on the retained real cohort."""
import argparse
from pathlib import Path
import numpy as np,h5py,anndata as ad
from safetensors.numpy import load_file
from wetlab import read,write,sha,require

def check(cohort,inputs,sources):
 receipt=read(cohort/'cohort.json');checks=[]
 animals=[]
 for source in receipt['sources']:
  chip=source['chip'];a=ad.read_h5ad(cohort/(chip+'.h5ad'));n=a.n_obs;neighbors=a.obsm['neighbor_indices']
  require(not (neighbors==np.arange(n)[:,None]).any(),'Self neighbor')
  require(all(len(set(row))==15 for row in neighbors),'Repeated neighbor')
  require(a.obs.barcode_negative_is_untreated.eq(False).all(),'Negative barcode relabelled untreated')
  require(a.obs.loc[~a.obs.reference_coverage,'projected_cell_type'].astype(str).eq('uncertain').all(),'Unsupported confident annotation')
  negative=a.obs.guide_assignment.astype(str).eq('barcode-negative').to_numpy();reference=a.obs.reference_admitted.to_numpy()
  noncontrol=np.flatnonzero((a.obs.guide_species_count.to_numpy()>0)&~a.obs.guide_assignment.astype(str).eq('sgrna_msafe').to_numpy())
  require(not reference[neighbors[noncontrol].ravel()].any(),'Exposed neighbor enters reference')
  file=next(sources.glob(source['accession']+'*.gef'))
  with h5py.File(file) as f:
   names=[g.decode() for g in f['cellBin/gene']['geneName']];cells=f['cellBin/cell'];exp=f['cellBin/cellExp'];lookup={g:i for i,g in enumerate(a.var_names)}
   for row in np.random.default_rng(314159).choice(n,64,replace=False):
    c=cells[row];e=exp[int(c['offset']):int(c['offset'])+int(c['geneCount'])];expected={lookup[names[int(x['geneID'])]]:int(x['count']) for x in e if names[int(x['geneID'])] in lookup};actual=a.X[row];require(dict(zip(actual.indices.tolist(),actual.data.tolist()))==expected,'RNA differs from source')
  queries=read(inputs/(chip+'-rows.json'))
  for q in queries:
   ids=q['referenceRows']+q['neighborhoodReferenceRows'];require(reference[ids].all(),'Nonreference covariate');require(not set(ids)&set(q['outcomeRows']),'Target leakage')
  arrays=load_file(str(inputs/'neighborhood'/(chip+'.safetensors')))
  if chip=='chip1':require(set(arrays)=={'context','descriptor','target','known'},'Test outcomes in prediction input')
  checks.append({'chip':chip,'cells':n,'allNeighborRowsChecked':n,'sourceRNARowsCompared':64,'queriesChecked':len(queries),'referenceCoverage':int(a.obs.reference_coverage.sum()),'passed':True});animals.append(set(source['animals']))
 require(all(not a&b for i,a in enumerate(animals) for b in animals[i+1:]),'Animal overlap')
 return {'format':'numivivo-spatial-cohort-qualification/v1','cohortSHA256':sha(cohort/'cohort.json'),'inputsSHA256':sha(inputs/'prepared.json'),'checks':checks,'animalOverlap':False,'softwareQualified':True,'biologicalValidation':False}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--cohort',type=Path,required=True);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();write(a.output,check(a.cohort,a.inputs,a.sources))
