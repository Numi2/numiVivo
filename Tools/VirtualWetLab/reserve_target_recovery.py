#!/usr/bin/env python3
"""One reserved experimental preparation after the target-learning diagnosis.

Uses the existing native learner and selection scorer. No expression from
assigned treated cells is read before the prediction seal.
"""
import argparse,json
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp,inventory
from prepare_intervention_design import fetch_descriptors,normalized
from train_intervention_design import native,VARIANTS
from design_objective import selection_score

def prepare(root,inputs):
 out=root/'prepared';require(not out.exists(),'Retain prior preparation');out.mkdir()
 protocol=read(root/'preregistration-v2.json');features=read(inputs/'features.json');a=ad.read_h5ad(root/'GSM2406677.h5ad',backed='r');original=pd.read_csv(root/'original-cell-identities.csv.gz')
 lookup={}
 for _,r in original.iterrows():
  key=(r['cell BC'].split('-')[0],r['guide identity'],int(r['read count']),int(r['UMI count']));lookup.setdefault(key,[]).append(str(r['cell BC']))
 gem=[];unresolved=[]
 for i,(barcode,r) in enumerate(a.obs.iterrows()):
  matches=lookup.get((barcode,str(r.perturbation),int(r['read count']),int(r['UMI count'])),[]) if pd.notna(r['read count']) and pd.notna(r['UMI count']) else []
  gem.append(matches[0].rsplit('-',1)[1] if len(matches)==1 else None)
  if len(matches)!=1:unresolved.append({'row':i,'matches':len(matches),'reason':'original barcode/guide/read/UMI identity not unique'})
 labels=a.obs.perturbation.astype(str).to_numpy();gem=np.asarray(gem);aliases={'ATF6_only_pMJ145':'ATF6','PERK_only_pMJ146':'EIF2AK3','IRE1_only_pMJ148':'ERN1'}
 desc=fetch_descriptors(sorted(aliases.values()),root/'target-sources');cols=[a.var_names.get_loc(g) for g in features];obj=features.index('HSPA5');rows=[];arrays={k:[] for k in ('context','descriptor','target','known')};excluded=[];coverage=[]
 for block in sorted(x for x in set(gem) if x is not None):
  controls=np.flatnonzero((gem==block)&np.char.startswith(labels.astype(str),'3x_neg_ctrl_pMJ144'))
  c=a.X[controls].tocsr();base=np.log1p(c[:,cols].toarray()/np.maximum(1,np.asarray(c.sum(1)).ravel())[:,None]*1e6).mean(0).astype(np.float32)
  detected=float((c[:,cols[obj]].toarray()>0).mean());coverage.append({'capture':block,'controls':len(controls),'HSPA5DetectedFraction':detected})
  for label,target in aliases.items():
   idx=np.flatnonzero((gem==block)&(labels==label));reasons=[]
   if len(idx)<100 or len(controls)<100:reasons.append('fewer than 100 treated or explicit control cells')
   if detected<.2:reasons.append('objective control detection below 20%')
   if target not in desc:reasons.append('source-bound target representation unavailable')
   if reasons:excluded.append({'capture':block,'target':target,'reasons':reasons});continue
   rows.append({'target':target,'label':label,'rows':idx.tolist(),'controls':controls.tolist(),'capture':block,'unit':'GSM2406677 pooled epistasis preparation','biologicalUnits':1,'captureMeaning':'original GEM group; not independent biological replicate','condition':{'1':'tunicamycin','2':'thapsigargin','3':'DMSO'}[block],'source':'GSE90546/GSM2406677','modality':'CRISPRi','context':'K562'})
   arrays['context'].append(base);arrays['descriptor'].append(desc[target]['descriptor']+[0.,1.,0.,0.,1.,0.,0.]);arrays['target'].append(0);arrays['known'].append([0.])
 write(out/'coverage.json',coverage);write(out/'exclusions.json',excluded);write(out/'identity-unresolved.json',unresolved);write(out/'rows.json',rows);write(out/'targets.json',desc);write(out/'features.json',features)
 require(len(rows)>=2,'Fewer than two supported reserved candidates');save_file({k:np.asarray(v,np.int32 if k=='target' else np.float32) for k,v in arrays.items()},str(out/'query.safetensors'))
 write(out/'registration.json',{'createdAt':timestamp(),'protocolSHA256':sha(root/'preregistration-v2.json'),'sourceSHA256':sha(root/'GSM2406677.h5ad'),'originalIdentitiesSHA256':sha(root/'original-cell-identities.csv.gz'),'sourceMetadataSHA256':sha(root/'GSM2406677-source.txt'),'outcomesOpened':False,'independence':'one new experimental preparation within an exposed study; three original captures cannot yield biological confidence','descriptorAliasSource':'https://github.com/thomasmaxwellnorman/perturbseq_demo/blob/master/perturbseq_demo.ipynb','authorNotebookSHA256':sha(root/'author-notebook.json'),'objective':'UPR-HSPA5-v2','files':inventory(out)})
 a.file.close()

def predict(root,campaign,binary):
 out=root/'predictions';require(not out.exists(),'Retain prediction attempt');out.mkdir();q=load_file(str(root/'prepared/query.safetensors'));reg=read(campaign/'registration.json');features=read(root/'prepared/features.json');results={};training=load_file(str(Path(reg['inputs'])/'training.safetensors'));tr=read(Path(reg['inputs'])/'training-rows.json');delta=training['observed']-training['context'];mask=np.asarray([r['modality']=='CRISPRi' for r in tr]);results['no-change']=q['context'];results['training-mean']=q['context']+delta[mask].mean(0)
 for variant in VARIANTS:
  folder=campaign/variant;query={k:v.copy() for k,v in q.items()}
  if variant=='target-ID':query['descriptor'][:,:149]=0
  if variant=='shuffled-target-descriptor':query['descriptor'][:,:149]=np.roll(query['descriptor'][:,:149],1,axis=0)
  save_file(query,str(out/(variant+'.safetensors')));step=read(folder/'selection.json')['selected']['step'];native(binary,'predict',folder/'plan.json',out/(variant+'.safetensors'),out/variant,folder/'training'/f'weights-{step}.safetensors');results[variant]=load_file(str(out/variant/'prediction.safetensors'))['mean']
 save_file(results,str(out/'all-predictions.safetensors'));write(out/'seal.json',{'createdAt':timestamp(),'objective':'UPR-HSPA5-v2','preparedRegistrationSHA256':sha(root/'prepared/registration.json'),'sourceCampaignSealSHA256':sha(campaign/'prediction-seal.json'),'binarySHA256':sha(binary),'files':inventory(out),'outcomesOpened':False})

def reveal(root):
 out=root/'evaluation';require(not out.exists(),'Already evaluated');out.mkdir();s=read(root/'predictions/seal.json')
 for p,h in s['files'].items():require(sha(root/'predictions'/p)==h,'Changed sealed prediction')
 reg=read(root/'prepared/registration.json');require(sha(root/'GSM2406677.h5ad')==reg['sourceSHA256'],'Changed source');require(sha(root/'prepared/registration.json')==s['preparedRegistrationSHA256'],'Changed registration')
 for p,h in reg['files'].items():require(sha(root/'prepared'/p)==h,'Changed prepared input')
 a=ad.read_h5ad(root/'GSM2406677.h5ad');rows=read(root/'prepared/rows.json');features=read(root/'prepared/features.json');cols=[a.var_names.get_loc(g) for g in features];q=load_file(str(root/'prepared/query.safetensors'));pred=load_file(str(root/'predictions/all-predictions.safetensors'));y=np.asarray([normalized(a,r['rows'],cols).mean(0) for r in rows]);save_file({'mean':y},str(out/'observations.safetensors'));j=features.index('HSPA5');scores={}
 for block in sorted({r['capture'] for r in rows}):
  ix=[i for i,r in enumerate(rows) if r['capture']==block];candidates=[rows[i]['target'] for i in ix]+['no-intervention'];observed=(-(y[ix,j]-q['context'][ix,j])).tolist()+[0.];baseline=(-(pred['training-mean'][ix,j]-q['context'][ix,j])).tolist()+[0.]
  scores[block]={'candidates':candidates,'observedUtility':observed,'models':{}}
  for model,p in pred.items():
   utility=(-(p[ix,j]-q['context'][ix,j])).tolist()+[0.];scores[block]['models'][model]={'selection':selection_score(candidates,utility,observed,baseline),'predictedUtility':utility,'endpointRMSE':float(np.sqrt(np.mean((p[ix]-y[ix])**2))),'objectiveResidual':(p[ix,j]-y[ix,j]).tolist()}
 write(out/'result.json',{'openedAt':timestamp(),'objective':'UPR-HSPA5-v2','scores':scores,'biologicalPromotion':False,'scope':'One reserved experimental preparation; captures are not biological replicates. GEM conditions resolved by author notebook: 1 tunicamycin, 2 thapsigargin, 3 DMSO; query uses condition-matched control RNA, no dose claim. No spatial validation.','predictionSealSHA256':sha(root/'predictions/seal.json')})
 write(out/'seal.json',{'files':inventory(out)})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','predict','reveal']);p.add_argument('--root',type=Path,required=True);p.add_argument('--inputs',type=Path);p.add_argument('--campaign',type=Path);p.add_argument('--binary',type=Path);a=p.parse_args();prepare(a.root,a.inputs) if a.mode=='prepare' else predict(a.root,a.campaign,a.binary) if a.mode=='predict' else reveal(a.root)
