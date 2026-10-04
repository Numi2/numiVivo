#!/usr/bin/env python3
"""Frozen ESM2 sequence features; no perturbation expression or model fitting."""
import argparse,hashlib,json,time
from pathlib import Path
import numpy as np
from safetensors.numpy import save_file
from wetlab import read,write,sha,require,timestamp

def run(targets,sources,out,reuse=None):
 from huggingface_hub import HfApi,snapshot_download
 import torch
 from transformers import AutoTokenizer,EsmModel
 require(not out.exists(),'Retain previous embedding attempt');out.mkdir()
 modelID='facebook/esm2_t33_650M_UR50D';revision=read(reuse/'registration.json')['revision'] if reuse else HfApi().model_info(modelID,files_metadata=True).sha
 write(out/'registration.json',{'createdAt':timestamp(),'model':modelID,'revision':revision,'modelURL':'https://huggingface.co/'+modelID,'license':'MIT per official model card','pretraining':'UniRef50/D 2021_04 protein sequences; no perturbation-response expression training claimed','targetManifestSHA256':sha(targets),'ownerSHA256':sha(__file__),'representation':'last hidden layer33; residue mean excluding special tokens; contiguous <=1022-residue chunks, length-weighted mean; no fine tuning','reduction':'not here; PCA fitted to unique training targets inside each nested fold','device':'Apple MPS float32','torch':torch.__version__})
 path=reuse/'model' if reuse else Path(snapshot_download(modelID,revision=revision,allow_patterns=['model.safetensors','config.json','tokenizer_config.json','special_tokens_map.json','vocab.txt','README.md'],local_dir=out/'model'))
 if reuse:
  from safetensors.numpy import load_file
  cached=load_file(str(reuse/'embeddings.safetensors'));provenance=read(reuse/'provenance.json')
  for f,h in provenance['modelFiles'].items():require(sha(path/f)==h,'Frozen model changed')
 tokenizer=AutoTokenizer.from_pretrained(path,local_files_only=True);model=EsmModel.from_pretrained(path,local_files_only=True,add_pooling_layer=False).eval().to('mps');model.requires_grad_(False)
 vectors={};records={};start=time.monotonic()
 for target,record in sorted(read(targets).items()):
  p=sources/(record['accession']+'.json');require(sha(p)==record['sourceSHA256'],'Sequence source changed');a=read(p);seq=a['sequence']['sequence'];require(hashlib.sha256(seq.encode()).hexdigest()==record['sequenceSHA256'],'Sequence changed');require(a['organism']['taxonomy']==9606,'Unsupported species')
  if reuse and provenance['records'][target]['sequenceSHA256']==record['sequenceSHA256']:
   vectors[target]=cached[target];records[target]=provenance['records'][target];continue
  means=[];lengths=[]
  for i in range(0,len(seq),1022):
   chunk=seq[i:i+1022];tokens=tokenizer(chunk,return_tensors='pt').to('mps')
   with torch.inference_mode():v=model(**tokens).last_hidden_state[0,1:len(chunk)+1].mean(0).cpu().numpy()
   means.append(v);lengths.append(len(chunk))
  vectors[target]=np.average(means,axis=0,weights=lengths).astype(np.float32);records[target]={'accession':record['accession'],'sourceURL':record['sourceURL'],'sourceSHA256':sha(p),'sequenceSHA256':record['sequenceSHA256'],'species':'human','length':len(seq),'chunks':lengths};print(target,len(seq),'embedded',flush=True)
 save_file(vectors,str(out/'embeddings.safetensors'));write(out/'provenance.json',{'modelDirectory':str(path.resolve()),'reusedFrom':str(reuse) if reuse else None,'records':records,'seconds':time.monotonic()-start,'modelFiles':{p.name:sha(p) for p in path.iterdir() if p.is_file()},'embeddingSHA256':sha(out/'embeddings.safetensors'),'biologicalValidation':False})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--targets',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--reuse',type=Path);a=p.parse_args();run(a.targets,a.sources,a.output,a.reuse)
