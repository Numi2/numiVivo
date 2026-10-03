#!/usr/bin/env python3
"""Frozen official Nicheformer on admitted reference cells, never test responses.

Uses source-pinned external model code, MGI one-to-one orthologs and NCBI
Ensembl identifiers. Stereo-seq has no pretrained technology token: this is
an explicitly out-of-domain representation experiment, not validated transfer.
"""
import argparse,sys,csv,gzip,json
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from sklearn.decomposition import PCA
from wetlab import read,write,sha,require

HF_REV='0ba4ba162f9ad3d6c8b713761954c31be2fa69ed'
def extract(root,research):
 root=Path(root);research=Path(research);out=root/'frozen-Nicheformer';out.mkdir(exist_ok=False)
 features=read(root/'features.json');vocab=read(research/'nicheformer/vocab.json')
 human={}
 with gzip.open(research/'Homo_sapiens.gene_info.gz','rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   ids=[x.split(':',1)[1] for x in row['dbXrefs'].split('|') if x.startswith('Ensembl:') and x.split(':',1)[1] in vocab]
   if len(ids)==1:human[row['GeneID']]=ids[0]
 groups={}
 with open(research/'MGI-HOM_MouseHumanSequence.rpt') as f:
  for row in csv.DictReader(f,delimiter='\t'):groups.setdefault(row['DB Class Key'],[]).append(row)
 mapping={}
 for group in groups.values():
  mouse=[r for r in group if r['NCBI Taxon ID']=='10090'];man=[r for r in group if r['NCBI Taxon ID']=='9606']
  if len(mouse)==len(man)==1 and man[0]['EntrezGene ID'] in human:mapping[mouse[0]['Symbol']]=human[man[0]['EntrezGene ID']]
 # Reject duplicate mapped Ensembl columns rather than double-counting a gene.
 mapped=[(i,mapping[g]) for i,g in enumerate(features) if g in mapping]
 frequencies={g:sum(x[1]==g for x in mapped) for g in set(x[1] for x in mapped)}
 mapped=[(i,g) for i,g in mapped if frequencies[g]==1];require(len(mapped)>5000,'Insufficient mapped reference gene coverage')
 write(out/'representation-admission.json',{'model':'theislab/Nicheformer','revision':HF_REV,'weightsSHA256':sha(research/'nicheformer/model.safetensors'),'license':'MIT','species':'mouse via one-to-one MGI mouse-human homology plus NCBI human Ensembl IDs','sourceFiles':{n:sha(research/n) for n in ['MGI-HOM_MouseHumanSequence.rpt','Homo_sapiens.gene_info.gz']},'mappedGenes':len(mapped),'assayStatus':'OUT OF DISTRIBUTION: Stereo-seq unavailable, unknown technology PAD token, no technology mean normalization','representation':'official last-layer get_embeddings(with_context=False), official padding-inclusive pooling retained','input':'admitted control-associated reference rows only','training':False})
 write(out/'gene-mapping.json',[{'source':features[i],'humanOrtholog':g,'token':vocab[g]} for i,g in mapped])
 # Code was fetched at HF_REV and inspected before this local import. No remote
 # code execution or automated tokenizer axis expansion is used.
 package=research/'nicheformer';(package/'__init__.py').touch();sys.path.insert(0,str(research))
 from nicheformer.configuration_nicheformer import NicheformerConfig
 from nicheformer.modeling_nicheformer import NicheformerForMaskedLM
 import torch
 torch.set_num_threads(4);device='cpu'  # bounded offline embedding; MPS masked attention stalled on this Mac
 torch.backends.mha.set_fastpath_enabled(False)
 model=NicheformerForMaskedLM.from_pretrained(str(package),local_files_only=True).to(device).eval()
 embedding={};columns=np.array([i for i,g in mapped]);tokens=np.array([vocab[g] for i,g in mapped])
 for chip in ('chip2','chip3','chip1'):
  raw=load_file(str(root/(chip+'-reference.safetensors')));counts=raw['counts'][:,columns];seq=np.zeros((len(counts),1500),np.int64);seq[:,0]=4;seq[:,1]=6
  for i,row in enumerate(counts):
   nz=np.flatnonzero(row>0);order=nz[np.argsort(-row[nz],kind='stable')[:1497]];seq[i,3:3+len(order)]=tokens[order]
  chunks=[]
  with torch.inference_mode():
   for start in range(0,len(seq),4):
    batch=torch.from_numpy(seq[start:start+4]).to(device);chunks.append(model.nicheformer.get_embeddings(batch,batch!=0,layer=-1,with_context=False).cpu().numpy())
    if start%40==0: print(chip,start,'/',len(seq),device,flush=True)
  embedding[chip]=np.concatenate(chunks);require(np.isfinite(embedding[chip]).all(),'Nonfinite frozen representation');np.save(out/(chip+'-reference-embeddings.npy'),embedding[chip]);print(chip,len(seq),'frozen reference embeddings',flush=True)
 pca=PCA(n_components=16,svd_solver='full').fit(embedding['chip2']);scale=np.maximum(np.sqrt(pca.explained_variance_),.1)
 np.savez(out/'frozen-pca.npz',mean=pca.mean_,components=pca.components_,scale=scale)
 for chip,emb in embedding.items():
  ref=load_file(str(root/(chip+'-reference.safetensors')));reduced=pca.transform(emb)/scale;local=(reduced[ref['localIndices']]*ref['localWeights'][:,:,None]).sum(1);nhood=local[ref['neighbors']].mean(1)
  arrays=load_file(str(root/'neighborhood'/(chip+'.safetensors')));rows=read(root/(chip+'-rows.json'));anchors=[r['anchor'] for r in rows];arrays['descriptor'][:,:16]=local[anchors];arrays['descriptor'][:,16:32]=nhood[anchors];save_file(arrays,str(out/(chip+'.safetensors')))
 write(out/'plan.json',read(root/'neighborhood/plan.json'))
 write(out/'representation-receipt.json',{'admissionSHA256':sha(out/'representation-admission.json'),'hardware':device,'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()}})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--inputs',required=True,type=Path);p.add_argument('--research',required=True,type=Path);a=p.parse_args();extract(a.inputs,a.research)
