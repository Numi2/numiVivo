#!/usr/bin/env python3
"""Execute the pinned Arc six-metric evaluator on public human CRISPRi counts.

This is a local Replogle-2020 no-change comparator, not VCC-2026 submission
conformance, a returned challenge score, or spatial-model qualification.
"""
import argparse,gzip,subprocess,re
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy.sparse import coo_matrix
from wetlab import read,write,sha,require,timestamp
SCORER='0ce25cf4c495e76e60dfca4ddc8f822507aa0e18'
MATRIX_SHA='0684ccfa61d8460ecd9765afe020441de152e4ae2b91b9a5c6000221a9b1ea1a'
def run(metadata,matrix,evaluator,out):
 out=Path(out);out.mkdir(exist_ok=False);metadata=Path(metadata);matrix=Path(matrix)
 require(sha(matrix)==MATRIX_SHA,'Original GEO counts hash mismatch')
 bars=pd.read_csv(metadata/'GSM4367979_exp1-5.barcodes.tsv.gz',header=None)[0].tolist();barindex={b:i for i,b in enumerate(bars)};ids=pd.read_csv(metadata/'GSM4367979_exp1-5.cell_identities.csv.gz');features=pd.read_csv(metadata/'GSM4367979_exp1-5.features.tsv.gz',sep='\t',header=None,names=['ensembl','symbol','featureType'])
 # Unique guide assignment, deposited QC, one experiment block. No response-
 # based target or gene selection. No gemgroup-as-animal interpretation.
 multiplicity=ids.groupby('cell_barcode').guide_identity.nunique();ids=ids[(ids.gemgroup==5)&ids.good_coverage&(ids.number_of_cells==1)&ids.cell_barcode.isin(barindex)&ids.cell_barcode.map(multiplicity).eq(1)].drop_duplicates('cell_barcode')
 support=ids.guide_identity.value_counts();targets=sorted(x for x,n in support.items() if n>=40 and x not in ('sgNegCtrl2','sgNegCtrl3') and re.fullmatch(r'sg[A-Za-z0-9]+',x))[:8];require(len(targets)>=4,'Too few supported local targets')
 controls=ids[ids.guide_identity.isin(['sgNegCtrl2','sgNegCtrl3'])].sort_values('cell_barcode').head(200);require(len(controls)>=40,'Too few deposited negative controls')
 selected=[*controls.cell_barcode];labels=['non-targeting']*len(selected)
 for t in targets:
  b=ids[ids.guide_identity==t].sort_values('cell_barcode').head(60).cell_barcode.tolist();selected+=b;labels += [t[2:]]*len(b)
 selectedindices={barindex[b]:i for i,b in enumerate(selected)}
 keep=(features.featureType=='Gene Expression')&features.ensembl.str.startswith('ENSG')&~features.symbol.duplicated(keep=False)&features.symbol.notna();geneindices={j:i for i,j in enumerate(np.flatnonzero(keep))};genes=features.loc[keep,'symbol'].tolist();require(len(genes)>10000,'RNA feature axis admission failed')
 write(out/'plan.json',{'source':'GSE146194/GSM4367979','sourceSHA256':MATRIX_SHA,'scorerRevision':SCORER,'profile':'vcc2026','targets':targets,'selection':'first eight lexicographic QC-supported guides, gemgroup5, at most60 measured cells each; unique RNA symbols','prediction':'deterministic resampling of measured negative controls; raw integral counts, no-change baseline','seed':314159,'geneCount':len(genes),'context':'K562 gemgroup5; not a new biological cell context','challengeConformance':'UNAVAILABLE: official 18533-gene panel/control bundle requires VCC login','returnedChallengeScore':None,'spatialQualification':False,'createdAt':timestamp()})
 rr=[];cc=[];vv=[]
 with gzip.open(matrix,'rt') as f:
  for line in f:
   if not line.startswith('%'):break
  for line in f:
   g,c,v=map(int,line.split());row=selectedindices.get(c-1);col=geneindices.get(g-1)
   if row is not None and col is not None:rr.append(row);cc.append(col);vv.append(v)
 counts=coo_matrix((np.asarray(vv,np.int32),(rr,cc)),shape=(len(selected),len(genes))).tocsr();real=ad.AnnData(counts,obs=pd.DataFrame({'target':labels,'context':'K562-gemgroup5'},index=selected),var=pd.DataFrame(index=genes))
 # Prediction written/sealed before real observation export and scorer execution.
 rng=np.random.default_rng(314159);sample=np.concatenate([rng.choice(len(controls),400,replace=True) for _ in targets]);pred=ad.AnnData(counts[sample],obs=pd.DataFrame({'target':np.repeat([x[2:] for x in targets],400),'context':'K562-gemgroup5'},index=['prediction-'+str(i) for i in range(len(sample))]),var=real.var.copy());pred.write_h5ad(out/'prediction.h5ad',compression='gzip');write(out/'prediction-seal.json',{'predictionSHA256':sha(out/'prediction.h5ad'),'planSHA256':sha(out/'plan.json'),'scorerRevision':SCORER,'model':'no-change measured-control resampling','newIndependentAnimals':0})
 real.write_h5ad(out/'observations.h5ad',compression='gzip')
 # Official local run API requires matching perturbation sets including the
 # real control label, unlike the challenge upload format which rejects controls.
 local=ad.concat([pred,real[real.obs.target=='non-targeting'].copy()],index_unique='-local');local.write_h5ad(out/'local-evaluator-input.h5ad',compression='gzip')
 command=[str(evaluator),'run','--preset','vcc2026','--adata-pred',str(out/'local-evaluator-input.h5ad'),'--adata-real',str(out/'observations.h5ad'),'--outdir',str(out/'official-evaluation'),'--set','num_threads=4']
 with (out/'official-evaluator.log').open('w') as log:result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
 write(out/'execution.json',{'command':command,'exitCode':result.returncode,'scorerRevision':SCORER,'sourceSHA256':MATRIX_SHA,'predictionSHA256':sha(out/'prediction.h5ad'),'observationsSHA256':sha(out/'observations.h5ad'),'localEvaluation':True,'challengeSubmissionConformance':False,'returnedChallengeScore':None,'biologicalContextTransferEstablished':False,'spatialQualification':False});require(result.returncode==0,'Official evaluator failed; log retained');print('Official six-metric local evaluation executed',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--metadata',type=Path,required=True);p.add_argument('--matrix',type=Path,required=True);p.add_argument('--evaluator',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.metadata,a.matrix,a.evaluator,a.output)
