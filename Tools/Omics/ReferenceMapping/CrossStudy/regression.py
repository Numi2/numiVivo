#!/usr/bin/env python3
"""Native feature-panel lifecycle checks with explicit nonbiological count fixtures."""
import argparse,copy,json,os,subprocess,sys
from pathlib import Path
import anndata as ad,numpy as np,pandas as pd
from scipy.sparse import csr_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'Logistic'))
from prepare import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--old-binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);commands=[]
 def run(label,args,code=0,binary=None):
  r=subprocess.run([str(binary or a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr);commands.append(dict(label=label,args=list(map(str,args)),exitCode=r.returncode,expectedExitCode=code));write(a.out/'commands.json',commands);assert r.returncode==code,(label,r.stderr)
 rng=np.random.default_rng(731);genes=np.array(['g'+str(i) for i in range(32)]);labels=np.array(['a','b','c']*20)
 counts=np.array([rng.poisson(20+np.arange(32)*2+((np.arange(32)%3)==i%3)*100) for i in range(60)],dtype=np.int32)
 train=ad.AnnData(csr_matrix(counts),obs=pd.DataFrame({'sample':['train']*60,'label':labels},index=['t'+str(i) for i in range(60)]),var=pd.DataFrame(index=genes));train.write_h5ad(a.out/'train.h5ad')
 qcounts=counts[:20].copy();query=ad.AnnData(csr_matrix(qcounts),obs=pd.DataFrame({'sample':['query']*20},index=['q'+str(i) for i in range(20)]),var=pd.DataFrame(index=genes));query.write_h5ad(a.out/'query-same.h5ad')
 def mapping(sample,group=None):
  d=dict(schemaVersion=1,id='synthetic-panel-'+sample,evidence='synthetic',sourceDescription='Controlled feature-panel lifecycle fixture; not biological observations.',countUnit='umiCount',matrixPath='X',sampleColumn='sample',samples=[dict(id=sample,biologicalReplicateID=sample,batchID='fixture-batch',condition='fixture',organism='synthetic-organism')],mitochondrialFeatureIDs=[])
  if group:d['groupColumn']=group
  return d
 default=dict(schemaVersion=1,mapping=mapping('train','label'),featureNamespace='synthetic-genes',labelProvenance='Synthetic class definitions',neighbors=3,reduction=dict(pca=dict(components=2,highlyVariableFeatures=8,maximumBasis=32,meanBins=3)))
 fit=copy.deepcopy(default);fit['reduction']['pca']['featurePanel']=genes[:16].tolist();fit['logistic']=dict(penalty=1,gradientTolerance=1e-8,maximumIterations=2000)
 qp=dict(schemaVersion=1,mapping=mapping('query'),featureNamespace='synthetic-genes');write(a.out/'fit.json',fit);write(a.out/'query.json',qp);write(a.out/'default.json',default)
 run('fit',['singlecell-reference-fit',a.out/'train.h5ad','--plan',a.out/'fit.json','--output',a.out/'reference']);run('fit-verify',['singlecell-reference-verify',a.out/'reference'])
 model=json.loads((a.out/'reference/model.json').read_text());selected=np.array(model['reduction']['selectedFeatureIndices']);unselected=next(g for i,g in enumerate(genes[:16]) if i not in selected)
 # Replace a nonpanel feature with a different gene and different counts. Keep all panel genes.
 varied=query.copy();varied.var_names=list(genes[:-1])+['query-extra'];matrix=varied.X.tolil();matrix[:,-1]=5000;varied.X=matrix.tocsr();varied.write_h5ad(a.out/'query-varied.h5ad')
 def map_run(label,file,plan='query.json',code=0,reference='reference'):
  run(label,['singlecell-reference-map',a.out/file,'--plan',a.out/plan,'--reference',a.out/reference,'--output',a.out/label],code)
  if code:assert not (a.out/label).exists()
 map_run('mapped','query-varied.h5ad');run('map-verify',['singlecell-reference-map-verify',a.out/'mapped'])
 report=json.loads((a.out/'mapped/report.json').read_text());x=varied.X.astype(float).tocsr();totals=np.asarray(x.sum(axis=1)).ravel();x=x.multiply((10000/totals)[:,None]).tocsr();x.data=np.log1p(x.data)
 direct=x[:,selected]@np.array(model['reduction']['loadings'])-np.array(model['reduction']['projectionCenters'])@np.array(model['reduction']['loadings'])
 np.testing.assert_allclose(direct,[c['scores'] for c in report['cells']],rtol=1e-11,atol=1e-11);np.testing.assert_array_equal(totals,[c['totalCounts'] for c in report['cells']])
 varied[:,np.arange(31,-1,-1)].copy().write_h5ad(a.out/'query-reordered.h5ad');map_run('reordered','query-reordered.h5ad')
 changed=json.loads((a.out/'reordered/report.json').read_text());np.testing.assert_allclose([c['classProbabilities'] for c in changed['cells']],[c['classProbabilities'] for c in report['cells']],rtol=1e-11,atol=1e-11)
 varied[:,varied.var_names!=unselected].copy().write_h5ad(a.out/'query-missing-panel.h5ad');map_run('missing-unselected-panel-rejected','query-missing-panel.h5ad',code=65)
 varied[:,varied.var_names!='query-extra'].copy().write_h5ad(a.out/'query-without-extra.h5ad');map_run('outside-panel-absence-accepted','query-without-extra.h5ad')
 empty=varied.copy();matrix=empty.X.tolil();matrix[0,:]=0;empty.X=matrix.tocsr();empty.X.eliminate_zeros();empty.write_h5ad(a.out/'query-empty.h5ad');map_run('empty','query-empty.h5ad');er=json.loads((a.out/'empty/report.json').read_text());assert er['cells'][0].get('candidateLabel') is None and er['cells'][0].get('classProbabilities') is None and er['classifierOperations']==19*3*3
 for label,change in [('label-leak',lambda p:p['mapping'].update(groupColumn='label')),('namespace',lambda p:p.update(featureNamespace='different')),('budget',lambda p:p.update(maximumClassifierOperations=1))]:
  plan=copy.deepcopy(qp);change(plan);write(a.out/(label+'.json'),plan);map_run(label+'-rejected','query-varied.h5ad',plan=label+'.json',code=65)
 # Same positive source/plan inputs through old and new native binaries: compare public result bytes.
 for prefix,binary in [('old',a.old_binary),('new',a.binary)]:
  run(prefix+'-default-fit',['singlecell-reference-fit',a.out/'train.h5ad','--plan',a.out/'default.json','--output',a.out/(prefix+'-reference')],binary=binary)
  run(prefix+'-default-map',['singlecell-reference-map',a.out/'query-same.h5ad','--plan',a.out/'query.json','--reference',a.out/(prefix+'-reference'),'--output',a.out/(prefix+'-mapped')],binary=binary)
 assert (a.out/'old-reference/model.json').read_bytes()==(a.out/'new-reference/model.json').read_bytes()
 assert (a.out/'old-mapped/report.json').read_bytes()==(a.out/'new-mapped/report.json').read_bytes()
 map_run('default-universe-rejected','query-varied.h5ad',code=65,reference='new-reference')
 # Oversized fit plans reject at the same two-MiB envelope used by stored plans.
 (a.out/'oversized-fit.json').write_text(' '*2097153)
 run('oversized-plan-rejected',['singlecell-reference-fit',a.out/'train.h5ad','--plan',a.out/'oversized-fit.json','--output',a.out/'oversized'],65);assert not (a.out/'oversized').exists()
 write(a.out/'checks.json',dict(status='passed',commands=len(commands),binarySHA256=sha(a.binary),oldBinarySHA256=sha(a.old_binary),fullLibraryProjectionExact=True,missingUnselectedPanelGeneRejected=True,nonpanelFeaturesMayDiffer=True,emptyLibraryUnmapped=True,featureOrderInvariant=True,defaultModelAndReportBytesExact=True,queryLabelsAndNamespaceAndBudgetRejected=True,oversizedPlanRejected=True,checkerSHA256=sha(__file__)))
 print(json.dumps(json.loads((a.out/'checks.json').read_text())),flush=True)
if __name__=='__main__':main()
