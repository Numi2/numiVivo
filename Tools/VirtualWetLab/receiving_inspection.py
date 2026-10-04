"""Read-only spatial development inspection of explicit receiving populations."""
from pathlib import Path
from functools import lru_cache
import numpy as np
import anndata as ad
from safetensors.numpy import load_file
from wetlab import read,sha,require
from prepare_intervention_design import V03,orthology

@lru_cache(maxsize=2)
def artifacts(root,seal_hash):
    root=Path(root);s=read(root/'artifact-seal.json')
    for f,h in s['files'].items():require(sha(root/f)==h,'Changed receiver artifact')
    return read(root/'comparison.json'),load_file(str(root/'development-results.safetensors'))

def inspect(root,gene=None,group=None):
    root=Path(root);comparison,tensors=artifacts(str(root.resolve()),sha(root/'artifact-seal.json'));rows=comparison['groups']
    result={'scope':'EXPOSED SPATIAL DEVELOPMENT','genes':read(root/'scored-features.json'),'groups':[{'id':str(i),'target':r['target'],'role':r['role'],'receivingType':r['receivingType'],'cells':r['cells']} for i,r in enumerate(rows)],
      'verifiedTissueEdges':False,'biologicalPromotion':False,'spatialBenefit':'UNAVAILABLE: cell-to-section identity unresolved'}
    registration=read(root/'registration.json');result['modelVersion']=registration.get('modelVersion','v04-legacy-receiver');result['runtimeSHA256']=registration.get('binarySHA256',read(root/'artifact-seal.json').get('binarySHA256'));result['artifactSealSHA256']=sha(root/'artifact-seal.json')
    result['nativeModels']={p.parent.name:read(p)['selected'] for p in root.glob('*/selection.json')}
    if gene is None:return result
    require(group is not None and group.isdigit() and int(group)<len(rows),'Receiver population required');i=int(group);r=rows[i]
    from intervention_design import bound
    if (root/'features.json').exists():features=read(root/'features.json')
    else:
        _,inputs,_=bound(Path(read(root/'registration.json')['pretrainedWeights']).parent)
        features=read(inputs/'features.json')
    mapping=orthology(V03/'research/MGI-HOM_MouseHumanSequence.rpt');canonical=mapping.get(gene,gene)
    require(canonical in features and canonical in read(root/'scored-features.json'),'Gene unavailable on the common measured model axis');j=features.index(canonical)
    source=V03/'cohort-v3/chip1.h5ad';require(sha(source)==r['sourceSHA256'],'Measured geometry source changed');a=ad.read_h5ad(source,backed='r');index={mapping.get(str(g)):k for k,g in enumerate(a.var_names) if mapping.get(str(g))};k=index[canonical];cells=r['outcomeRows'];x=a.X[cells].tocsr();raw=x[:,k].toarray().ravel();values=np.log1p(raw/np.maximum(1,np.asarray(x.sum(1)).ravel())*1e6);a.file.close()
    prediction=load_file(str(root/'receiver-pretrained/development-prediction/prediction.safetensors'));control=float(tensors['control'][i,j]);pred=float(tensors['receiver-pretrained'][i,j]);obs=float(tensors['observed'][i,j]);variance=float(prediction['variance'][i,j])
    result.update(gene=gene,canonicalGene=canonical,population=r['receivingType'],populationRows=cells,sourceSHA256=r['sourceSHA256'],support=r,
      modelErrors={m:v['groupRMSE'][i] for m,v in comparison['models'].items()},
      field={'gene':gene,'units':'log1p(CPM)','rawUnits':'UMI counts','measuredInSource':True,'modeled':True,'referenceRows':[],'referenceRaw':[],'referenceNormalized':[],'populationEffects':[],
       'arms':[{'id':str(i),'target':r['target'],'role':'neighbor' if r['role']=='receiver' else 'direct','control':control,'predicted':pred,'observed':obs,'variance':variance,'residual':pred-obs,'observedRows':cells,'observedValues':values.tolist(),'observedRaw':raw.astype(float).tolist(),'evidence':'MODEL INFERENCE'}]})
    return result
