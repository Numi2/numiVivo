#!/usr/bin/env python3
"""Independent source, numerical and experiment-lifecycle qualification."""
import argparse,copy,json,shutil
from pathlib import Path
import numpy as np,anndata as ad,h5py
from wetlab import read,write,sha,require
from molecular_adapter import MolecularPerturbationAdapter,score
p=argparse.ArgumentParser();p.add_argument('--assay',type=Path,required=True);p.add_argument('--sources',type=Path,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
write(a.output/'plan.json',{'checks':['exact GEF projection','native aggregate oracle','independent ridge and all baselines','pre-reveal isolation','held-out RMSE','exact native replay','unsupported requests','tamper rejection'],'tolerance':1e-9,'stoppingRule':'one complete check; preserve failures','purpose':'software plus retrospective development evidence; no biological qualification'})
assay=read(a.assay);root=a.assay.parent;lineage=read(root/'lineage.json');x=ad.read_h5ad(root/'counts.h5ad');adapter=MolecularPerturbationAdapter();before=adapter.summary(a.run)
assert not before['revealed'] and 'comparison' not in before
assert all(not m['evidence']['heldOut'] for m in before['design']['measurements'])
assert all(e['annotations'].get('condition')!='sgrna_clu' for e in before['design']['entities'] if e['id'] in before['design']['measurements'][0]['entityIDs'])
# Re-read original GEF rows independently of preparation; all admitted cells/genes.
lookup={n:i for i,n in enumerate(x.var_names)}
for source in lineage['files']:
    file=next(a.sources.glob(source['id']+'*.gef'));assert sha(file)==source['sha256']
    with h5py.File(file) as f:
        genes=[n.decode() for n in f['cellBin/gene']['geneName']];cells=f['cellBin/cell'][:]
        for i,row in enumerate(lineage['mapping']):
            if row['donor']!=source['id']:continue
            c=cells[row['sourceIndex']];assert int(c['id'])==row['sourceCellID'];assert [int(c['x']),int(c['y'])]==[row['x'],row['y']]
            values=f['cellBin/cellExp'][int(c['offset']):int(c['offset'])+int(c['geneCount'])];expected={}
            for entry in values:
                name=genes[entry['geneID']]
                if name in lookup:expected[lookup[name]]=expected.get(lookup[name],0)+int(entry['count'])
            native={int(j):int(v) for j,v in zip(x.X[i].indices,x.X[i].data)};assert native==expected
bulk=read(root/'source/report.json')['pseudobulk'];m=bulk['matrix'];counts=np.zeros((len(bulk['groups']),len(bulk['featureIDs'])),dtype=np.uint64)
for i,g in enumerate(bulk['groups']):
    counts[i,m['featureIndices'][m['rowOffsets'][i]:m['rowOffsets'][i+1]]]=m['counts'][m['rowOffsets'][i]:m['rowOffsets'][i+1]]
    np.testing.assert_array_equal(counts[i],np.asarray(x.X[g['sourceCellIndices']].sum(axis=0)).ravel())
logs=np.log1p(counts/counts.sum(axis=1)[:,None]*1e6);groups=bulk['groups'];train=assay['trainingUnits']
c=np.array([next(i for i in assay['trainingRows'] if groups[i]['donorID']==d and groups[i]['condition']=='sgrna_msafe') for d in train]);tr=np.array([next(i for i in assay['trainingRows'] if groups[i]['donorID']==d and groups[i]['condition']=='sgrna_clu') for d in train]);response=logs[tr]-logs[c];mean=response.mean(axis=0);pair=counts[c]+counts[tr];selected=(pair.sum(axis=0)>=10)&((pair>0).sum(axis=0)>=2);context=logs[c][:,selected];center=context.mean(axis=0);scale=context.std(axis=0);constant=np.all(context==context[0],axis=0);center[constant]=context[0,constant];scale[constant]=1;z=(context-center)/scale/np.sqrt(selected.sum());oracles={}
for region in assay['regions']:
    q=region['controlRow'];query=(logs[q,selected]-center)/scale/np.sqrt(selected.sum());ridge=mean+(query@z.T)@np.linalg.solve(z@z.T+np.eye(len(c)),response-mean)
    effects={'contextRidge':ridge,'noChange':np.zeros(len(mean)),'meanResponse':mean,'medianResponse':np.median(response,axis=0)}
    pred=next(r for r in before['regionalPredictions'] if r['regionID']==region['id']);oracles[region['id']]={}
    for e in pred['estimates']:
        expected=np.maximum(0,logs[q]+effects[e['baseline']]);np.testing.assert_allclose(e['predictedTreated'],expected,rtol=1e-9,atol=1e-9);oracles[region['id']][e['baseline']]=expected
# First held-out scoring occurs only after predictions and instrument comparisons.
result=adapter.reveal(a.run,{})
for region in result['regions']:
    observedRow=next(r['observedRow'] for r in assay['regions'] if r['id']==region['regionID'])
    for metric in region['metrics']:
        expected=np.sqrt(np.mean((oracles[region['regionID']][metric['model']]-logs[observedRow])**2));np.testing.assert_allclose(metric['rmse'],expected,rtol=1e-9,atol=1e-9)
adapter.verify(a.run,{})
selection=read(a.run/'registration.json')['plan'];base={'specimen':assay['specimen'],'target':'Clu','regionIDs':selection['regionIDs'],'timepoint':'study-endpoint'}
rejected=[]
for key,value in [('target','Srf'),('timepoint','48h'),('regionIDs',['region-0-0']),('dose',1)]:
    request={**base,key:value}
    try:adapter.compile(a.assay,request)
    except ValueError:rejected.append(key)
    else:raise AssertionError('unsupported '+key)
for name in ['prediction/folds/000/prediction.json','specimen.json','comparison.json']:
    dest=a.output/('tamper-'+name.split('/')[0]);shutil.copytree(a.run,dest);(dest/name).write_text('{}')
    try:adapter.summary(dest)
    except (ValueError,KeyError):rejected.append(name)
    else:raise AssertionError('tamper '+name)
write(a.output/'checks.json',{'passed':True,'cells':x.n_obs,'features':x.n_vars,'heldOutBiologicalUnits':1,'admittedRegions':len(assay['regions']),'excludedRegions':len(assay['excludedRegions']),'sourceCountsExact':True,'ridgeAndBaselinesMatch':True,'beforeRevealNoTargetValues':True,'exactReplay':True,'rejections':rejected,'result':{k:v for k,v in result.items() if k not in ['regions','geneRMSE']}})
print(json.dumps({'passed':True,'metrics':result['metrics'],'verdict':result['verdict']}))
