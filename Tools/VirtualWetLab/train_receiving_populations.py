#!/usr/bin/env python3
"""Test the receiver-reference hypothesis on permanently exposed spatial data.

Original nearest-coordinate associations are retained as hypotheses. They are
not promoted to verified tissue edges: the source lacks cell-to-section identity.
The spatial-edge and edge-shuffle comparisons therefore fail admission explicitly.
"""
import argparse
from pathlib import Path
import numpy as np
import anndata as ad
from safetensors.numpy import load_file,save_file
from wetlab import read,write,sha,require,timestamp,inventory
from prepare_intervention_design import orthology,axis,normalized,V03
from train_intervention_design import native,validation_score


def run(inputs,campaign,binary,out):
    require(not out.exists(),'Receiver experiment already exists');out.mkdir();features=read(inputs/'features.json');targetrecords=read(inputs/'targets.json');targets=sorted(targetrecords);mapping=orthology(V03/'research/MGI-HOM_MouseHumanSequence.rpt')
    write(out/'registration.json',{'createdAt':timestamp(),'protocolSHA256':sha(inputs/'protocol.json'),'scope':'EXPOSED SPATIAL DEVELOPMENT; no new generalization evidence',
      'hypothesis':'Receiver-specific control state improves matched population endpoints versus the original sender control state',
      'variants':['receiver-pretrained','anchor-reference-pretrained','receiver-spatial-only'],'selection':'validation chip3 equal-target RMSE, fixed three checkpoints',
      'relationshipGate':{'verifiedSectionIdentity':False,'edgeModel':'UNAVAILABLE','edgeShuffle':'UNAVAILABLE','reason':'Source does not bind individual cells to sections. Do not invent biological edges across unresolved boundaries.'},
      'outcomeAssociations':'Original v0.3 15-nearest chip-coordinate membership is retained for development comparison only. No physical-distance or propagation claim.',
      'receivers':'Separate each reference-projected cell type; barcode-negative cells are not assumed unperturbed; remove detected coexposure to another noncontrol target',
      'pretrainedWeights':str(campaign/'target-descriptor'),'ownerSHA256':sha(__file__)})
    arrays={};metadata={};geometry=[]
    original=read(V03/'cohort-v3/cohort.json');eligible={x['guideID']:x['target'] for x in original['targets'] if x['eligible']}
    for chip in ('chip2','chip3','chip1'):
        a=ad.read_h5ad(V03/'cohort-v3'/f'{chip}.h5ad');ax=axis(a,mapping,'mouse');mask=np.asarray([g in ax and bool(a.var.measured_in_source.iloc[ax[g]]) for g in features],np.float32);cols=[ax.get(g,0) for g in features]
        labels=a.obs.guide_assignment.astype(str).to_numpy(dtype=str);types=a.obs.projected_cell_type.astype(str).to_numpy(dtype=str);refs=np.flatnonzero(a.obs.reference_admitted.to_numpy());lookup={int(r):i for i,r in enumerate(refs)};refx=normalized(a,refs,cols)*mask;ri=a.obsm['reference_indices'];rw=a.obsm['reference_weights'];neighbors=a.obsm['neighbor_indices'];ri_local=np.asarray([[lookup[int(j)] for j in r] for r in ri]);reference=(refx[ri_local]*rw[:,:,None]).sum(1)
        entries=[];own=[];anchorctx=[];descriptors=[];tids=[];ys=[];vs=[];contamination={}
        for guide,target in eligible.items():
            canonical=mapping.get(target)
            if canonical not in targets:continue
            anchors=np.flatnonzero(labels==guide);receivers=np.unique(neighbors[anchors].ravel());receivers=receivers[labels[receivers]=='barcode-negative']
            clean=[]
            for r in receivers:
                local=set(labels[neighbors[r]])-{'barcode-negative','sgrna_msafe',guide}
                if not local:clean.append(int(r))
            contamination[target]={'candidateReceivers':len(receivers),'retained':len(clean),'excludedDetectedCoexposure':len(receivers)-len(clean)}
            for role,cells in [('direct',anchors),('receiver',np.asarray(clean,int))]:
                for celltype in sorted(set(types[cells])):
                    rows=cells[types[cells]==celltype]
                    if len(rows)<5:continue
                    # A separate response distribution and reference for each
                    # receiving type; no pooled heterogeneous neighbor target.
                    control=reference[rows].mean(0);sender=reference[anchors].mean(0)
                    y=normalized(a,rows,cols)*mask
                    d=np.r_[targetrecords[canonical]['descriptor'],1.,0.,0.,1.,0.,0.,0.,float(role=='direct'),float(role=='receiver')].astype(np.float32)
                    own.append(control);anchorctx.append(sender);descriptors.append(d);tids.append(targets.index(canonical));ys.append(y.mean(0));vs.append(y.var(0))
                    entries.append({'source':chip,'target':canonical,'role':role,'receivingType':celltype,'outcomeRows':rows.tolist(),
                       'referenceRows':np.unique(ri[rows]).tolist(),'controlCells':len(np.unique(ri[rows])),'cells':len(rows),'sourceSHA256':sha(V03/'cohort-v3'/f'{chip}.h5ad'),
                       'annotationConfidence':float(a.obs.annotation_confidence.iloc[rows].mean()),'referenceCoverage':float(a.obs.reference_coverage.iloc[rows].mean()),
                       'identityEvidence':'control-projected uncertain identity, not measured pre-intervention state','verifiedSection':None,'coordinateUnits':'uncalibrated original chip coordinate',
                       'detectedCoexposure':contamination[target],'xy':a.obsm['spatial'][rows].astype(float).tolist()})
        n=len(own);ar={'context':np.asarray(own,np.float32),'descriptor':np.asarray(descriptors,np.float32),'target':np.asarray(tids,np.int32),'known':np.zeros((n,1),np.float32),'observed':np.asarray(ys,np.float32),'observedVariance':np.asarray(vs,np.float32),'mask':np.tile(mask,(n,1)),'stratum':np.arange(n,dtype=np.int32),'prior':np.zeros((len(targets),len(features)),np.float32)}
        arrays[chip]=(ar,np.asarray(anchorctx,np.float32));metadata[chip]=entries;write(out/(chip+'-rows.json'),entries);write(out/(chip+'-coexposure.json'),contamination)
        print(chip,n,'separate target/receiving-type distributions',flush=True);del a,reference,refx
    plan=read(inputs/'plan.json');plan['descriptorCount']=158
    selected=read(campaign/'target-descriptor/selection.json')['selected']['step'];weights=load_file(str(campaign/'target-descriptor/training'/f'weights-{selected}.safetensors'))
    key='descriptorProjection.weight';require(weights[key].shape[1]==156,'Pretrained descriptor identity');weights[key]=np.pad(weights[key],((0,0),(0,2)));save_file(weights,str(out/'warm-start.safetensors'))
    results={};predictions={};common=np.logical_and.reduce([v[0]['mask'][0].astype(bool) for v in arrays.values()]);require(common.any(),'No common measured spatial features');write(out/'scored-features.json',[g for g,keep in zip(features,common) if keep])
    for variant in ('receiver-pretrained','anchor-reference-pretrained','receiver-spatial-only'):
        folder=out/variant;folder.mkdir();write(folder/'plan.json',plan)
        for chip,(raw,anchor) in arrays.items():
            ar={k:v.copy() for k,v in raw.items()}
            if variant=='anchor-reference-pretrained':ar['context']=anchor.copy()
            if chip=='chip1':ar={k:v for k,v in ar.items() if k in ('context','descriptor','target','known')}
            save_file(ar,str(folder/(chip+'.safetensors')))
        native(binary,'train',folder/'plan.json',folder/'chip2.safetensors',folder/'training',None if variant=='receiver-spatial-only' else out/'warm-start.safetensors')
        scores=[]
        for step in plan['steps']:
            dest=folder/f'validation-{step}';native(binary,'predict',folder/'plan.json',folder/'chip3.safetensors',dest,folder/'training'/f'weights-{step}.safetensors');p=load_file(str(dest/'prediction.safetensors'))
            scores.append({'step':step,'score':validation_score(p['mean'][:,common],arrays['chip3'][0]['observed'][:,common],metadata['chip3'])})
        best=min(scores,key=lambda x:(x['score'],x['step']));write(folder/'selection.json',{'scores':scores,'selected':best})
        native(binary,'predict',folder/'plan.json',folder/'chip1.safetensors',folder/'development-prediction',folder/'training'/f"weights-{best['step']}.safetensors")
        p=load_file(str(folder/'development-prediction/prediction.safetensors'));predictions[variant]=p['mean'];error=np.sqrt(np.mean((p['mean'][:,common]-arrays['chip1'][0]['observed'][:,common])**2,axis=1));results[variant]={'groupRMSE':error.tolist(),'meanGroupRMSE':float(error.mean()),'selected':best};print(variant,results[variant]['meanGroupRMSE'],flush=True)
    reference=arrays['chip1'][0];results['no-change']={'groupRMSE':np.sqrt(np.mean((reference['context'][:,common]-reference['observed'][:,common])**2,axis=1)).tolist()};results['no-change']['meanGroupRMSE']=float(np.mean(results['no-change']['groupRMSE']))
    save_file({**predictions,'control':reference['context'],'observed':reference['observed'],'observedVariance':reference['observedVariance']},str(out/'development-results.safetensors'))
    write(out/'comparison.json',{'scope':'permanently exposed spatial development','models':results,'groups':metadata['chip1'],'biologicalPromotion':False,'spatialInformationBenefit':'UNAVAILABLE: section identity not resolved; no edge model executed','receiverImprovement':results['receiver-pretrained']['meanGroupRMSE']<results['anchor-reference-pretrained']['meanGroupRMSE'],'pretrainingImprovement':results['receiver-pretrained']['meanGroupRMSE']<results['receiver-spatial-only']['meanGroupRMSE']})
    write(out/'artifact-seal.json',{'createdAt':timestamp(),'files':inventory(out),'binarySHA256':sha(binary)})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.inputs,a.campaign,a.binary,a.output)
