#!/usr/bin/env python3
"""Multi-target GEF ingestion retaining ALL measured cells and neighborhood RNA.

No prediction or effect scoring happens here. Target support uses guide counts
only. Chip/animal pooling is resolved from supplementary data, never file count.
"""
import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,h5py,anndata as ad
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree
from wetlab import sha,read,write,require

SOURCE_ROOT=Path(__file__).parent
# Supplementary Data 3 cell totals uniquely match deposits to paper chip IDs.
IDENTITIES={
 'GSM8449354':{'chip':'chip1','animals':['paper-mouse1'],'sections':2,'sectionAssignment':'unresolved per cell; same animal','split':'test','expectedCells':39862},
 'GSM9659883':{'chip':'chip2','animals':['paper-mouse2','paper-mouse3'],'sections':None,'sectionAssignment':'pooled; Data 3 says two, Fig S7 says three; do not split','split':'training','expectedCells':123921},
 'GSM9659884':{'chip':'chip3','animals':['paper-mouse4'],'sections':2,'sectionAssignment':'unresolved per cell; same animal','split':'validation','expectedCells':65992}}
# Preserve exact deposited barcode labels. Do not silently repair ambiguous aliases.
CANONICAL={'sgrna_C9orf72':'C9orf72','sgrna_Cfap410':'Cfap410','sgrna_clu':'Clu','sgrna_fasn':'Fasn','sgrna_flcn':'Flcn','sgrna_gfap':'Gfap','sgrna_lrrk2':'Lrrk2','sgrna_rraga':'Rraga','sgrna_sh3gl2':'Sh3gl2','sgrna_srf':'Srf','sgrna_stk39':'Stk39','sgrna_tbk1':'Tbk1','sgrna_trem2':'Trem2','sgrna_msafe':'mSafe'}
# Frozen marker panel from study Figure 1; annotations are uncertain inferences.
MARKERS={'neuron':['Snap25','Syt1','Rbfox3'],'oligodendrocyte':['Mbp','Plp1','Mog'],'astrocyte':['Gfap','Aqp4','Slc1a3'],'microglia':['C1qa','C1qb','P2ry12'],'endothelial':['Flt1','Kdr','Cldn5'],'choroid-plexus':['Ttr','Klotho'],'erythroid':['Hba-a1','Hbb-bs']}

def protocol():
 return {'format':'numivivo-spatial-response-study/v1','source':'GSE274447','identities':IDENTITIES,
 'sourceIdentityEvidence':'Supplementary Data 3 cell totals and section/animal row; Fig S7; GEO titles retained as conflicting labels',
 'eligibility':{'minControlCellsAppliesTo':'isolated mSafe controls; insufficient controls forbid biological qualification, not exploratory execution','minimumUniqueGuideCells':{'training':10,'validation':5,'test':5},'canonicalTargetRequired':True,'minControlCells':5,'excludedPriorRevealedTarget':'Clu remains development/regression only'},
 'controls':'mSafe positive only for direct reference; barcode-negative neighbors of mSafe are separately labelled control-associated, not proven untreated; exclude ALL reference cells including mSafe within 15-neighbor sets of non-mSafe or ambiguous guide anchors',
 'geometry':'fixed endpoint GEF coordinates; no anatomical or per-animal section assignment inferred',
 'neighbors':15,'referenceAnchors':3,'referenceCoverage':'nearest admitted control within ten times the chip median 15-neighbor radius; otherwise annotation uncertain','annotation':{'markers':MARKERS,'input':'reference controls only; projected onto query cells using geometry','uncertainty':'reference-label probabilities; uncalibrated; not measured cell types'},
 'features':'union of named non-guide RNA; source presence mask; missing is not zero; all measured RNA retained',
 'splitPolicy':'pooled chip2 training; chip3 validation; chip1 locked test; no animal crosses splits; no within-chip independence claim',
 'preprocessing':'log1p(CPM), training-only variance-ranked 256-gene encoder panel, training-reference PCA16; whole-gene masked decoder and scores',
 'models':['native-ridge-v02','no-change','training-mean','cell-type-matched','learned-no-neighborhood','learned-neighborhood','learned-shuffled-neighborhood','learned-frozen-Nicheformer'],
 'learning':{'hiddenWidth':64,'batchSize':16,'learningRate':0.001,'weightDecay':0.0001,'seed':314159,'checkpointSteps':[240,720,1440],'selection':'validation equal-target mean response RMSE; no test model selection','loss':'Gaussian population negative log likelihood in log1p(CPM); uncertainty not calibrated'},
 'unseenPerturbations':'separate leave-one-target-out run for first eligible non-Clu target in lexicographic order; no identity-specific response prior or embedding at inference',
 'metrics':['whole-gene response RMSE','response MAE','direction agreement outside absolute 0.1 log1p(CPM) deadband','effect-size bias','Gaussian NLL','1D Wasserstein on training-selected 64-gene panel'],
 'reportAxes':['target','projected cell type','reference neighborhood','role direct/neighbor','chip/animal pool'],
 'promotion':'not established with this small pooled cohort; spatial benefit requires improvement over no-neighborhood and shuffle for each eligible held-out target, role and adequately supported population; no genes/cells used as independent animals',
 'performanceTargets':{'specimenCells':123921,'warmFeatureSwitchP95Ms':100,'pickP95Ms':20,'panFrameP95Ms':33.4},
 'stoppingRule':'one frozen campaign; any code/instrument failure retained; no outcome-driven eligibility or test tuning'}

def ingest(source,output):
 output=Path(output);output.mkdir(parents=True,exist_ok=False);write(output/'protocol.json',protocol())
 manifest=read(SOURCE_ROOT/'spatial_perturb_sources.json');axis=set();panels={}
 for entry in manifest:
  file=Path(source)/entry['file'];require(sha(file)==entry['sha256'],'Pinned GEF mismatch')
  with h5py.File(file) as f:names=[x.decode() for x in f['cellBin/gene']['geneName']]
  kept=[x for x in names if x and not x.startswith('sgrna_')];require(len(set(kept))==len(kept),'Duplicate named gene');axis.update(kept);panels[entry['file']]=names
 features=sorted(axis);lookup={g:i for i,g in enumerate(features)};types=list(MARKERS);summary=[];supports={};source_receipts=[]
 for entry in manifest:
  file=Path(source)/entry['file'];acc=file.name.split('_')[0];identity=IDENTITIES[acc];names=panels[file.name]
  with h5py.File(file) as f:
   cells=f['cellBin/cell'][:];require(len(cells)==identity['expectedCells'],'Chip identity cell count mismatch')
   exp=f['cellBin/cellExp'][:];offsets=np.r_[cells['offset'].astype(np.int64),len(exp)];require(np.array_equal(np.diff(offsets),cells['geneCount']),'GEF contiguous row offsets')
   raw=csr_matrix((exp['count'].astype(np.uint32),exp['geneID'].astype(np.int32),offsets),shape=(len(cells),len(names)))
  bars=sorted(n for n in names if n.startswith('sgrna_'));bc=raw[:,[names.index(n) for n in bars]].toarray();positive=(bc>0).sum(axis=1)
  labels=np.array(['barcode-negative']*len(cells),dtype=object);labels[positive>1]='ambiguous-multiple';one=positive==1;labels[one]=np.asarray(bars)[bc[one].argmax(axis=1)]
  source_columns=np.array([i for i,n in enumerate(names) if n in lookup]);local=raw[:,source_columns].tocsr();local.indices=np.asarray([lookup[names[source_columns[j]]] for j in local.indices],np.int32);local._shape=(len(cells),len(features));local.sum_duplicates();local.sort_indices()
  del raw,exp
  xy=np.c_[cells['x'],cells['y']].astype(np.float64);tree=cKDTree(xy);dist,neighbors=tree.query(xy,k=16,workers=2);keep=neighbors!=np.arange(len(cells))[:,None]; neighbors=np.asarray([row[mask][:15] for row,mask in zip(neighbors,keep)],np.int32);dist=np.asarray([row[mask][:15] for row,mask in zip(dist,keep)])
  safe=np.flatnonzero(labels=='sgrna_msafe');noncontrol=np.flatnonzero((positive>0)&(labels!='sgrna_msafe'))
  contaminated=np.zeros(len(cells),bool);contaminated[neighbors[noncontrol].ravel()]=True
  reference=np.zeros(len(cells),bool);reference[safe]=True;reference[neighbors[safe].ravel()]=True;reference &= ((labels=='sgrna_msafe')|(positive==0)) & ~contaminated
  # Annotate only admitted reference cells. Query annotations are geometry-based
  # projections, never classifiers run on held-out perturbed expression.
  refs=np.flatnonzero(reference);library=np.asarray(local[refs].sum(axis=1)).ravel();marker_scores=[]
  for kind in types:
   ix=[lookup[g] for g in MARKERS[kind] if g in lookup];counts=local[refs][:,ix].toarray();marker_scores.append(np.log1p(counts/np.maximum(1,library[:,None])*1e6).mean(axis=1))
  scores=np.asarray(marker_scores).T;prob=np.exp(scores-scores.max(axis=1,keepdims=True));prob/=prob.sum(axis=1,keepdims=True)
  ref_tree=cKDTree(xy[refs]);rd,ri=ref_tree.query(xy,k=min(3,len(refs)));weights=1/np.maximum(rd,1);weights/=weights.sum(axis=1,keepdims=True);projected=(prob[ri]*weights[:,:,None]).sum(axis=1);confidence=projected.max(axis=1);ctype=np.asarray(types,dtype=object)[projected.argmax(axis=1)];coverage=rd[:,0]<=10*np.median(dist[:,-1]);ctype[(confidence<.45)|~coverage]='uncertain'
  obs=pd.DataFrame({'source_cell_id':cells['id'],'accession':acc,'chip':identity['chip'],'animal_pool':'|'.join(identity['animals']),'section':'unresolved-within-'+identity['chip'],'partition':identity['split'],'guide_assignment':labels,'guide_species_count':positive,'reference_admitted':reference,'projected_cell_type':ctype,'annotation_confidence':confidence,'reference_coverage':coverage,'nearest_reference_distance':rd[:,0],'barcode_negative_is_untreated':False},index=[acc+':'+str(i) for i in cells['id']])
  measured=np.array([g in names for g in features]);var=pd.DataFrame({'measured_in_source':measured},index=features)
  a=ad.AnnData(local,obs=obs,var=var);a.obsm['spatial']=xy;a.obsm['guide_counts']=bc;a.obsm['projected_type_probabilities']=projected.astype(np.float32);a.obsm['neighbor_indices']=neighbors;a.obsm['neighbor_distances']=dist.astype(np.float32);a.obsm['reference_indices']=refs[ri].astype(np.int32);a.obsm['reference_weights']=weights.astype(np.float32);a.uns['guide_ids']=bars;a.uns['projected_type_names']=types;a.uns['coordinate_unit']='GEF-coordinate';a.uns['identity_json']=json.dumps(identity);a.uns['source_sha256']=entry['sha256'];a.uns['evidence']='RNA MEASURED; guide assignment and projected cell type MODEL INFERENCE; negative barcode exposure UNKNOWN'
  a.write_h5ad(output/(identity['chip']+'.h5ad'),compression='gzip')
  support={g:int((labels==g).sum()) for g in bars};supports[identity['split']]=support
  summary.append({'accession':acc,**identity,'sourceSHA256':entry['sha256'],'cells':len(cells),'measuredGenes':int(measured.sum()),'barcodes':support,'barcodeNegative':int((positive==0).sum()),'ambiguous':int((positive>1).sum()),'isolatedSafeControls':int(((labels=='sgrna_msafe')&reference).sum()),'referenceCells':len(refs),'projectedTypes':{str(t):int((ctype==t).sum()) for t in np.unique(ctype)}})
  print(identity['chip'],len(cells),'cells',len(refs),'reference cells',flush=True)
  del a,local
 targets=[]
 for guide in sorted(set.union(*(set(s) for s in supports.values()))-{'sgrna_msafe'}):
  reasons=[]
  if guide not in CANONICAL:reasons.append('canonical target alias unresolved in deposit')
  for split,minimum in protocol()['eligibility']['minimumUniqueGuideCells'].items():
   if supports[split].get(guide,0)<minimum:reasons.append(split+' below guide support minimum')
  targets.append({'guideID':guide,'target':CANONICAL.get(guide),'counts':{s:v.get(guide,0) for s,v in supports.items()},'eligible':not reasons,'exclusions':reasons,'previouslyExposed':guide=='sgrna_clu'})
 write(output/'cohort.json',{'format':'numivivo-spatial-cohort/v1','protocolSHA256':sha(output/'protocol.json'),'builderSHA256':sha(__file__),'sources':summary,'targets':targets,'featureCount':len(features),'files':{p.name:sha(p) for p in output.glob('*.h5ad')},'biologicalUnits':4,'independentSplitBlocks':3,'limitations':['chip2 animals remain pooled','section boundaries unresolved; neighbors are fixed chip-coordinate KNN, not validated section assignments','controls are sparse','query types and reference neighborhoods are inferred solely from control-associated RNA; not measured pre-intervention state']})
 return output
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();ingest(a.sources,a.output)
