#!/usr/bin/env python3
"""Source-bound GSE274447 ingestion; author a real spatial perturbation assay.

Counts and coordinates come from the deposited GEF files. The native Omics
owner aggregates and predicts. This script performs format/annotation mapping.
"""
import argparse, json
from pathlib import Path
import h5py, numpy as np, pandas as pd, anndata as ad
from scipy.sparse import csr_matrix, vstack
from wetlab import write,read,sha,require,native

SOURCE_IDS=['GSM8449354','GSM9659883','GSM9659884']
TARGET='sgrna_clu';CONTROL='sgrna_msafe'

def extract(path):
    with h5py.File(path) as f:
        cells=f['cellBin/cell'][:];genes=f['cellBin/gene'][:]
        names=[n.decode() for n in genes['geneName']]
        require(len(set(n for n in names if n))==sum(bool(n) for n in names),'Duplicate named GEF gene identities')
        barcodes=[i for i,n in enumerate(names) if n.startswith('sgrna_')]
        count=np.zeros(len(cells),dtype=np.uint16);which=np.full(len(cells),-1,dtype=np.int32)
        id_lookup={int(c['id']):i for i,c in enumerate(cells)}
        for j in barcodes:
            g=genes[j];rows=f['cellBin/geneExp'][int(g['offset']):int(g['offset'])+int(g['cellCount'])]
            for row in rows:
                if row['count']>0:
                    i=id_lookup[int(row['cellID'])];count[i]+=1;which[i]=j
        selected=np.flatnonzero((count==1)&np.isin(which,[names.index(TARGET),names.index(CONTROL)]))
        offsets=[0];columns=[];values=[]
        for i in selected:
            c=cells[i];row=f['cellBin/cellExp'][int(c['offset']):int(c['offset'])+int(c['geneCount'])]
            columns.extend(row['geneID'].tolist());values.extend(row['count'].tolist());offsets.append(len(values))
        x=csr_matrix((np.array(values,dtype=np.uint32),columns,offsets),shape=(len(selected),len(names)));x.sum_duplicates();x.sort_indices()
        obs=pd.DataFrame({'sourceIndex':selected,'sourceCellID':cells['id'][selected],
            'x':cells['x'][selected],'y':cells['y'][selected],'condition':[names[which[i]] for i in selected]})
        # Background is geometry only, selected by fixed stride, with exact source IDs.
        background=cells[::max(1,len(cells)//2400)]
        return names,x,obs,background,{'unnamedFeaturesExcluded':sum(not n for n in names),'allCells':len(cells),'selected':len(selected),'ambiguousBarcodes':int((count>1).sum()),
            'bounds':[int(cells['x'].min()),int(cells['x'].max()),int(cells['y'].min()),int(cells['y'].max())],'resolutionAttribute':int(f.attrs['resolution'][0]),'offsetX':int(f.attrs['offsetX'][0]),'offsetY':int(f.attrs['offsetY'][0])}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True)
    p.add_argument('--binary',type=Path,required=True);p.add_argument('--tissue-binary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    manifest=read(Path(__file__).with_name('spatial_perturb_sources.json'))
    files=[a.sources/entry['file'] for entry in manifest]
    require([entry['file'].split('_')[0] for entry in manifest]==SOURCE_IDS,'Source accession order')
    for file,entry in zip(files,manifest):require(sha(file)==entry['sha256'],'Source differs from pinned public GEF: '+file.name)
    # The analysis and selection rules are frozen before native predictions or scores.
    protocol={'id':'spatial-clu-heldout-mouse3-v1','target':TARGET,'control':CONTROL,'trainingUnits':SOURCE_IDS[:2],
        'heldOutUnit':SOURCE_IDS[2],'barcodeRule':'exactly one positive gRNA species; count > 0',
        'regionRule':'2 by 2 grid on full held-out specimen coordinate bounds; no anatomical alignment',
        'minimumCellsPerArm':3,'featureRule':'intersection of nonempty measured gene IDs; exclude unnamed features and gRNA barcodes; no expression-based selection',
        'primaryModel':'contextRidge','primaryMetric':'equal-region all-common-gene log1p(CPM) response RMSE',
        'successCriterion':'primary RMSE strictly lower than both noChange and meanResponse, for every admitted region',
        'stoppingRule':'one frozen run; retain failures; no tuning after held-out reveal',
        'purpose':'retrospective development benchmark; three mice, cells and regions are subsamples',
        'limits':['No same-cell pretreatment measurements','Mixed cell composition and sparse safe-harbour controls confound effects',
            'Barcode assignment is inferred from measured barcode RNA; functional knockout per cell unverified',
            'Known intervention transfer; not unseen-context or unseen-perturbation qualification']}
    write(a.output/'protocol.json',protocol)
    raw=[extract(f) for f in files];features=sorted(set.intersection(*(set(r[0]) for r in raw))-{n for r in raw for n in r[0] if n.startswith('sgrna_') or not n})
    matrices=[];obs=[];samples=[];provenance=[];entities=[]
    for k,(file,(names,x,frame,bg,meta)) in enumerate(zip(files,raw)):
        sid=SOURCE_IDS[k];lookup={n:i for i,n in enumerate(names)};x=x[:,[lookup[n] for n in features]].tocsr();matrices.append(x)
        frame=frame.copy();frame['donor']=sid;frame['sample']=sid+'__'+frame.condition;frame['cellGroup']='whole'
        if k==2:
            xmin,xmax,ymin,ymax=meta['bounds']
            frame['cellGroup']=['region-'+str(min(1,int((xx-xmin)*2/max(1,xmax-xmin))))+'-'+str(min(1,int((yy-ymin)*2/max(1,ymax-ymin)))) for xx,yy in zip(frame.x,frame.y)]
        frame.index=[sid+':'+str(i) for i in frame.sourceCellID];obs.append(frame)
        for condition in (CONTROL,TARGET):samples.append({'id':sid+'__'+condition,'biologicalReplicateID':sid,'donorID':sid,'condition':condition,'batchID':sid,'organism':'NCBITaxon:10090'})
        provenance.append({'id':sid,'uri':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+sid,'sha256':sha(file),'description':'Deposited Stereo-seq cell-segmented GEF; '+json.dumps(meta,sort_keys=True)})
        if k==2:
            evidence={'state':'MEASURED','sourceID':sid,'modelID':None,'assumptions':['GEF segmentation-derived cell centroids; no histology image'], 'validationDomain':None,'uncertainty':'segmentation error not quantified','heldOut':False}
            entities=[{'id':'specimen','level':'specimen','parentID':None,'memberIDs':[],'position':None,'biologicalUnitID':sid,'annotations':{},'evidence':evidence}]
            for region in sorted(set(frame.cellGroup)):
                entities.append({'id':region,'level':'region','parentID':'specimen','memberIDs':[], 'position':None,'biologicalUnitID':sid,'annotations':{'definition':'coordinate-grid, not anatomical annotation'},'evidence':{**evidence,'state':'HYPOTHESIS','assumptions':['Authored coordinate grid']}})
            rows={int(r.sourceCellID):r for _,r in frame.iterrows()}
            for c in sorted({int(c['id']):(int(c['x']),int(c['y'])) for c in bg}.items()):
                if c[0] not in rows:
                    entities.append({'id':sid+':'+str(c[0]),'level':'cell','parentID':'specimen','memberIDs':[],'position':list(c[1]),'biologicalUnitID':sid,'annotations':{'role':'geometry-background','sourceCellID':str(c[0])},'evidence':evidence})
            for cell,r in frame.iterrows():
                entities.append({'id':cell,'level':'cell','parentID':r.cellGroup,'memberIDs':[],'position':[float(r.x),float(r.y)],'biologicalUnitID':sid,'annotations':{'condition':r.condition,'sourceCellID':str(r.sourceCellID),'sourceIndex':str(r.sourceIndex),'cellType':'unannotated','interventionAssignment':'barcode-derived hypothesis'},'evidence':evidence})
    data=vstack(matrices,format='csr');observations=pd.concat(obs);adata=ad.AnnData(X=data,obs=observations[['sample','cellGroup']],var=pd.DataFrame(index=features));adata.write_h5ad(a.output/'counts.h5ad',compression='gzip')
    write(a.output/'source-plan.json',{'schemaVersion':1,'mapping':{'schemaVersion':1,'id':'GSE274447-clu','evidence':'measured','sourceDescription':'GSE274447 exact integer counts for uniquely Clu/mSafe barcode-positive cells, shared gene intersection. Regions remain within mouse biological units.','countUnit':'umiCount','matrixPath':'X','samples':samples,'sampleColumn':'sample','groupColumn':'cellGroup'},'contrasts':[]})
    native(a.binary,['singlecell-h5ad-pseudobulk',a.output/'counts.h5ad','--plan',a.output/'source-plan.json','--output',a.output/'source'],a.output/'native-ingest-log.json')
    bulk=read(a.output/'source/report.json')['pseudobulk'];groups=bulk['groups'];regions=[];excluded=[]
    for region in sorted(set(obs[2].cellGroup)):
        rows={g['condition']:i for i,g in enumerate(groups) if g['donorID']==SOURCE_IDS[2] and g['cellGroup']==region}
        counts={c:len(groups[rows[c]]['sourceCellIndices']) if c in rows else 0 for c in (CONTROL,TARGET)}
        entry={'id':region,'counts':counts}
        if min(counts.values())>=3:regions.append({**entry,'controlRow':rows[CONTROL],'observedRow':rows[TARGET]})
        else:excluded.append({**entry,'status':'UNAVAILABLE','reason':'fewer than 3 control or target cells'})
    require(regions,'No region has adequate support; do not lower the registered threshold')
    # Shared representation keeps measured sparse RNA separate from visual geometry.
    q=matrices[2];qids=obs[2].index.tolist()
    measurement={'id':'endpoint-rna','modality':'rna','unit':'umiCount','timepoint':'study-endpoint','entityIDs':qids,'featureIDs':features,'rowOffsets':q.indptr.tolist(),'featureIndices':q.indices.tolist(),'values':q.data.tolist(),'absentValue':'zero','evidence':{'state':'MEASURED','sourceID':SOURCE_IDS[2],'modelID':None,'assumptions':['Integer counts restricted to common gene universe; barcode features excluded'], 'validationDomain':None,'uncertainty':'sampling noise and segmentation uncertainty not quantified','heldOut':True}}
    specimen={'format':'numivivo-tissue-specimen/v1','id':SOURCE_IDS[2],'title':'Mouse brain · Clu perturbation','organism':'NCBITaxon:10090','coordinateSystem':{'id':SOURCE_IDS[2]+'-chip','unit':'GEF-coordinate','axes':['x','y'],'description':'Original GEF coordinates; no cross-section registration. Resolution attribute retained in source metadata.','micrometresPerUnit':None},'sources':provenance,'entities':entities,'measurements':[measurement],'limitations':protocol['limits']+['Background geometry is a fixed-stride display sample; all eligible intervention/control cells retained','Only a destructive endpoint is available; no temporal interpolation']}
    write(a.output/'specimen.json',specimen);native(a.tissue_binary,[a.output/'specimen.json'],a.output/'tissue-validation.json')
    # Public design contains coordinates and control RNA only; target values are reveal-only.
    controlRows=[i for i,r in obs[2].iterrows() if r.condition==CONTROL];ix=[qids.index(i) for i in controlRows];qc=q[ix].tocsr()
    public={**specimen,'measurements':[{**measurement,'entityIDs':controlRows,'rowOffsets':qc.indptr.tolist(),'featureIndices':qc.indices.tolist(),'values':qc.data.tolist(),'evidence':{**measurement['evidence'],'heldOut':False}}]}
    write(a.output/'design-specimen.json',public)
    write(a.output/'lineage.json',{'files':provenance,'mapping':observations.reset_index().rename(columns={'index':'cellID'}).to_dict(orient='records'),'geneIntersection':features,'scriptSHA256':sha(__file__),'countsSHA256':sha(a.output/'counts.h5ad'),'excludedRegions':excluded})
    config={'format':'numivivo-molecular-spatial-assay/v1','id':'spatial-clu-mouse3','family':'molecular-perturbation','title':'Clu knockout · mouse brain','specimen':SOURCE_IDS[2],'intervention':{'kind':'gene-knockout','target':'Clu','barcode':TARGET,'control':CONTROL},'timepoints':['study-endpoint'],'regions':regions,'excludedRegions':excluded,'trainingRows':[i for i,g in enumerate(groups) if g['donorID'] in SOURCE_IDS[:2]],'trainingUnits':SOURCE_IDS[:2], 'sourceCitation':'https://doi.org/10.1038/s41467-026-69677-6','files':{n:sha(a.output/n) for n in ['source/report.json','source/receipt.json','specimen.json','design-specimen.json','protocol.json','lineage.json']},'runtime':{'binary':str(a.binary.resolve()),'tissueBinary':str(a.tissue_binary.resolve())},'limits':protocol['limits']}
    write(a.output/'assay.json',config)
    print(json.dumps({'assay':str(a.output/'assay.json'),'regions':regions,'excluded':excluded,'cells':data.shape[0],'genes':len(features)}))
if __name__=='__main__':main()
