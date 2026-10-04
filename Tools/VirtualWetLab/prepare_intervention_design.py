#!/usr/bin/env python3
"""Bounded source admission and target-aware inputs for the existing native learner.

The reserved response rows are not accessed by prepare(). reveal_design.py is
the only stage that summarizes them. Accession/batch is never inferred from RNA.
"""
import argparse, gzip, hashlib, io, json, re, tarfile
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import requests
from scipy import sparse
from scipy.io import mmread
from safetensors.numpy import save_file
from wetlab import read, write, sha, require, timestamp

SEED=271828
V03=Path('/Users/home/virtual-wet-lab-v03-20261003')


def orthology(path):
    table=pd.read_csv(path,sep='\t');mapping={}
    for _,g in table.groupby('DB Class Key'):
        m=g[g['NCBI Taxon ID']==10090];h=g[g['NCBI Taxon ID']==9606]
        if len(m)==len(h)==1:mapping[str(m.Symbol.iloc[0])]=str(h.Symbol.iloc[0])
    return mapping


def mouse(path):
    with tarfile.open(path) as t:
        genes=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile('features.tsv.gz').read())),sep='\t',header=None)
        bars=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile('barcodes.tsv.gz').read())),header=None)[0].astype(str).tolist()
        x=mmread(gzip.GzipFile(fileobj=t.extractfile('matrix.mtx.gz'))).T.tocsr()
    guide=np.flatnonzero(genes[1].str.endswith('_gRNA'));g=x[:,guide].toarray();n=(g>0).sum(1)
    targets=np.asarray(['barcode-negative']*len(bars),object);targets[n>1]='multiple-guides'
    names=genes.loc[guide,1].str.replace('_gRNA','',regex=False).to_numpy()
    targets[n==1]=names[g[n==1].argmax(1)];targets[targets=='mSafe']='control'
    targets[np.isin(targets,['Myrf','Ndufaf','Rbfox'])]='ambiguous-source-alias'
    keep=genes[0].str.startswith('ENSMUSG')&~genes[1].duplicated(keep=False)
    x=x[:,keep.to_numpy()];symbols=genes.loc[keep,1].tolist()
    obs=pd.DataFrame({'target':targets,'original_barcode':bars},index=bars)
    return ad.AnnData(x,obs=obs,var=pd.DataFrame(index=symbols))


def source_metadata(root):
    if (root/'sources-manifest.json').exists():return read(root/'sources-manifest.json')
    specs=[]
    for path in sorted((root/'sources').glob('GSM844277*.tar.gz')):
        tag=path.name.split('_filtered')[0];batch=re.search(r'B([123])_',tag)[1]
        specs.append({'id':tag,'path':str(path),'kind':'mouse','split':{'1':'training','2':'validation','3':'reserved'}[batch],
          'study':'GSE274058','unit':'GSE274058:batch'+batch,'context':'mouse-hippocampus-mixed','modality':'CRISPR-KO','species':'mouse',
          'accession':path.name.split('_')[0],'batch':batch,'independence':'batch blocks; five named mice versus three pooled experiments conflict retained'})
    for name,study,split,context,modality in [
      ('AdamsonWeissman2016_GSM2406675_10X001.h5ad','GSE90546','training','K562','CRISPRi'),
      ('DixitRegev2016.h5ad','GSE90063','training','K562-high-MOI','CRISPR-KO'),
      ('DatlingerBock2017.h5ad','GSE92872','reserved','Jurkat','CRISPR-KO')]:
        specs.append({'id':study,'path':str(root/'sources'/name),'kind':'harmonized','split':split,'study':study,'unit':study,
          'context':context,'modality':modality,'species':'human','accession':study,
          'independence':'one source study; capture runs are not independent biological replicates'})
    return specs


def load_source(s):
    if s['kind']=='mouse':return mouse(Path(s['path']))
    a=ad.read_h5ad(s['path'])
    if s['study']=='GSE90546':
        labels=a.obs.perturbation.astype(str);a.obs['target']=labels.str.split('_').str[0]
        a.obs.loc[labels.str.contains('(mod)',regex=False),'target']='control'
        a.obs.loc[labels=='*','target']='unassigned'
    elif s['study']=='GSE90063':
        a.obs['target']=a.obs.target.astype(str)
        a.obs.loc[a.obs.target.str.startswith('INTERGENIC') & ~a.obs.target.str.contains(' + ',regex=False),'target']='control'
        a.obs.loc[(a.obs.nperts!=1)|a.obs.target.str.contains(' + ',regex=False),'target']='multiple-guides'
    elif s['study']=='GSE168620':
        a.obs['target']=a.obs.perturbation.astype(str).str.split('_').str[0]
    else:
        a.obs['target']=a.obs.target.astype(str)
        a.obs.loc[a.obs.perturbation.astype(str)=='control','target']='control'
    if s.get('captureFilter'):
        capture=s['captureFilter'];mask=a.obs.replicate.astype(str).eq(str(capture['replicate']))
        a=a[mask if capture['include'] else ~mask].copy()
    return a


def contexts(s,a):
    if s['study'] in ('GSE92872','GSE168620'):return a.obs.perturbation_2.astype(str).fillna('unassigned').to_numpy(dtype=str)
    return np.asarray([s['context']]*a.n_obs)


def axis(a,mapping,species):
    names=[mapping.get(str(g)) if species=='mouse' else str(g) for g in a.var_names]
    counts=pd.Series(names).value_counts();return {g:i for i,g in enumerate(names) if g is not None and counts[g]==1}


def normalized(a,rows,cols):
    x=a.X[rows].tocsr();total=np.maximum(1,np.asarray(x.sum(1)).ravel())
    return np.log1p(x[:,cols].toarray()/total[:,None]*1e6).astype(np.float32)


def fetch_descriptors(targets,output):
    output.mkdir(exist_ok=True);records={};aa='ACDEFGHIKLMNPQRSTVWY'
    # EBI serves the same UniProt records; its API is usable when the primary
    # endpoint returns 503. MGI supplies pinned primary accession identities.
    from concurrent.futures import ThreadPoolExecutor
    mgi=pd.read_csv(V03/'research/MGI-HOM_MouseHumanSequence.rpt',sep='\t')
    accession={str(r.Symbol):str(r['SWISS_PROT IDs']).split(',')[0] for _,r in mgi[mgi['NCBI Taxon ID']==9606].iterrows() if pd.notna(r['SWISS_PROT IDs'])}
    def fetch(gene):
        if gene not in accession:return gene,None
        acc=accession[gene];path=output/(acc+'.json')
        if not path.exists():
            r=requests.get('https://www.ebi.ac.uk/proteins/api/proteins/'+acc,headers={'Accept':'application/json'},timeout=60);r.raise_for_status();path.write_text(r.text)
        return gene,path
    with ThreadPoolExecutor(max_workers=4) as pool:
        for gene,path in pool.map(fetch,targets):
            if path is None:continue
            entry=read(path);require(entry['organism']['taxonomy']==9606,'Target species mismatch')
            seq=entry['sequence']['sequence'];go=sorted(x['id'] for x in entry.get('dbReferences',[]) if x['type']=='GO')
            candidate={'accession':entry['accession'],'sequenceSHA256':hashlib.sha256(seq.encode()).hexdigest(),'go':go,'sourceSHA256':sha(path),'sourceURL':'https://www.ebi.ac.uk/proteins/api/proteins/'+entry['accession'],'recordVersion':entry['info'],'annotationEvidence':'GO evidence retained in source JSON; not a causal response mechanism'}
            vec=np.zeros(149,np.float32);vec[:20]=[seq.count(c)/len(seq)*10 for c in aa];vec[20]=np.log1p(len(seq))/10
            for term in go:vec[21+int(hashlib.sha256(term.encode()).hexdigest()[:8],16)%128]+=1
            vec[21:]/=max(1,np.linalg.norm(vec[21:]));candidate['descriptor']=vec.tolist();records[gene]=candidate
    manifest=output/('targets-'+hashlib.sha256('\n'.join(targets).encode()).hexdigest()[:16]+'.json')
    if manifest.exists():require(read(manifest)==records,'Changed cached target records')
    else:write(manifest,records)
    return records


def prepare(root,out):
    require(not out.exists(),'Preparation already exists');out.mkdir()
    protocol=read(root/'preregistration.json');write(out/'protocol.json',protocol)
    mapping=orthology(V03/'research/MGI-HOM_MouseHumanSequence.rpt');specs=source_metadata(root)
    require(len(specs)==8 or (root/'sources-manifest.json').exists(),'All five mouse samples and three human studies required')
    hashes=[sha(s['path']) for s in specs]
    for digest in set(hashes):
        matches=[s for s,h in zip(specs,hashes) if h==digest]
        if len(matches)>1:
            require(len(matches)==2 and all(s.get('captureFilter') for s in matches),'Duplicate source content')
            a,b=[s['captureFilter'] for s in matches]
            require(a['replicate']==b['replicate'] and a['include']!=b['include'],'Overlapping source row filters')
    for s,h in zip(specs,hashes):s['sha256']=h
    write(out/'reservation.json',{'createdAt':timestamp(),'sources':specs,'pretrainedExpressionSources':[],
        'orthologySHA256':sha(V03/'research/MGI-HOM_MouseHumanSequence.rpt'),
        'priorExposures':protocol.get('priorExposures',['GSE274447 all spatial outcomes','GSE146194 v0.3 Arc local assay']),
        'reserved':[s['id']+' all treated RNA' for s in specs if s['split']=='reserved'],
        'protocolSHA256':sha(out/'protocol.json'),'testOutcomeUse':'not read by preparation; counts/guide labels only for admission'})
    axes=[];targetset=set();groups=[];exclusions=[];varsum={};all_fingerprints={};duplicate_rows=[]
    for s in specs:
        a=load_source(s);ax=axis(a,mapping,s['species']);axes.append(set(ax));labels=a.obs.target.astype(str).fillna('unassigned').to_numpy(dtype=str);ctx=contexts(s,a)
        s['cells']=a.n_obs;s['featureCount']=a.n_vars
        for context in sorted(set(ctx)):
            controls=np.flatnonzero((ctx==context)&(labels=='control'))
            for target in sorted(set(labels[ctx==context])-{'control'}):
                idx=np.flatnonzero((ctx==context)&(labels==target));canonical=mapping.get(target) if s['species']=='mouse' else target
                reasons=[]
                if target in ('barcode-negative','multiple-guides','ambiguous-source-alias','unassigned','nan') or not canonical:reasons.append('unassigned, ambiguous or unmapped target')
                if len(idx)<5:reasons.append('fewer than five assigned cells')
                if len(controls)<5:reasons.append('fewer than five explicit controls')
                row={'source':s['id'],'context':context,'target':canonical,'sourceTarget':target,'rows':idx.tolist(),'controls':controls.tolist(),'split':s['split'],'unit':s['unit'],'modality':s['modality'],'species':s['species'],'reasons':reasons}
                if reasons:exclusions.append({k:v for k,v in row.items() if k not in ('rows','controls')}|{'cells':len(idx),'controlCells':len(controls)})
                else:groups.append(row);targetset.add(canonical)
            if s['split']=='training' and len(controls):
                names=sorted(ax);v=normalized(a,controls,[ax[g] for g in names]).var(0)
                for g,value in zip(names,v):varsum[g]=varsum.get(g,0)+float(value)
        del a
        print(s['id'],'admitted metadata',flush=True)
    descriptors=fetch_descriptors(sorted(targetset),root/'target-sources')
    admitted=[]
    for g in groups:
        if g['target'] in descriptors:admitted.append(g)
        else:exclusions.append({k:v for k,v in g.items() if k not in ('rows','controls')}|{'reasons':['No reviewed source-bound target descriptor']})
    groups=admitted;common=set.intersection(*axes);ranked=sorted(common,key=lambda g:(-varsum.get(g,0),g));features=ranked[:512]
    for o in protocol['objectives']:
        for g in o['genes']:
            if g in common and g not in features:features.append(g)
    write(out/'features.json',features);write(out/'groups.json',groups);write(out/'exclusions.json',exclusions)
    targets=sorted({g['target'] for g in groups});write(out/'targets.json',{t:descriptors[t] for t in targets});tindex={t:i for i,t in enumerate(targets)}
    # Biological target block149; modality3; species2; context2. Unknown dose/time
    # are absent, not zero-dose or zero-hour treatment claims.
    arrays={split:{k:[] for k in ['context','descriptor','target','known','observed','observedVariance','mask','stratum']} for split in ('training','validation','reserved')};meta={k:[] for k in arrays};rng=np.random.default_rng(SEED)
    for s in specs:
        a=load_source(s);ax=axis(a,mapping,s['species']);cols=[ax[g] for g in features]
        for group in [g for g in groups if g['source']==s['id']]:
            split=s['split'];control=normalized(a,group['controls'],cols);base=control.mean(0)
            d=np.r_[descriptors[group['target']]['descriptor'],float(s['modality']=='CRISPR-KO'),float(s['modality']=='CRISPRi'),0.,float(s['species']=='mouse'),float(s['species']=='human'),float(group['context']=='stimulated'),float(group['context']=='unstimulated')].astype(np.float32)
            # Reserved treated rows are never normalized or summarized here.
            y=None if split=='reserved' else normalized(a,group['rows'],cols)
            chunks=[np.arange(len(group['rows']))] if split!='training' else [rng.choice(len(group['rows']),min(32,len(group['rows'])),replace=False) for _ in range(8)]
            for indices in chunks:
                ar=arrays[split];ar['context'].append(base);ar['descriptor'].append(d);ar['target'].append(tindex[group['target']]);ar['known'].append([0.]);ar['stratum'].append(len(meta[split])//8 if split=='training' else len(meta[split]));ar['mask'].append(np.ones(len(features),np.float32))
                if y is not None:ar['observed'].append(y[indices].mean(0));ar['observedVariance'].append(y[indices].var(0))
                meta[split].append({**group,'controlCellCount':len(control),'outcomeCellCount':len(group['rows']),'referenceState':'matched explicit control population; not a cell trajectory'})
        del a
    plan={'featureCount':len(features),'targetCount':len(targets),'descriptorCount':156,'hiddenWidth':64,'seed':SEED,'steps':[240,720,1440],'batchSize':16,'learningRate':.001,'weightDecay':.0001}
    write(out/'plan.json',plan)
    for split,ar in arrays.items():
        data={k:np.asarray(v,np.int32 if k in ('target','stratum') else np.float32) for k,v in ar.items() if v}
        if split=='training':data['prior']=np.zeros((len(targets),len(features)),np.float32)
        if split=='reserved':data={k:v for k,v in data.items() if k in ('context','descriptor','target','known')}
        require(data['descriptor'].shape[1]==156,'Descriptor axis');save_file(data,str(out/(split+'.safetensors')));write(out/(split+'-rows.json'),meta[split])
    write(out/'prepared.json',{'createdAt':timestamp(),'sources':specs,'features':len(features),'targets':len(targets),'groups':{k:len(v) for k,v in meta.items()},'exclusions':len(exclusions),'descriptorCount':156,'protocolSHA256':sha(out/'protocol.json'),'files':{p.name:sha(p) for p in out.iterdir() if p.is_file()},'overlapAudit':'Raw source hashes distinct; original studies distinct. No expression pretrained weights. Shared study mice conservatively split by batch. Spatial corpus excluded from these model fits.'})
    print('Prepared', {k:len(v) for k,v in meta.items()},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();prepare(a.root,a.output)
