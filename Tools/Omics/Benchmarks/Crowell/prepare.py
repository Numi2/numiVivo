#!/usr/bin/env python3
"""Create a complete deposited-count AnnData view and freeze eight animal contrasts."""
import argparse,gzip,hashlib,json,zipfile,xml.etree.ElementTree as ET
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(r/'source/Crowell19_4vs4.Rda')=='3193a33fa7acb650fdb1e682822705d7c68057aba93d84dac63a364e40fdb1fe'
m=json.loads(gzip.decompress((r/'wire/metadata.json.gz').read_bytes()))
arrays=[np.frombuffer(gzip.decompress((r/'wire'/name).read_bytes()),dtype='<i4').copy() for name in ['values.i32.gz','indices.i32.gz','indptr.i32.gz']]
x=csr_matrix(tuple(arrays),shape=tuple(reversed(m['dim'])));assert x.has_canonical_format and x.nnz==m['nonzeros'] and int(x.sum())==m['umiTotal'] and (x.data>0).all()
def frame(columns,index):
 d=pd.DataFrame(index=pd.Index(index,dtype=object))
 for name,c in columns.items():
  if c['type']=='factor':v=pd.Categorical(c['values'],categories=c['levels'] if isinstance(c['levels'],list) else [c['levels']],ordered=c['ordered'])
  else:v=pd.Series(c['values'],index=d.index,dtype=object if c['type']=='character' else None)
  d[name]=v
 return d
obs=frame(m['obs'],m['obsNames']);var=frame(m['var'],m['varNames']);assert var.ENSEMBL.is_unique and obs.index.is_unique
# Read deposited spreadsheet cell strings without interpreting formatted dates.
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
with zipfile.ZipFile(r/'source/metadata.xlsx') as z:
 strings=[''.join(e.itertext()) for e in ET.fromstring(z.read('xl/sharedStrings.xml'))]
 cells={}
 for c in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('.//m:c',ns):
  v=c.find('m:v',ns)
  if v is not None:cells[c.attrib['r']]=strings[int(v.text)] if c.attrib.get('t')=='s' else v.text
 samples={cells[f'D{i}']:('Vehicle' if cells[f'F{i}']=='WT' else cells[f'F{i}']) for i in range(2,10)}
assert set(samples)==set(obs.sample_id.astype(str)) and all(obs.group_id.astype(str)==obs.sample_id.astype(str).map(samples))
assert list(samples.values()).count('Vehicle')==list(samples.values()).count('LPS')==4
obj=ad.AnnData(x,obs=obs,var=var);obj.uns['numivivo_source']=dict(sourceRdaSHA256=sha(r/'source/Crowell19_4vs4.Rda'),source='muscData/Crowell19_4vs4.Rda',assaySelection=m['assaySelection'],qualification='Complete deposited post-QC count assay and row/column metadata; normalized assay and reductions remain in the original Rda. This is not a complete SCE conversion.')
out=r/'counts.h5ad';assert not out.exists();obj.write_h5ad(out,compression='gzip',compression_opts=6,convert_strings_to_categoricals=False)
restored=ad.read_h5ad(out);assert (restored.X!=x).nnz==0
# Pandas 3 reconstructs character axes/columns with StringDtype. Compare the
# source R character semantics after explicit dtype normalization; factor levels,
# ordering, numerical and logical columns still require exact dtype/value equality.
def normalize_characters(df,columns):
 df=df.copy();df.index=pd.Index(df.index,dtype=object);df.columns=pd.Index(df.columns,dtype=object)
 for name,c in columns.items():
  if c['type']=='character':df[name]=df[name].astype(object)
 return df
pd.testing.assert_frame_equal(normalize_characters(restored.obs,m['obs']),normalize_characters(obs,m['obs']))
pd.testing.assert_frame_equal(normalize_characters(restored.var,m['var']),normalize_characters(var,m['var']))
# Independent Scanpy QC; no additional cell/feature selection.
qc=obj.copy();qc.var['mitochondrial']=qc.var.SYMBOL.str.lower().str.startswith('mt-').to_numpy();sc.pp.calculate_qc_metrics(qc,qc_vars=['mitochondrial'],percent_top=None,log1p=False,inplace=True)
qc.obs[['total_counts','n_genes_by_counts','total_counts_mitochondrial','pct_counts_mitochondrial']].to_csv(r/'scanpy-qc.tsv.gz',sep='\t',index_label='sourceCellID',compression={'method':'gzip','mtime':0})
map_plan=dict(schemaVersion=1,id='crowell2019-cortex',evidence='measured',countUnit='umiCount',matrixPath='X',sampleColumn='sample_id',barcodeColumn='barcode',groupColumn='cluster_id',featureIDColumn='ENSEMBL',featureNameColumn='SYMBOL',mitochondrialFeatureIDs=var.loc[qc.var.mitochondrial,'ENSEMBL'].tolist(),sourceDescription='Crowell et al. 2020, doi:10.1038/s41467-020-19894-4; complete author-deposited post-QC counts from muscData Crowell19_4vs4. Four independent vehicle and four independent LPS mice; snRNA-seq cortex. Source cell labels and upstream QC are author annotations; logcounts and reductions are retained only in the original Rda.',samples=[dict(id=s,biologicalReplicateID=s,donorID=s,condition=samples[s],batchID='unreported',organism='NCBITaxon:10090') for s in sorted(samples)])
cases=[]
for i,group in enumerate(obs.cluster_id.cat.categories):
 name=f'{i:02d}';d=r/name;d.mkdir(exist_ok=False)
 c=dict(id=f'LPS-vs-Vehicle-{i:02d}',model='negativeBinomial',design='independentReplicates',controlCondition='Vehicle',treatmentCondition='LPS',cellGroup=group,minimumCellsPerPseudobulk=10,minimumReplicatesPerCondition=3,minimumFeatureCounts=10,minimumExpressingPseudobulks=3,minimumReferenceFeatures=10,sizeFactors='medianRatio',adjustForBatch=False,negativeBinomialOptions=dict(trend='gammaParametric',effectPriorStandardDeviationLog2=1))
 plan=dict(schemaVersion=1,mapping=map_plan,contrasts=[c]);(d/'plan.json').write_text(json.dumps(plan,sort_keys=True)+'\n')
 sizes=obs.loc[obs.cluster_id==group].groupby('sample_id',observed=False).size();valid=[s for s,n in sizes.items() if n>=10];rep={condition:sum(samples[s]==condition for s in valid) for condition in ['Vehicle','LPS']}
 cases.append(dict(id=name,cellGroup=group,cells=int(sizes.sum()),sampleCellCounts={str(k):int(v) for k,v in sizes.items()},replicatesAfterCellGate=rep,expectedReplicationGate='ready' if min(rep.values())>=3 else 'insufficient',planSHA256=sha(d/'plan.json')))
protocol=dict(sourceRdaSHA256=sha(r/'source/Crowell19_4vs4.Rda'),sourceH5ADSHA256=sha(out),sourceH5ADBytes=out.stat().st_size,cells=obj.n_obs,features=obj.n_vars,nonzeros=x.nnz,umiTotal=int(x.sum()),countsRoundTripExact=True,allRowAndColumnMetadataExact=True,cases=cases,expectedBiology=dict(source='https://pmc.ncbi.nlm.nih.gov/articles/PMC7705760/',hypothesis='LPS-associated immune response across endothelial and glial populations; source analysis reports stronger effects outside neurons.',predeclaredGeneSymbols=['Stat1','Irf7','Isg15','B2m','H2-D1'],panelRole='Descriptive positive-direction hypotheses; all missing, unavailable and opposing effects retained. Neither author labels nor reference calls are independent biological truth.'),scope='All 25224 deposited post-QC nuclei and 11076 deposited genes; no new gene/cell subset. CPE replication failures retained; no pairing or batch inferred. Only pseudobulk is dense.',software=dict(anndata=ad.__version__,scanpy=sc.__version__,numpy=np.__version__,pandas=pd.__version__))
(r/'protocol.json').write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n');print(json.dumps(protocol,indent=2))
