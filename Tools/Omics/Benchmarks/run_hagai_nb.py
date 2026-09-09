#!/usr/bin/env python3
"""Full original Hagai mouse LPS6 paired comparison; no cell/gene subsetting.

Requires the pinned preparation from prepare_hagai_stream.py and original
SDRF/supplement. Source prefixes identify the paired individuals; Table 2
establishes three mouse individuals but does not supply a prefix-to-row map.
"""
import argparse,csv,hashlib,json,platform,subprocess,time,warnings
from importlib.metadata import version
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse,stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

EXPECTED={'Nfkb2':'ENSMUSG00000025225','Nfkbia':'ENSMUSG00000021025','Cxcl10':'ENSMUSG00000034855','Isg15':'ENSMUSG00000035692','Tnf':'ENSMUSG00000024401'}
SUPPLEMENT_SHA='dde588efff9da9febb2d882713759ccb71c75642ae4d3fb83902cc23c6787974'
SDRF_SHA='79e2a16fbe5d3248343bcb493c03f787a8a740dbedc546a72c1f98adf2b9fd63'
p=argparse.ArgumentParser(description=__doc__)
for name in ['binary','prepared','sdrf','supplement','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
def save(name,value):
    path=a.out/name;path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');return path
def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
assert sha(a.supplement)==SUPPLEMENT_SHA
assert sha(a.sdrf)==SDRF_SHA
preparation=json.loads((a.prepared/'preparation.json').read_text())
source=a.prepared/'prepared.h5ad';assert sha(source)==preparation['h5adSHA256']
assert (preparation['cells'],preparation['features'],preparation['nonzeros'])==(13863,22048,32848185)
with a.sdrf.open() as f:rows=list(csv.reader(f,delimiter='\t'))
header=rows[0];by_name={row[0].lower():row for row in rows[1:]}
mapping=json.loads((a.prepared/'plan.json').read_text())['mapping'];identities=[]
for sample in mapping['samples']:
    sid=sample['id'];row=by_name[sid]
    for key,value in [('Characteristics[organism]','Mus musculus'),('Characteristics[strain]','C57BL/6'),('Characteristics[sex]','female'),('Characteristics[age]','8 weeks')]:assert row[header.index(key)]==value
    donor=sid.split('_')[0];assert donor in ['mouse1','mouse2','mouse3']
    sample.update(biologicalReplicateID=donor,donorID=donor)
    identities.append(dict(sampleID=sid,donorID=donor,sourceName=row[0],enaSample=row[header.index('Comment[ENA_SAMPLE]')]))
mapping['sourceDescription']='E-MTAB-6754 original author QC cluster0; mouse1/2/3 paired time courses, supported by SDRF naming and Supplementary Table 2 three individuals'
contrast=dict(id='LPS6-vs-unstimulated',controlCondition='unstimulated',treatmentCondition='LPS6',design='pairedDonors',sizeFactors='medianRatio',variance='empiricalBayes',adjustForBatch=False,minimumCellsPerPseudobulk=10,minimumReplicatesPerCondition=3,minimumFeatureCounts=10,minimumExpressingPseudobulks=3,minimumReferenceFeatures=10,priorCount=0.5,model='negativeBinomial',negativeBinomialOptions=dict(trend='parametric'))
plan=save('plan.json',dict(schemaVersion=1,mapping=mapping,contrasts=[contrast]))
save('design-and-acceptance.json',dict(scope='All six author-filtered mouse source files, LPS6 versus unstimulated',sourcePreparation=preparation,
    donorEvidence=dict(supplementURL='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6347972/supplementaryFiles',supplementSHA256=SUPPLEMENT_SHA,pdfPage=12,table=2,sdrfSHA256=sha(a.sdrf),samples=identities,
    interpretation='Three mouse individuals explicitly listed; pairing inferred from deposited mouse1/2/3 time-course prefixes. No claimed mapping to table row numbers 4/5/6; no batch assignment inferred.'),
    preregisteredBeforeFit=dict(expectedPositiveGenes=EXPECTED,minimumPositiveInBoth=4,minimumEffectSpearman=0.95,maximumNativeNumericalFailures=0),
    boundary='A six-hour contrast, not the four-hour paper DE reproduction; two residual degrees of freedom and one study do not establish production calibration'))
commands=[]
def run(label,*args):
    command=[str(a.binary.resolve()),*map(str,args)];start=time.perf_counter()
    result=subprocess.run(command,capture_output=True,text=True)
    (a.out/(label+'.log')).write_text(result.stdout+result.stderr)
    commands.append(dict(command=command,returncode=result.returncode,seconds=time.perf_counter()-start));save('commands.json',commands)
    assert result.returncode==0,(label,result.stderr)
bundle=a.out/'native';run('native','singlecell-h5ad-pseudobulk',source,'--plan',plan,'--output',bundle)
run('verify','singlecell-h5ad-pseudobulk-verify',bundle)
native=json.loads((bundle/'report.json').read_text());nb=native['pseudobulk'];matrix=nb['matrix']
bulk=sparse.csr_matrix((np.array(matrix['counts'],dtype=np.int64),matrix['featureIndices'],matrix['rowOffsets']),shape=(6,22048)).toarray()
reference=np.load(a.prepared/'reference.npz');genes=reference['genes'].tolist();assert nb['featureIDs']==genes
for i,g in enumerate(nb['groups']):
    assert len(g['sampleIDs'])==1
    k=reference['sample_ids'].tolist().index(g['sampleIDs'][0]);np.testing.assert_array_equal(bulk[i],reference['bulk'][k])
np.testing.assert_array_equal([q['totalCounts'] for q in native['quality']],reference['totals'])
np.testing.assert_array_equal([q['detectedFeatures'] for q in native['quality']],reference['detected'])
assert native['canonicalNonzeros']==32848185 and len(native['metadata']['cells'])==13863
groups=nb['groups'];metadata=pd.DataFrame(dict(donor=[g['donorID'] for g in groups],condition=[g['condition'] for g in groups]),index=[g['sampleIDs'][0] for g in groups]).astype('category')
assert metadata.groupby('donor',observed=True).size().tolist()==[2,2,2]
eligible=(bulk.sum(axis=0)>=10)&((bulk>0).sum(axis=0)>=3)
counts=pd.DataFrame(bulk[:,eligible],index=metadata.index,columns=np.array(genes)[eligible])
counts.to_csv(a.out/'reference-pseudobulk.tsv',sep='\t');metadata.to_csv(a.out/'reference-design.tsv',sep='\t')
start=time.perf_counter()
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    dds=DeseqDataSet(counts=counts,metadata=metadata,design='~ donor + condition',refit_cooks=False,n_cpus=2,quiet=True)
    dds.deseq2();ds=DeseqStats(dds,contrast=['condition','LPS6','unstimulated'],cooks_filter=False,independent_filter=False,n_cpus=2,quiet=True);ds.summary()
reference_warnings=[dict(category=w.category.__name__,message=str(w.message)) for w in caught]
dds.var.to_csv(a.out/'reference-gene-diagnostics.tsv',sep='\t')
ref=ds.results_df;ref.to_csv(a.out/'pydeseq2.tsv',sep='\t')
de=native['contrasts'][0];table=pd.DataFrame(de['features']).set_index('featureID');table.to_csv(a.out/'native-de.tsv',sep='\t')
joined=table.join(ref,rsuffix='_reference');finite=joined[(joined.status=='tested')&np.isfinite(joined.log2FoldChange)&np.isfinite(joined.log2FoldChange_reference)]
def number(x):return float(x) if pd.notna(x) and np.isfinite(x) else None
biology={name:dict(featureID=gene,nativeStatus=str(table.loc[gene,'status']),nativeLog2Effect=number(table.loc[gene,'log2FoldChange']),referenceLog2Effect=number(ref.loc[gene,'log2FoldChange']),nativeBH=number(table.loc[gene,'adjustedPValue']),referenceBH=number(ref.loc[gene,'padj'])) for name,gene in EXPECTED.items()}
positive=sum(x['nativeStatus']=='tested' and x['nativeLog2Effect'] is not None and x['referenceLog2Effect'] is not None and x['nativeLog2Effect']>0 and x['referenceLog2Effect']>0 for x in biology.values())
correlation=float(stats.spearmanr(finite.log2FoldChange,finite.log2FoldChange_reference).statistic)
failures=int((table.status=='numericalFailure').sum())
passed=positive>=4 and correlation>=0.95 and failures==0
report=dict(status='passed-predeclared-descriptive-comparison' if passed else 'failed-predeclared-comparison',cells=13863,genes=22048,nonzeros=32848185,donors=3,pseudobulks=6,residualDF=de['design']['residualDegreesOfFreedom'],
    binarySHA256=sha(a.binary),sourceSHA256=sha(source),versions={k:version(k) for k in ['numpy','pandas','scipy','pydeseq2']},platform=platform.platform(),exactSourceQCAndPseudobulks=True,
    nativeStatuses={str(k):int(v) for k,v in table.status.value_counts().items()},eligibleGenes=int(eligible.sum()),commonFiniteGenes=len(finite),effectSpearman=correlation,effectSignAgreement=float(np.mean(np.sign(finite.log2FoldChange)==np.sign(finite.log2FoldChange_reference))),top50BHOverlap=len(set(finite.nsmallest(50,'adjustedPValue').index)&set(finite.nsmallest(50,'padj').index)),
    nativeTrend=de['negativeBinomial']['trend'],nativeNumericalFailures=failures,referenceWarnings=reference_warnings,referenceSeconds=time.perf_counter()-start,
    expectedBiology=dict(positiveInBoth=positive,required=4,genes=biology),qualification='Predeclared descriptive real-study comparison; no causal, calibrated-FDR, multi-study competitiveness or clinical claim; source-prefix donor pairing inference explicit in design')
save('report.json',report);print(json.dumps(report,indent=2));raise SystemExit(0 if passed else 1)
