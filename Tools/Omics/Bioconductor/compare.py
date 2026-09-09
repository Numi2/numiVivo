#!/usr/bin/env python3
"""Describe method disagreements without using any method as biological ground truth."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
meta=json.loads((a.input/'input.json').read_text());native=pd.read_csv(a.input/'native.tsv',sep='\t',float_precision='round_trip').set_index('featureID')
eligible=native.index[native.eligible].tolist();runs=json.loads((a.reference/'runs.json').read_text());results={}
def check_bh(t):
 mask=np.isfinite(t.pValue);values=t.loc[mask,'pValue'].to_numpy();assert ((values>=0)&(values<=1)).all()
 np.testing.assert_allclose(multipletests(values,method='fdr_bh')[1],t.loc[mask,'adjustedPValue'],rtol=2e-13,atol=1e-15)
check_bh(native)
for name,run in runs.items():
 assert run['status']=='completed',run
 path=a.reference/(name+'.tsv');ref=pd.read_csv(path,sep='\t',float_precision='round_trip').set_index('featureID');assert ref.index.tolist()==eligible
 check_bh(ref)
 joined=native.join(ref,rsuffix='_reference');use=(joined.status=='tested')&np.isfinite(joined.log2FoldChange)&np.isfinite(joined.log2FoldChange_reference)
 j=joined.loc[use];delta=j.log2FoldChange_reference-j.log2FoldChange
 common_p=j.loc[np.isfinite(j.adjustedPValue)&np.isfinite(j.adjustedPValue_reference)]
 k=min(100,len(common_p));ntop=set(common_p.nsmallest(k,'adjustedPValue').index);rtop=set(common_p.nsmallest(k,'adjustedPValue_reference').index)
 expected={}
 for gene in meta['expectedGenes']:
  if gene not in j.index:expected[gene]=dict(status='not jointly testable',nativeStatus=str(native.loc[gene,'status']) if gene in native.index else 'absent');continue
  row=j.loc[gene];expected[gene]=dict(status='jointly testable',nativeEffect=float(row.log2FoldChange),referenceEffect=float(row.log2FoldChange_reference),nativeBH=float(row.adjustedPValue),referenceBH=float(row.adjustedPValue_reference))
 joint_native_bh=multipletests(common_p.pValue,method='fdr_bh')[1]
 joint_reference_bh=multipletests(common_p.pValue_reference,method='fdr_bh')[1]
 unavailable=ref.index.difference(j.index)
 largest=delta.abs().nlargest(20).index
 results[name]=dict(referenceFeatures=len(ref),nativeMultiplicityFeatures=int(native.pValue.notna().sum()),referenceMultiplicityFeatures=int(ref.pValue.notna().sum()),jointUniverseRecomputedBH=dict(features=len(common_p),nativeBelow005=int((joint_native_bh<.05).sum()),referenceBelow005=int((joint_reference_bh<.05).sum()),qualification='Descriptive sensitivity on the intersection, not the native or reference published multiplicity family'),referenceSignificantAmongNativeUnavailable=int((ref.loc[unavailable,'adjustedPValue']<.05).sum()),nativeStatusCounts=native.status.value_counts().to_dict(),jointlyTestedEffects=len(j),referenceMissingPValues=int(ref.pValue.isna().sum()),referenceMissingEffects=int(ref.log2FoldChange.isna().sum()),
  effectSpearman=float(spearmanr(j.log2FoldChange,j.log2FoldChange_reference).statistic),effectSignAgreement=float(np.mean(np.sign(j.log2FoldChange)==np.sign(j.log2FoldChange_reference))),medianAbsoluteEffectDifference=float(np.median(np.abs(delta))),
  nativeBHBelow005=int((native.adjustedPValue<.05).sum()),referenceBHBelow005=int((ref.adjustedPValue<.05).sum()),nativeZeroPValues=int((native.pValue==0).sum()),referenceZeroPValues=int((ref.pValue==0).sum()),topK=k,topKBHOverlap=len(ntop&rtop),expectedGenes=expected,
  largestEffectDisagreements=[dict(featureID=g,nativeEffect=float(j.loc[g,'log2FoldChange']),referenceEffect=float(j.loc[g,'log2FoldChange_reference']),nativeMeanNormalizedCount=float(j.loc[g,'meanNormalizedCount'])) for g in largest],
  referenceBetaNonconvergence=int((ref.betaConverged==False).sum()) if 'betaConverged' in ref else None,referenceDispersionOutliers=int(ref.dispersionOutlier.sum()) if 'dispersionOutlier' in ref else None,referenceSHA256=hashlib.sha256(path.read_bytes()).hexdigest())
 if name.endswith('DESeq2'):
  default=pd.read_csv(a.reference/(name+'-default-results.tsv'),sep='\t',float_precision='round_trip')
  results[name]['DESeq2DefaultResultPolicy']=dict(missingPValues=int(default.pvalue.isna().sum()),missingAdjustedPValues=int(default.padj.isna().sum()),BHBelow005=int((default.padj<.05).sum()))
out=dict(status='completed-descriptive-comparison',input=meta,results=results,independentBHVerified=True,qualification='Same real pseudobulk observations, prefilter and paired donor design. Native support rejections reduce its tested multiplicity family; raw and intersection-recomputed BH summaries are distinct. Different estimators and tests are retained, not forced into numerical equality. No calibrated FDR, causal truth or general competitive claim.')
a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:{f:v[f] for f in ['effectSpearman','effectSignAgreement','nativeBHBelow005','referenceBHBelow005','topKBHOverlap','referenceBetaNonconvergence']} for k,v in results.items()},indent=2))
