#!/usr/bin/env python3
"""Descriptive comparison of changed-coverage NB results with pinned R tables.

R fits keep all donor observations and use distinct estimation/testing methods.
No R method is treated as biological truth or native numerical reference.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--report',type=Path,required=True);p.add_argument('--reference',type=Path,required=True)
p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert not a.out.exists()
c=json.loads(a.report.read_text())['contrasts'][0]
native=pd.DataFrame(c['features']).set_index('featureID')
recovered=[f['featureID'] for f,d in zip(c['features'],c['negativeBinomial']['features']) if f['status']=='tested' and 'supportResolution' in d]
runs=json.loads((a.reference/'runs.json').read_text());out={}
for name,run in runs.items():
    assert run['status']=='completed'
    path=a.reference/(name+'.tsv');ref=pd.read_csv(path,sep='\t',float_precision='round_trip').set_index('featureID')
    j=native.loc[native.status=='tested'].join(ref,rsuffix='_reference')
    j=j.loc[j.pValue_reference.notna() & j.log2FoldChange_reference.notna()]
    def compare(t):
        if len(t)==0:return dict(features=0)
        qn=multipletests(t.pValue,method='fdr_bh')[1];qr=multipletests(t.pValue_reference,method='fdr_bh')[1]
        return dict(features=len(t),effectSpearman=float(spearmanr(t.log2FoldChange,t.log2FoldChange_reference).statistic),
            effectSignAgreement=float(np.mean(np.sign(t.log2FoldChange)==np.sign(t.log2FoldChange_reference))),
            medianAbsoluteEffectDifference=float(np.median(abs(t.log2FoldChange-t.log2FoldChange_reference))),
            nativePublishedBHBelow005=int((t.adjustedPValue<.05).sum()),referencePublishedBHBelow005=int((t.adjustedPValue_reference<.05).sum()),
            subsetRecomputedBHBelow005=dict(native=int((qn<.05).sum()),reference=int((qr<.05).sum())))
    out[name]=dict(joint=compare(j),recovered=compare(j.loc[j.index.intersection(recovered)]),
        referenceSignificantAmongNativeUnavailable=int((ref.loc[ref.index.difference(j.index),'adjustedPValue']<.05).sum()),
        nativePublishedFamily=int(native.pValue.notna().sum()),referencePublishedFamily=int(ref.pValue.notna().sum()),
        referenceSHA256=hashlib.sha256(path.read_bytes()).hexdigest())
summary=dict(status='completed-descriptive-R-comparison',results=out,
    sourceReportSHA256=hashlib.sha256(a.report.read_bytes()).hexdigest(),
    qualification='Distinct active-donor native profiles and full-donor R methods. Published-family BH and descriptive subset-recomputed BH are separate; no calibrated FDR or biological truth claim.')
a.out.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))
