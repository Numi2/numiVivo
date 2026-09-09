#!/usr/bin/env python3
"""Compare two native trend choices on identical experimental observations."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
p=argparse.ArgumentParser(description=__doc__)
for name in ['baseline','candidate','reference','out']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--expected',nargs='+',required=True)
a=p.parse_args();left=json.loads(a.baseline.read_text());right=json.loads(a.candidate.read_text())
for field in ['metadata','quality','pseudobulk','canonicalNonzeros']:assert left[field]==right[field],field
assert len(left['contrasts'])==len(right['contrasts'])==1
l=left['contrasts'][0];r=right['contrasts'][0];assert l['design']==r['design']
lr=json.loads(json.dumps(l['request']));rr=json.loads(json.dumps(r['request']))
del lr['negativeBinomialOptions']['trend'];del rr['negativeBinomialOptions']['trend'];assert lr==rr
reference=pd.read_csv(a.reference,sep='\t',index_col=0)
def describe(c):
    table=pd.DataFrame(c['features']).set_index('featureID');j=table.join(reference,rsuffix='_reference')
    finite=j[(j.status=='tested')&np.isfinite(j.log2FoldChange_reference)]
    expected={gene:dict(nativeStatus=str(j.loc[gene,'status']),nativeEffect=float(j.loc[gene,'log2FoldChange']),referenceEffect=float(j.loc[gene,'log2FoldChange_reference'])) for gene in a.expected}
    return dict(method=c['negativeBinomial']['trend']['method'],testedGenes=c['testedFeatures'],numericalFailures=int((table.status=='numericalFailure').sum()),
        effectSpearman=float(spearmanr(finite.log2FoldChange,finite.log2FoldChange_reference).statistic),effectSignAgreement=float(np.mean(np.sign(finite.log2FoldChange)==np.sign(finite.log2FoldChange_reference))),
        top50BHOverlap=len(set(finite.nsmallest(50,'adjustedPValue').index)&set(finite.nsmallest(50,'padj').index)),nativeBHBelow005=int((table.adjustedPValue<.05).sum()),
        nativeZeroPValues=int((table.pValue==0).sum()),expectedGenes=expected)
before=describe(l);after=describe(r)
assert after['numericalFailures']==0
assert all(v['nativeStatus']=='tested' and v['nativeEffect']>0 and v['referenceEffect']>0 for v in after['expectedGenes'].values())
result=dict(status='completed-controlled-trend-comparison',cells=len(right['metadata']['cells']),features=len(right['metadata']['features']),nonzeros=right['canonicalNonzeros'],identicalCountsQCMetadataAndDesign=True,baseline=before,candidate=after,
    inputs={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [a.baseline,a.candidate,a.reference]},
    boundary='Same observed data and experimental design; significance counts and reference agreement are descriptive, not FDR calibration or evidence for promoting either method.')
with a.out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
print(json.dumps({key:result[key] for key in ['status','cells','features','nonzeros']}))
