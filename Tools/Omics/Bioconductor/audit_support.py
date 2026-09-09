#!/usr/bin/env python3
"""Describe native support rejections; this does not authorize reduced-model inference."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
samples=pd.read_csv(a.input/'samples.tsv',sep='\t');counts=pd.read_csv(a.input/'counts.tsv',sep='\t',index_col=0)
native=pd.read_csv(a.input/'native.tsv',sep='\t').set_index('featureID')
assert counts.columns.tolist()==samples.sampleID.tolist()
donors=samples.donor.astype(str).to_numpy();conditions=samples.condition.astype(str).to_numpy();levels=sorted(set(conditions));assert len(levels)==2
for donor in set(donors):assert sorted(conditions[donors==donor])==levels
histogram=Counter();rows=[]
for gene in native.index[native.status=='rankDeficientSupport']:
 y=counts.loc[gene].to_numpy();active=sorted(d for d in set(donors) if y[donors==d].sum()>0)
 use=np.isin(donors,active);local_donors=donors[use];local_y=y[use]
 x=np.column_stack([np.ones(use.sum()),(conditions[use]==levels[1]).astype(float)]+[(local_donors==d).astype(float) for d in active[1:]])
 positive=x[local_y>0];rank=int(np.linalg.matrix_rank(positive)) if len(positive) else 0
 full=rank==x.shape[1]
 histogram[(len(active),full)]+=1
 rows.append(dict(featureID=gene,activeDonors=len(active),zeroTotalDonors=len(set(donors))-len(active),reducedPositiveRank=rank,reducedColumns=x.shape[1],fullReducedPositiveRank=full))
result=dict(status='descriptive-support-audit',rejectedGenes=len(rows),categories=[dict(activeDonors=n,fullReducedPositiveRank=full,genes=count) for (n,full),count in sorted(histogram.items())],genes=rows,qualification='Algebraic support classification only. Donor removal, boundary nuisance treatment, dispersion profiling, residual degrees of freedom and uncertainty have not been fitted or qualified. This does not bypass the native inference gate.')
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='genes'},indent=2))
