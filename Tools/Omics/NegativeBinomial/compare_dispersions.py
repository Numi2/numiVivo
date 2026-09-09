#!/usr/bin/env python3
"""Describe native/PyDESeq2 dispersion differences without equating estimators."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__)
for name in ['native-report','reference-diagnostics','out']:p.add_argument('--'+name,type=Path,required=True)
a=p.parse_args();r=json.loads(a.native_report.read_text());c=r['contrasts'][0]
ref=pd.read_csv(a.reference_diagnostics,sep='\t',index_col=0)
features=c['features'];diags=c['negativeBinomial']['features'];rows=[]
for f,d in zip(features,diags):
    if f['status']!='tested' or f['featureID'] not in ref.index:continue
    other=ref.loc[f['featureID']]
    rows.append(dict(featureID=f['featureID'],nativeFinal=d['finalDispersion'],referenceFinal=float(other['dispersions']),nativeTrend=d['trendDispersion'],referenceTrend=float(other['fitted_dispersions'])))
t=pd.DataFrame(rows).set_index('featureID');assert len(t)>0 and np.isfinite(t.to_numpy()).all() and (t.to_numpy()>0).all()
ratios=t.nativeFinal/t.referenceFinal
markers=['ENSMUSG00000025225','ENSMUSG00000021025','ENSMUSG00000034855','ENSMUSG00000035692','ENSMUSG00000024401']
result=dict(status='descriptive-estimator-differences',commonTestedGenes=len(t),nativeOverReferenceFinalDispersionQuantiles={str(q):float(ratios.quantile(q)) for q in [0,.1,.5,.9,1]},
    responseGenes={g:t.loc[g].to_dict() for g in markers if g in t.index},
    referenceNonconverged={col:int((~ref[col]).sum()) for col in ['_genewise_converged','_MAP_converged','_LFC_converged']},
    sourceSHA256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [a.native_report,a.reference_diagnostics]},
    boundary='Differences in fitted dispersions explain a source of Wald uncertainty differences; neither agreement nor disagreement establishes calibration. No estimator retuned to these results.')
with a.out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
