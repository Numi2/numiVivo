#!/usr/bin/env python3
"""Post-score diagnostic of every native sham hit; never changes fits or filters."""
import argparse,gzip,hashlib,json,subprocess
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--results',type=Path,required=True);p.add_argument('--source-root',type=Path,default=Path('/Users/home'));p.add_argument('--out',type=Path,required=True)
a=p.parse_args();assert not a.out.exists();hits=[]
for study in ['kang','hagai']:
 source=a.source_root/('numivivo-'+study+'-r-native-current-20260909')/'original.h5ad'
 prep=json.loads((a.root/study/'preparation.json').read_text())
 with source.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==prep['sourceSHA256']
 obj=ad.read_h5ad(source);x=obj.X.tocsc();base=json.loads((a.root/study/'1/memberships.json').read_text());ids=base['features']
 for seed in range(1,11):
  table=pd.read_csv(a.results/f'{study}-{seed}-native-default.tsv.gz',sep='\t',float_precision='round_trip');selected=table.loc[table.adjustedPValue<0.05]
  if selected.empty:continue
  d=a.root/study/str(seed);members=json.loads((d/'memberships.json').read_text());ref=np.load(d/'reference.npz')['counts']
  ar=json.loads((d/'default/archive.json').read_text());entry=next(e for e in ar['entries'] if e['logicalPath']=='report.json')
  remote=Path('/Users/n/numivivo-null-benchmark-20260910')/study/str(seed)/'default/report.json.gz'
  b=subprocess.check_output(['ssh','macmini','cat '+str(remote)],timeout=60);assert hashlib.sha256(b).hexdigest()==entry['sha256']
  r=json.loads(gzip.decompress(b));c=r['contrasts'][0]
  for row in selected.to_dict('records'):
   j=ids.index(row['featureID']);diag=c['negativeBinomial']['features'][j];fit=diag['finalFit'];arms=[]
   for k,sid in enumerate(members['observationSampleIDs']):
    rows=members['groups'][sid]['sourceRows'];values=np.asarray(x[rows,j].toarray()).ravel();total=int(values.sum());assert total==int(ref[k,j])
    arms.append(dict(sampleID=sid,cells=len(rows),expressingCells=int((values>0).sum()),counts=total,largestCellCounts=int(values.max()),largestCellFraction=float(values.max()/total) if total else None))
   hits.append(dict(study=study,seed=seed,**row,arms=arms,finalDispersion=diag['finalDispersion'],geneWiseDispersion=diag.get('geneWiseDispersion'),trendDispersion=diag.get('trendDispersion'),
    dispersionOutlier=diag.get('dispersionOutlier'),maximumCooksDistance=max(fit['cooksDistances']) if fit.get('cooksDistances') else None,
    coefficients=fit['coefficients'],reportedMeans=fit['means'],reportSHA256=entry['logicalSHA256']))
a.out.write_text(json.dumps(dict(status='completed-all-native-default-sham-hit-diagnostics',postScoreDiagnosticOnly=True,fittingOrFilteringChanged=False,hits=hits),indent=2)+'\n')
print(json.dumps(dict(hits=len(hits),distinctGenes=len({(h['study'],h['featureID']) for h in hits}),maximumCooksDistance=max(h['maximumCooksDistance'] for h in hits))))
