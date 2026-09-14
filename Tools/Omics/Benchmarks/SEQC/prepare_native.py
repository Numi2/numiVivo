#!/usr/bin/env python3
"""Bind complete technical-library counts to the unchanged native inference owner."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
a = p.parse_args(); r = a.root
out = r / 'native-inputs'; out.mkdir(exist_ok=False)
h = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
manifest = []
for site in ['AGR','BGI','CNL','COH','MAY','NVS']:
    directory = r / 'counts' / site
    feature = pd.read_csv(directory/'features.tsv', sep='\t', dtype=str, keep_default_na=False)
    counts = pd.read_csv(directory/'counts.tsv.gz', sep='\t')
    libraries = pd.read_csv(directory/'libraries.tsv', sep='\t', dtype=str)
    assert counts.columns.tolist() == libraries.library.tolist()
    x = counts.to_numpy(); assert np.isfinite(x).all() and (x>=0).all() and (x==np.floor(x)).all() and (x<2**53).all()
    assert len(feature)==len(x)==25794
    ids = [f'SEQC:RefSeq:row:{i}' for i in range(len(feature))]
    samples=[]; cells=[]; groups=[]; offsets=[0]; columns=[]; values=[]
    for i, library in enumerate(libraries.itertuples(index=False)):
        sample = f'SEQC:{site}:{library.library}'; replicate = 'technical-library:'+sample
        samples.append(dict(id=sample,biologicalReplicateID=replicate,condition=library.condition,batchID=site,organism='Homo sapiens'))
        cells.append(dict(barcode=sample,sampleID=sample))
        groups.append(dict(biologicalReplicateID=replicate,condition=library.condition,organism='Homo sapiens',
                           sampleIDs=[sample],batchIDs=[site],sourceCellIndices=[i]))
        indices=np.flatnonzero(x[:,i]);columns.extend(indices.tolist());values.extend(x[indices,i].astype(np.uint64).tolist());offsets.append(len(values))
    assert int(sum(values))==int(x.sum())
    metadata=dict(id='SEQC-'+site,evidence='measured',countUnit='fragmentCount',samples=samples,cells=cells,
                  sourceDescription='Bulk RNA technical library preparations. Legacy cells and biologicalReplicateID fields identify matrix rows and technical units only; no single-cell or biological-donor replication claim.',
                  features=[dict(id=ids[i],name=row.Symbol if row.Symbol not in ('','NA') else ids[i],mitochondrial=False) for i,row in enumerate(feature.itertuples(index=False))])
    bulk=dict(method='SEQC-technical-library-aggregate-v1',countUnit='fragmentCount',groups=groups,featureIDs=ids,
              matrix=dict(cellCount=len(samples),featureCount=len(ids),rowOffsets=offsets,featureIndices=columns,counts=values))
    payload=dict(metadata=metadata,pseudobulk=bulk)
    dest=out/(site+'.json.gz');dest.write_bytes(gzip.compress(json.dumps(payload,sort_keys=True,separators=(',',':')).encode(),mtime=0))
    request=dict(id='SEQC-'+site+'-B-vs-A',model='negativeBinomial',controlCondition='A',treatmentCondition='B',
                 design='independentReplicates',sizeFactors='medianRatio',adjustForBatch=False,minimumCellsPerPseudobulk=1,
                 minimumReplicatesPerCondition=3,minimumFeatureCounts=10,minimumExpressingPseudobulks=3,
                 negativeBinomialOptions=dict(trend='gammaParametric',minimumTrendGenes=20,minimumPriorVariance=.25,
                    outlierStandardDeviations=2,effectPriorEstimation='weightedUpperQuantile'))
    req=out/(site+'-request.json');req.write_text(json.dumps(request,sort_keys=True,indent=2))
    manifest.append(dict(site=site,features=len(ids),technicalLibraries=len(samples),fragmentCounts=int(x.sum()),
                         sourceFiles={p.name:h(p) for p in directory.iterdir()},inputSHA256=h(dest),requestSHA256=h(req)))
(out/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
