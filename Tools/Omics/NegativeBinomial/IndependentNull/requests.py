#!/usr/bin/env python3
"""Create every frozen native request without inspecting differential-expression results."""
import argparse,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
a=p.parse_args();root=a.root;cohorts=json.loads((root/'cohorts.json').read_text())
assert cohorts['protocolSHA256']==hashlib.sha256(Path(__file__).with_name('PROTOCOL.md').read_bytes()).hexdigest()
directory=root/'requests';directory.mkdir(exist_ok=False);manifest=[]
for cohort in cohorts['cohorts']:
    for name,method in [('wald',None),('lrt','likelihoodRatio'),('ql','quasiLikelihoodAdjusted')]:
        options=dict(trend='gammaParametric',minimumTrendGenes=20,minimumPriorVariance=.25,outlierStandardDeviations=2)
        if method:options['testMethod']=method
        request=dict(id='independent-null-'+cohort['id']+'-'+name,model='negativeBinomial',controlCondition='shamA',
            treatmentCondition='shamB',cellGroup='B cell',design='independentReplicates',sizeFactors='medianRatio',
            adjustForBatch=True,minimumCellsPerPseudobulk=10,minimumReplicatesPerCondition=3,minimumFeatureCounts=10,
            minimumExpressingPseudobulks=3,negativeBinomialOptions=options,
            includedDonorIDs=sorted(d['donorID'] for d in cohort['donors']))
        raw=(json.dumps(request,sort_keys=True,indent=2)+'\n').encode();filename=cohort['id']+'-'+name+'.json'
        (directory/filename).write_bytes(raw);manifest.append(dict(cohort=cohort['id'],method=name,path=filename,requestSHA256=hashlib.sha256(raw).hexdigest()))
(directory/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n');print('Prepared',len(manifest),'fixed requests')
