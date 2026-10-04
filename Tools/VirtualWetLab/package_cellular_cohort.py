#!/usr/bin/env python3
"""Verify native composite corpora and register measured exploration, not a model."""
import argparse,json,shutil,collections
from pathlib import Path
from wetlab import read,write,sha,require,inventory
from prepare_cellular_cohort import native,TASKS

def package(registration,cohort,binary,metadata):
 root=Path(cohort).resolve();reg=read(registration);(root/'source-metadata').mkdir();sources=[];specimens=[];populations=[];conditions=[];feature_names=set();masks={}
 for s in reg['sources']:
  metadata_path=Path(metadata)/(s['id']+'.json');require(sha(metadata_path)==s['metadataSHA256'],'Frozen source metadata changed');shutil.copy2(metadata_path,root/'source-metadata'/metadata_path.name)
  d=root/s['id'];report=read(d/'admission.json');provenance=read(d/'source.json');store_metadata=read(d/'store/metadata.json');rows=read(d/'rows.json')
  require(provenance['registrationSHA256']==sha(registration),'Admission registration differs')
  observed={f['id'] for f in store_metadata['features']};masks[s['id']]=[f in observed for f in reg['featureAxis']];feature_names.update(f['name'] for f in store_metadata['features'])
  src={**s,'directory':s['id'],'populationID':s['id']+'|'+s['context'],'conditionID':s['conditionID'],'sha256':provenance['verifiedSHA256'],'admission':report};sources.append(src)
  specimens.append({'id':s['id'],'title':s['study']+' · '+s['context'],'sourceSHA256':src['sha256'],'inferenceSupported':False})
  populations.append({'id':src['populationID'],'title':s['context']+' · '+str(len(rows))+' admitted cells','specimenID':s['id'],'conditionID':s['conditionID'],'controlCells':report['controls'],'biologicalUnits':None})
  conditions.append({'id':s['conditionID'],'title':s['context']+' · day '+str(s['daysPostTransduction'])+' post-transduction'})
 require(len({s['sha256'] for s in sources})==len(sources),'Duplicate source bytes cannot become independent sources')
 eligibility=[]
 for s in sources:
  counts=collections.Counter(read(root/'source-metadata'/(s['id']+'.json'))['obs']['gene'])
  eligibility.extend({'source':s['id'],'target':target,'sourceCells':counts[target],'eligible':counts[target]>=reg['support']['minimumCells'],'reason':'metadata support met' if counts[target]>=reg['support']['minimumCells'] else 'below registered per-source cell support'} for target in reg['targets'])
 write(root/'source-eligibility.json',eligibility)
 reports={}
 for task in TASKS:
  manifest={'schemaVersion':1,'format':'numivivo-cell-response-cli-sources/v1','sources':[{'corpus':str(root/s['id']/('corpus-'+task)),'store':str(root/s['id']/'store')} for s in sources]}
  write(root/(task+'-sources.json'),manifest);reports[task]=native(binary,'cell-response-sources-verify',root/(task+'-sources.json'),output=root/(task+'-verified.json'))
 shutil.copy2(registration,root/'registration.json');write(root/'feature-masks.json',masks);write(root/'features.json',reg['featureAxis'])
 write(root/'execution.json',{'nativeBinarySHA256':sha(binary),'owners':{name:sha(Path(__file__).with_name(name)) for name in ['prepare_cellular_cohort.py','package_cellular_cohort.py','measured_cellular_cohort.py']},'trainingExecuted':False})
 access=[{'source':s['id'],'rowsSHA256':sha(root/s['id']/'rows.json')} for s in sources];write(root/'observation-access.json',access)
 # Exclude invocation path manifests from scientific identity; they are location
 # maps and can be regenerated after relocation without rewriting native receipts.
 files=inventory(root);files={k:v for k,v in files.items() if not k.endswith('-sources.json')}
 names=sorted(feature_names);a={'format':'numivivo-measured-cellular-cohort/v1','id':'replogle-nadig-engineering-256','title':'Replogle–Nadig · measured engineering cohort','family':'measured-cellular-cohort','presentation':'population','geometry':None,'sources':sources,'specimens':specimens,'populations':populations,'conditions':conditions,'targets':[{'target':t} for t in reg['targets']],'features':names[:128],'featuresTotal':len(names),'datasetFingerprint':sha(root/'registration.json'),'observationAccessFingerprint':sha(root/'observation-access.json'),'evidence':'MEASURED','biologicalPromotion':False,'limits':[reg['splitScope'],reg['sourcePanelLimitation'],'No trained compatible predictor; constant corpus descriptor is admission-only, not biological target information','Development-only reservation, not independent biological validation'],'files':files}
 write(root/'assay.json',a)
 write(root/'qualification.json',{'nativeCompositeVerified':list(reports),'sources':[{k:s['admission'][k] for k in ('source','cells','controls','targets','features')} for s in sources],'unionFeatures':len(reg['featureAxis']),'structurallyMissingFeatures':{sid:len(v)-sum(v) for sid,v in masks.items()},'sampledSourceRowsExactlyReconstructed':sum(len(read(root/s['id']/'row-reconstruction.json')['samples']) for s in sources),'biologicalPromotion':False,'trainingExecuted':False,'sourceIdentitiesResolvedTo':'original study + cell context + recorded capture; independent biological units unknown'})
 return root/'assay.json'
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--registration',required=True);p.add_argument('--cohort',required=True);p.add_argument('--binary',required=True);p.add_argument('--metadata',required=True);a=p.parse_args();print(package(**vars(a)))
