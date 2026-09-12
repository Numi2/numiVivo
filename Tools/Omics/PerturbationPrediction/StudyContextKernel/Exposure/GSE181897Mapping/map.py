from pathlib import Path
import json,hashlib
import anndata
r=Path('/Users/n/numivivo-gse181897-exposure-mapping-20260912');r.mkdir(exist_ok=False)
p=Path('/Users/n/numivivo-gse181897-20260911');source=p/'handoff/B-source-codes-RNA.h5ad';roles=p/'prediction-inputs/primary-condition-roles.json';meta=Path('/Users/n/numivivo-parse-context-evaluation-20260912/source-metadata.json')
a=anndata.read_h5ad(source,backed='r');obs=a.obs;mapping=json.loads(roles.read_text());assert mapping['mapping']=={'B':'IFNB','C':'control'} and mapping['conditionMeaningQualified']
m=json.loads(meta.read_text());donors=m['studies']['GSE181897'];eligible={'GSE181897:exp_id:'+str(d) for d,g in obs.groupby('exp_id',observed=True) if {'B','C'}<=set(g['cond'].astype(str))};assert eligible==set(donors) and len(donors)==62
rows=[]
for donor in donors:
 d=donor.split('exp_id:')[1];conditions=[]
 for code in ['B','C']:
  g=obs[(obs['exp_id'].astype(str)==d)&(obs['cond'].astype(str)==code)];assert len(g)>0
  conditions.append({'sourceCode':code,'role':mapping['mapping'][code],'cells':len(g),'nativeSampleIDs':sorted(g['native_sample'].astype(str).unique().tolist()),'batchIDs':sorted(g['batch'].astype(str).unique().tolist()),'poolCodes':sorted(g['pool_code'].astype(str).unique().tolist()),'primaryPoolAccessions':None})
 rows.append({'modelDonorID':donor,'sourceExperimentID':d,'conditions':conditions})
h=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
out={'status':'all-eligible-donor-condition-joins-verified-primary-pool-accession-mapping-open','sourceHashes':{str(f):h(f) for f in [source,roles,meta]},'allSourceBCells':len(obs),'retainedPairedDonors':len(rows),'pairedCells':sum(c['cells'] for x in rows for c in x['conditions']),'rows':rows,'exposureCovariatesAdmitted':False}
(r/'mapping.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS',len(rows),'donors',out['pairedCells'],'paired cells')
