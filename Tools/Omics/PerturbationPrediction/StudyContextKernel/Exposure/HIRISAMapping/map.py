from pathlib import Path
import json,hashlib
r=Path('/Users/n/numivivo-hirisa-exposure-mapping-20260912');r.mkdir(exist_ok=False)
p=Path('/Users/n/numivivo-hirisa-20260910/design.json');m=Path('/Users/n/numivivo-parse-context-evaluation-20260912/source-metadata.json')
a=json.loads(p.read_text());meta=json.loads(m.read_text());samples={s['accession']:s for s in a['samples']};matches=[c for c in a['comparisons'] if c['population']=='Bcell' and c['treatment']=='IFNb'];assert len(matches)==1
pairs=matches[0]['pairs'];assert {'HIRISA:'+v['donor'] for v in pairs}==set(meta['studies']['HIRISA']) and len(pairs)==5
out=[]
for v in pairs:
 row=dict(v,modelDonorID='HIRISA:'+v['donor'],dose=100,doseUnit='units/mL',durationHours=21,vendor='R&D',catalogue='8499-IF-010/CF');evidence={}
 for role,condition in [('control','none'),('treated','IFNb')]:
  s=samples[v[role]]
  for key,expected in [('subject id',v['donor']),('cell type','Bcell'),('treatment',condition),('treatment time','21 hours'),('batch id',v['batch']),('pool id',v['pool'])]:assert s[key]==[expected],(v,role,key)
  text=' '.join(s['treatment_protocol_ch1']);assert 'IFN-b (100 units/mL, R&D #8499-IF-010/CF)' in text
  evidence[role]={'accession':s['accession'],'sourceURL':'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+s['accession'],'fields':{k:s[k] for k in ['subject id','cell type','treatment','treatment time','batch id','pool id']}}
 row['evidence']=evidence;out.append(row)
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
result={'scope':'Retained deposited-metadata join for five HIRISA B-cell IFNb/control donor pairs; no count reprocessing or prediction','inputs':{str(p):sha(p),str(m):sha(m)},'upstreamSOFTSHA256':a['sourceMetadataSHA256'],'rows':out,'allFiveDonorsAndTenLibrariesMatched':True,'unitConversionPerformed':False}
(r/'mapping.json').write_text(json.dumps(result,indent=2)+'\n');print('PASS: all five model donors and ten source libraries matched')
