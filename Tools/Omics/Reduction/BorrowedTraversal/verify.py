from pathlib import Path
import hashlib,json,statistics
p=Path(__file__).resolve().parent
repo=p.parents[3]
m=json.loads((p/'runtime-manifest.json').read_text())
for name in ['VivoH5ADCountStore.swift','VivoH5ADReduction.swift']:
 source=repo/'Sources/NumiVivoKit/Omics'/name
 expected=m['/Users/n/numivivo-pca-borrowed-20260912/source/Sources/NumiVivoKit/Omics/'+name]['sha256']
 assert hashlib.sha256(source.read_bytes()).hexdigest()==expected,name
for cohort,file,prefix in [('Kang','checks.json','pca-'),('Hagai','hagai-checks.json','hagai-pca-')]:
 checks=json.loads((p/file).read_text())
 assert [x['backend'] for x in checks]==['old','new','new','old','old','new']
 assert all(x['exitCode']==0 and all(x['matches'].values()) for x in checks)
 for name in ['scores.bin','loadings.bin','metadata.json']:
  hashes={m['/Users/n/numivivo-pca-borrowed-20260912/'+prefix+str(i)+'/'+name]['sha256'] for i in range(6)}
  assert len(hashes)==1
 old,new=[statistics.median(x['seconds'] for x in checks if x['backend']==b) for b in ['old','new']]
 print(cohort,old,new,'lower elapsed fraction',1-new/old)
print('PASS: source hashes, twelve command checks and numerical artifact hashes')
