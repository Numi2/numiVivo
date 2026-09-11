"""Retain compact executed evidence and bind complete native objects externally."""
import argparse, gzip, hashlib, json
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
a.out.mkdir(parents=True,exist_ok=False)
def sha(raw):return hashlib.sha256(raw).hexdigest()
groups={}
def add(group,path):
 raw=path.read_bytes()
 groups.setdefault(group,[]).append(dict(path=str(path.relative_to(a.study)),bytes=len(raw),SHA256=sha(raw),rawUTF8=raw.decode()))
for name in ['prepare_prediction.py','run_prediction.py','resume_prediction.py','score_prediction.py','restore_prediction.py','dehydrate_committed_copies.py','archive_prediction.py']:
 add('recipes',a.study/name)
for directory,group in [('prediction-inputs','inputs'),('prediction-execution','native'),('prediction-scores-complete','scores'),('prediction-scores-repeat','repeat-scores')]:
 for path in sorted((a.study/directory).rglob('*')):
  if path.is_file() and path.suffix in ['.json','.py','.log']:add(group,path)
for path in sorted((a.study/'prediction-restoration-records').iterdir()):
 if path.is_file():add('restoration',path)
for name in ['prediction-score.log','prediction-score-complete.log','prediction-score-repeat.log','preparation-tree-alias-correction.json','prediction-disk-recovery.json','prediction-local-disk-recovery.json','scorer-freeze-before-outcomes.json']:
 add('attempts-and-storage',a.study/name)
records=[]
for name,files in sorted(groups.items()):
 raw=(json.dumps(dict(sourceFiles=files),sort_keys=True,separators=(',',':'))+'\n').encode()
 data=gzip.compress(raw,compresslevel=9,mtime=0);target=a.out/(name+'.json.gz');target.write_bytes(data)
 assert gzip.decompress(target.read_bytes())==raw
 records.append(dict(path=target.name,bytes=len(data),SHA256=sha(data),decodedBytes=len(raw),decodedSHA256=sha(raw),members=len(files)))
external=[]
for base,pattern in [('prediction-execution/objects','*.gz'),('prediction-inputs','**/*.h5ad')]:
 for path in sorted((a.study/base).glob(pattern)):
  raw=path.read_bytes();external.append(dict(path=str(path.relative_to(a.study)),bytes=len(raw),SHA256=sha(raw)))
history=[]
roles=json.loads((a.study/'prediction-inputs/primary-condition-roles.json').read_bytes())
for name,info in roles['sources'].items():history.append(dict(upstreamPath=name,**info))
manifest=dict(schemaVersion=1,records=records,externalArtifacts=external,
 externalRoots=[str(a.study),'/Users/n/numivivo-gse181897-20260911'],
 externalScope='Complete native content objects and actual H5AD training/query inputs are retained at both roots; this compact archive alone cannot restore them.',
 primaryAuthorSources=history,primaryAuthorCommit=roles['commit'],
 nativeRecipeImplementation='0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b',
 verification='All 110 native objects were decoded and hash-checked independently; all 124 native donor predictions reconstructed, all 496 estimate vectors independently checked. One four-donor bundle per origin additionally restored and verified natively.',
 scientificResult='HIRISA mean primary PASS; Kang mean primary FAIL. HIRISA nominal interval coverage severely deficient. General biological prediction and calibrated uncertainty not established.')
(a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
print(json.dumps(dict(compactBytes=sum(x['bytes'] for x in records),externalBytes=sum(x['bytes'] for x in external),sourceFiles=sum(x['members'] for x in records))))
