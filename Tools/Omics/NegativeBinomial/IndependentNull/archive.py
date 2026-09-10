#!/usr/bin/env python3
"""Archive complete bounded evidence, retaining failures and excluding the external H5AD."""
import argparse,gzip,hashlib,json
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();root=a.root;out=a.out
load=lambda p:json.loads(p.read_text());sha=lambda b:hashlib.sha256(b).hexdigest()
assert load(root/'summary.json')['status']=='complete-inventory'
assert load(root/'native-model-check-summary.json')['status']=='passed'
assert load(root/'native-fixed-model-check-summary.json')['status']=='passed'
assert load(root/'repair-check.json')['status']=='passed'
commands=load(root/'native-commands.json');assert [c['returnCode'] for c in commands['commands']]==[0,0]
assert '42 tests in 9 suites passed' in (root/'swift-tests-influence-fix.log').read_text()
assert 'downstreamSnapshotsAcceptTheStreamedSourceByteBudget() passed' in (root/'swift-tests-snapshot-chain.log').read_text()
out.mkdir(parents=True,exist_ok=False);files=[]
def add(source,relative):
    raw=source.read_bytes();target=Path(relative)
    if source.suffix not in ['.gz','.npz'] and len(raw)>262144:
        target=Path(str(target)+'.gz');stored=gzip.compress(raw,mtime=0)
    else:stored=raw
    destination=out/target;destination.parent.mkdir(parents=True,exist_ok=True);destination.write_bytes(stored)
    files.append(dict(path=target.as_posix(),sourcePath=str(source.relative_to(root)),sourceBytes=len(raw),
        storedBytes=len(stored),sourceSHA256=sha(raw),storedSHA256=sha(stored)))
for source in sorted(root.iterdir()):
    if source.is_file() and (source.suffix in ['.json','.npz','.gz','.log','.py','.tsv'] or source.name.endswith('.attempt-1')):
        if source.name.startswith(('available-','archive','publish')) or source.name in ['model-check-available.json','model-check-available.log','native-model-check-available.json','native-fixed-model-check-available.json','repair-check-available.json']:continue
        add(source,source.name)
for folder in ['requests','reference-inputs','reference','native','native-fixed','native-driver-attempt-1-output','build','build-influence-fix']:
    for source in sorted((root/folder).rglob('*')):
        if not source.is_file():continue
        assert source.name!='driver.lock','Cannot archive a running family'
        if folder=='native-driver-attempt-1-output' and source.name=='output.json':continue
        add(source,source.relative_to(root))
manifest=dict(schemaVersion=1,files=files,fileCount=len(files),storedBytes=sum(f['storedBytes'] for f in files),
    sourceH5AD=dict(sha256='2f019ba5ae48fdf61e02f510aac66cf6035b44aa4c582d52efce9ebccad90ee3',bytes=1199390735,
        url='https://datasets.cellxgene.cziscience.com/b6986a7f-981e-4e04-93ed-57f444749b8b.h5ad',
        retainedExternallyOnBothHosts=True),
    scope='Complete original and separately identified repair families; all diagnostics, probabilities, independent aggregate inputs and failures. Full H5AD remains external.')
(out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
for f in files:
    stored=(out/f['path']).read_bytes();assert sha(stored)==f['storedSHA256']
    raw=gzip.decompress(stored) if f['path'].endswith('.gz') and not f['sourcePath'].endswith('.gz') else stored
    assert sha(raw)==f['sourceSHA256']
print(json.dumps(dict(status='verified',files=len(files),storedBytes=manifest['storedBytes']),indent=2))
