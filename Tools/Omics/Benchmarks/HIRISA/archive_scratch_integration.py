#!/usr/bin/env python3
"""Archive integration scratch retirement, native replay and full-shape lifetime evidence."""
import argparse,errno,gzip,hashlib,json,os,shutil
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reuse-archive',type=Path,action='append',default=[])
    p.add_argument('--reuse-legacy-archive',type=Path,action='append',default=[]);a=p.parse_args()
    roots={'integration-scratch':a.root}
    def read(path):return json.loads(path.read_text())
    root=a.root
    native=read(root/'admission-validation-complete.json');fixtures=read(root/'admission-fixture-complete.json')
    real=read(root/'real-complete.json');complete=read(root/'qualification-complete.json')
    assert all(v['status']=='passed' for v in [native,fixtures,real,complete])
    assert native['binarySHA256']==fixtures['binarySHA256']==real['binarySHA256']
    assert all(c['returnCode']==0 for v in [native,fixtures,real] for c in v['commands'])
    assert '16 tests in 4 suites passed' in (root/'admission-tests.log').read_text()
    assert "Build of product 'numivivo' complete!" in (root/'admission-release-build.log').read_text()
    for name in ['admission-fixtures','mnn-fixtures']:
        check=read(root/name/'checks.json');assert check['status']=='passed' and check['binarySHA256']==native['binarySHA256']
    assert {(c['cohort'],c['seed']) for c in real['cases']}=={(name,seed) for name in ['kang','hagai'] for seed in [7,19,41]}
    assert all(c['status']=='passed' and c['nativeReplay'] and c['singleWindowDirectAccess'] for c in real['cases'])
    probe=read(root/'lifetime-probe-v2/results.json');run=read(root/'lifetime-probe-v2/run-status.json')
    assert probe['status']=='passed' and probe['cells']==1612594 and probe['clusters']==100
    assert probe['allOutputBytesIndependentlyRehashed'] and probe['baselineAndRetiredFingerprintsExact']
    assert run['returnCode']==0 and run['scratchAndOutputsRemoved']
    baseline,retired=probe['results'];assert baseline['mode']=='baseline' and retired['mode']=='retired'
    assert baseline['fingerprints']==retired['fingerprints']
    assert next(p['logicalBytes'] for p in baseline['phases'] if p['phase']=='after-assignment-scores')==8566099328
    assert next(p['logicalBytes'] for p in retired['phases'] if p['phase']=='before-serialization')==2476944384
    assert read(root/'lifetime-probe/observer-failure.json')['status']=='observer-failed'
    env=read(root/'runtime-release/environment.json');assert env['binarySHA256']==native['binarySHA256']
    assert not any(k in env['hardware'] for k in ['Serial Number','Hardware UUID','Provisioning UDID'])
    assert not any(p.name.startswith('.integration-matrix-') for p in root.rglob('*'))
    cache={};reused={}
    for archive in a.reuse_archive:
        reused[archive.name]=sha(archive/'manifest.json')
        for item in read(archive/'manifest.json')['records']:
            if item['gzipEncoded']:cache[(item['sourceSHA256'],item['sourceBytes'])]=(archive/item['storedPath'],item['storedBytes'],item['storedSHA256'])
    for archive in a.reuse_legacy_archive:
        reused[archive.name]=sha(archive/'archive-manifest.json')
        for item in read(archive/'archive-manifest.json'):
            if item['encoding']=='gzip':cache[(item['uncompressedSHA256'],item['uncompressedBytes'])]=(archive/item['path'],item['storedBytes'],item['storedSHA256'])
    a.output.mkdir(parents=True,exist_ok=False);records=[];full={};external=[]
    for label,root in roots.items():
        for source in sorted(root.rglob('*')):
            if not source.is_file():continue
            rel=source.relative_to(root)
            if '__pycache__' in rel.parts or any(v.startswith(('.numivivo-','.integration-matrix-')) for v in rel.parts):continue
            name=label+'/'+str(rel)
            if source.name in ['numivivo','lifetime-probe'] or (source.suffix=='.h5ad' and source.stat().st_size>16777216):
                external.append(dict(sourcePath=name,bytes=source.stat().st_size,SHA256=sha(source)));continue
            assert not source.is_symlink(),source
            size=source.stat().st_size;before=sha(source);full[name]=dict(bytes=size,SHA256=before)
            with source.open('rb') as f:
                offset=0;part=0
                while offset<size or (size==0 and part==0):
                    block=f.read(min(67108864,size-offset));h=hashlib.sha256(block).hexdigest()
                    stored=name+('.part-%04d'%part if size>67108864 else '')+'.gz';dest=a.output/stored;dest.parent.mkdir(parents=True,exist_ok=True)
                    prior=cache.get((h,len(block)))
                    if prior:
                        path,n,digest=prior;assert path.stat().st_size==n and sha(path)==digest
                        try:os.link(path,dest)
                        except OSError as e:
                            if e.errno!=errno.EXDEV:raise
                            shutil.copyfile(path,dest)
                    else:
                        with dest.open('xb') as out:
                            with gzip.GzipFile(fileobj=out,mode='wb',filename='',mtime=0,compresslevel=6) as z:z.write(block)
                    encoded=sha(dest);cache[(h,len(block))]=(dest,dest.stat().st_size,encoded)
                    records.append(dict(sourcePath=name,sourceOffset=offset,sourceBytes=len(block),sourceSHA256=h,storedPath=stored,storedBytes=dest.stat().st_size,storedSHA256=encoded,gzipEncoded=True))
                    offset+=len(block);part+=1
            assert source.stat().st_size==size and sha(source)==before
    manifest=dict(schemaVersion=2,records=records,fullSources=full,externalPayloads=external,archiveScriptSHA256=sha(Path(__file__)),reuseManifests=reused,
        scope='Early integration scratch retirement, 16 native tests including write-failure and cancellation ownership, full Kang/Hagai exact numerical regressions and replay, fitted/query lifecycle and tampering controls. Complete HIRISA dimensions covered by synthetic storage lifetime measurement with independently rehashed output equality. First observer race retained. Logical file extents, APFS allocated blocks and sampled transient peaks are distinct. No full-HIRISA integration, biological preservation, Metal or controlled end-to-end performance qualification.')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(members=len(records),files=len(full),storedBytes=sum(r['storedBytes'] for r in records),manifestSHA256=sha(a.output/'manifest.json'))))
if __name__=='__main__':main()
