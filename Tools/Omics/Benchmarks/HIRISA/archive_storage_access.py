#!/usr/bin/env python3
"""Archive qualified CSR reader changes and the separate integration storage probe."""
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
    p.add_argument('--reuse-archive',type=Path,action='append',default=[]);a=p.parse_args()
    roots={'clustering':a.root/'clustering-storage-adaptive','first-candidate':a.root/'clustering-storage','integration-prototype':a.root/'integration-storage'}
    def read(path):return json.loads(path.read_text())
    for key in ['clustering','first-candidate']:
        root=roots[key];status=read(root/'small-native-complete.json')
        assert status['status']=='passed' and all(c['returnCode']==0 for c in status['commands'])
        assert read(root/'probe-results.json')['allEdgeRecordsBitExact']==34707084
        assert read(root/'probe-input-identity.json')['status']=='passed'
        assert "Build of product 'numivivo' complete!" in (root/'release-build.log').read_text()
        assert ('12 tests' if key=='clustering' else '11 tests')+' in 3 suites passed' in (root/'native-tests.log').read_text()
        for name in ['baron-exact.json','hagai-exact.json']:
            check=read(root/name);assert check['status']=='passed' and check['allFrozenResultBytesExact'] and check['graphRecordsBitExact']
        for name in ['graph-fixtures','clustering-fixtures']:
            assert read(root/name/'checks.json')['status']=='passed'
        env=read(root/'runtime-release/environment.json');assert env['binarySHA256']==status['binarySHA256']
        assert not any(k in env['hardware'] for k in ['Serial Number','Hardware UUID','Provisioning UDID'])
    assert read(roots['integration-prototype']/'probe-results.json')['status']=='passed'
    assert not any(p.name.startswith('.integration-matrix-') for p in roots['integration-prototype'].iterdir())
    cache={};reused={}
    for archive in a.reuse_archive:
        reused[archive.name]=sha(archive/'manifest.json')
        for item in read(archive/'manifest.json')['records']:
            if item['gzipEncoded']:cache[(item['sourceSHA256'],item['sourceBytes'])]=(archive/item['storedPath'],item['storedBytes'],item['storedSHA256'])
    a.output.mkdir(parents=True,exist_ok=False);records=[];full={};external=[]
    for label,root in roots.items():
        for source in sorted(root.rglob('*')):
            if not source.is_file():continue
            rel=source.relative_to(root)
            if '__pycache__' in rel.parts or any(v.startswith(('.numivivo-','.integration-matrix-')) for v in rel.parts):continue
            name=label+'/'+str(rel)
            if source.name in ['numivivo','storage-probe','storage-probe-tiled'] or (source.suffix=='.h5ad' and source.stat().st_size>16777216):
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
        scope='All original HIRISA graph edge records checked bit for bit for the new bounded reader; native full Baron/Hagai results and lifecycles preserved. Original million-cell clustering is a separate ongoing frozen run. Integration data are a synthetic storage-only full-shape probe; no complete integration, biology or controlled end-to-end performance claim.')
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(members=len(records),files=len(full),storedBytes=sum(r['storedBytes'] for r in records),manifestSHA256=sha(a.output/'manifest.json'))))
if __name__=='__main__':main()
