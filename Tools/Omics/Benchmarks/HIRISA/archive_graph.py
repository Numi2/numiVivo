#!/usr/bin/env python3
"""Archive complete HIRISA graph evidence only after replay and reference gates."""
import argparse
import errno
import gzip
import hashlib
import json
import os
import shutil
from pathlib import Path
from acquire import digest

def read(path):
    return json.loads(path.read_text())

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reuse-archive',type=Path,action='append',default=[],
                   help='Reuse verified compressed chunks with identical decoded SHA256 and length.')
    a=p.parse_args();base=a.root/'graph';out=a.output
    independent=read(base/'independent-check/checks.json')
    assert read(base/'collection-complete.json')['status']=='passed'
    protocol=read(base/'protocol.json')
    assert independent['status']=='passed' and independent['cells']==1612594
    assert independent['protocolSHA256']==digest(base/'protocol.json')
    native=read(base/'native-complete.json')
    assert native['status']=='passed' and native['protocolSHA256']==digest(base/'protocol.json')
    assert all(c['returnCode']==0 for c in native['commands'])
    assert read(base/'full-graph-publish-status.json')['returnCode']==read(base/'full-graph-verify-status.json')['returnCode']==0
    assert read(base/'full-pca-regression.json')['status']=='passed'
    assert '12 tests in 4 suites passed' in (base/'scale-tests.log').read_text()
    assert "Build of product 'numivivo' complete!" in (base/'release-build.log').read_text()
    for path in ['hnsw-gate/checks.json','graph-store-gate/checks.json','baron-regression/checks.json',
                 'hagai-regression/checks.json','baron-independent-check-final/checks.json']:
        assert read(base/path)['status']=='passed',path
    for name,record in read(base/'full-output-freeze.json').items():
        path=base/'full'/name
        assert path.stat().st_size==record['bytes'] and digest(path)==record['SHA256'],name
    for name in ['scores','metadata','quality']:
        path=base/'full/input'/(name+('.bin' if name=='scores' else '.json'))
        assert digest(path)==protocol[name+'SHA256']
    assert digest(base/'reference_graph.py')==read(base/'reference-source-freeze.json')['scriptSHA256']
    # Runtime hardware was redacted before freezing; never publish device IDs.
    hardware=read(base/'runtime-release/environment.json')['hardware']
    assert not any(label in hardware for label in ['Serial Number','Hardware UUID','Provisioning UDID'])
    paths=[]
    for path in sorted(base.rglob('*')):
        if not path.is_file():continue
        relative=path.relative_to(base)
        if any(part.startswith('.numivivo-') or part=='__pycache__' for part in relative.parts):continue
        if path.name=='numivivo' or path.suffix=='.npy':continue
        if path.suffix=='.h5ad' and path.stat().st_size>16_777_216:continue
        assert not path.is_symlink(),path
        paths.append(path)
    cached={};reused=0;reuse_manifests={}
    for archive in a.reuse_archive:
        reuse_manifests[archive.name]=digest(archive/'manifest.json')
        for record in read(archive/'manifest.json')['records']:
            if record['gzipEncoded']:
                path=archive/record['storedPath']
                assert path.resolve().is_relative_to(archive.resolve())
                cached[(record['sourceSHA256'],record['sourceBytes'])]=(path,record['storedBytes'],record['storedSHA256'])
    out.mkdir(parents=True,exist_ok=False)
    records=[];sources={};chunk_bytes=64*1024*1024
    for source in paths:
        name=str(source.relative_to(base));size=source.stat().st_size;before=digest(source)
        sources[name]=dict(bytes=size,SHA256=before)
        with source.open('rb') as src:
            offset=0;part=0
            while offset<size or (size==0 and part==0):
                length=min(chunk_bytes,size-offset)
                stored=name+('.part-%04d'%part if size>chunk_bytes else '')+'.gz'
                destination=out/stored;destination.parent.mkdir(parents=True,exist_ok=True)
                # At most one 64 MiB archive chunk is resident. Identical PCA
                # witnesses need neither recompression nor duplicate disk blocks.
                block=src.read(length);assert len(block)==length
                raw_sha=hashlib.sha256(block).hexdigest();key=(raw_sha,length)
                if key in cached:
                    prior,stored_bytes,stored_sha=cached[key]
                    assert prior.stat().st_size==stored_bytes and digest(prior)==stored_sha
                    try:os.link(prior,destination)
                    except OSError as error:
                        if error.errno!=errno.EXDEV:raise
                        with prior.open('rb') as src_gzip,destination.open('xb') as dst_gzip:
                            shutil.copyfileobj(src_gzip,dst_gzip,1048576)
                    reused+=1
                else:
                    with destination.open('xb') as dst,gzip.GzipFile(fileobj=dst,mode='wb',filename='',mtime=0) as zipped:
                        zipped.write(block)
                    cached[key]=(destination,destination.stat().st_size,digest(destination))
                records.append(dict(sourcePath=name,sourceOffset=offset,sourceBytes=length,sourceSHA256=raw_sha,
                    storedPath=stored,storedBytes=destination.stat().st_size,storedSHA256=digest(destination),gzipEncoded=True))
                offset+=length;part+=1
            assert not src.read(1)
        assert source.stat().st_size==size and digest(source)==before
    manifest=dict(schemaVersion=2,records=records,fullSources=sources,sourceSHA256=protocol['sourceSHA256'],
        archiveScriptSHA256=digest(Path(__file__)),protocolSHA256=digest(base/'protocol.json'),nativeReplayPassed=True,
        reusedCompressedChunks=reused,reusedArchiveManifests=reuse_manifests,
        independentComparison=independent,
        externalPayloads='Large source H5AD files and immutable native executables remain external; source and executable identities are retained. Small H5AD fixtures, all binary graph records, input PCA artifacts, query references, checks, plans and execution/failure logs are archived.',
        scope='Complete-source approximate graph numerical qualification and frozen sampled exact recall. Original HNSW index and cell metadata remain resident. No clustering, biological preservation or controlled native/scverse/Metal performance claim.')
    (out/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    print(json.dumps(dict(members=len(records),files=len(sources),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=digest(out/'manifest.json'))))

if __name__=='__main__':main()
