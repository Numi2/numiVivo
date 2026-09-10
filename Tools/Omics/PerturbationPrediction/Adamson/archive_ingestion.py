#!/usr/bin/env python3
"""Retain bounded native-ingestion evidence and verify every compressed payload."""
import argparse
import gzip
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil
from prepare_ingestion import sha, write

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=False)
paths=['protocol.sha256','native-environment.log','release-build.log','focused-tests.log','regression-tests.log',
       'annotation-receipt.json','annotation-time.log','pseudobulk-receipt.json','pseudobulk-time.log',
       'pseudobulk-verify.json','pseudobulk-verify-time.log','old-annotation-rejection.log',
       'annotation-regression.log','file-bounds.log','ingestion-check.json','ingestion-check.log',
       'ingestion-check-strict.json','ingestion-check-strict.log','prepare.log','initial-whitespace-check.log',
       'prepared/annotation-plan.json','prepared/pseudobulk-plan.json','prepared/identities.json',
       'prepared/reference.npz','prepared/audit.json','pseudobulk/plan.json','pseudobulk/report.json','pseudobulk/receipt.json',
       'source/source.json','source/scverse-headers.txt','source/download.log',
       'source/metadata-retrieval.log','source/metadata-export.log','source/header-retrieval.log',
       'source/file-headers.txt','source/paper-metadata.json','source/elsevier-article.xml','source/supplementary.zip',
       'annotation-regression/report.json','file-bounds/report.json','file-bounds/rejection.log','file-bounds/plan.json']
paths+=['annotation-regression/'+f.name for f in sorted((a.work/'annotation-regression').glob('*.log'))]
entries=[]
for relative in paths:
    source=a.work/relative
    assert source.is_file(),relative
    # Raw compiler diagnostics and HTTP headers retain whitespace and CRLF.
    # Compress those bytes instead of editing evidence to satisfy patch lint.
    compressed=source.suffix in ('.log','.txt') or (source.stat().st_size>100_000 and source.suffix!='.npz')
    destination=a.out/(relative+('.gz' if compressed else ''))
    destination.parent.mkdir(parents=True,exist_ok=True)
    if compressed:
        with source.open('rb') as inp,destination.open('wb') as output:
            with gzip.GzipFile(filename='',mode='wb',fileobj=output,mtime=0) as stream: shutil.copyfileobj(inp,stream)
        with gzip.open(destination,'rb') as stream: assert hashlib.file_digest(stream,'sha256').hexdigest()==sha(source)
    else: shutil.copyfile(source,destination)
    entries.append(dict(path=str(destination.relative_to(a.out)),sha256=sha(destination),bytes=destination.stat().st_size,
                        logicalPath=relative,logicalSHA256=sha(source),logicalBytes=source.stat().st_size,compression='gzip' if compressed else None))
external=[]
for relative in ['source/original.h5ad','annotated.h5ad','numivivo']:
    source=a.work/relative
    external.append(dict(path=relative,sha256=sha(source),bytes=source.stat().st_size,retainedInGit=False))
captures=[dict(path=f.name,sha256=sha(f),bytes=f.stat().st_size) for f in sorted((a.work/'source').iterdir()) if f.is_file() and f.name!='original.h5ad']
write(a.out/'manifest.json',dict(schemaVersion=1,status='passed-native-ingestion-prediction-not-qualified',files=entries,
    externalFiles=external,sourceCaptureInventory=captures,
    python=dict(version=platform.python_version(),packages={n:version(n) for n in ['numpy','scipy','anndata','h5py','pandas']}),
    sourceCodeBase='6a83c08d0708d3236db7f3136ddd4ab92119b18b',
    sourceEdits=['Sources/NumiVivoKit/Artifacts/VivoRootedFileStore.swift','Sources/NumiVivoKit/Omics/VivoH5ADAnnotations.swift'],
    qualificationScope='Complete Adamson UPR source preservation and exact guide-group counts; no control assignment or prediction scoring.',
    retention='Large original/output H5AD and executable retained externally with hashes; native report, independent count reference, plans, receipts, validation and retrieval failure logs included. Source-code captures from external repositories are inventoried, not redistributed.'))
for e in entries: assert sha(a.out/e['path'])==e['sha256']
print(json.dumps(dict(status='verified',files=len(entries),storedBytes=sum(e['bytes'] for e in entries)),indent=2))
