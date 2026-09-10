#!/usr/bin/env python3
"""Archive the corrected identities, native results and immutable GO captures."""
import argparse
import gzip
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil
from prepare_ingestion import sha,write

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--work',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
paths=['protocol.sha256','native-environment.log','control-lookup.json','superseded-duplicate-cleanup.json',
    'source/geo-identity-retrieval.json','source/GSM2406681_10X010_barcodes.tsv.gz','source/GSM2406681_10X010_cell_identities.csv.gz',
    'descriptor-capture.log','descriptor-check.json','descriptor-check.log','descriptor-geo-compatibility.json',
    'descriptor-integrity-repeat.log','geo-restoration.log',
    'geo-restoration/annotation-plan.json','geo-restoration/annotation-receipt.json','geo-restoration/annotation-time.log',
    'geo-restoration/pseudobulk-plan.json','geo-restoration/pseudobulk-receipt.json','geo-restoration/pseudobulk-time.log',
    'geo-restoration/pseudobulk/plan.json','geo-restoration/pseudobulk/report.json','geo-restoration/pseudobulk/receipt.json',
    'geo-restoration/verify.json','geo-restoration/verify-time.log','geo-restoration/check.json','geo-restoration/check.log',
    'geo-restoration/identities.json','geo-restoration/identity-audit.json','geo-restoration/reference.npz',
    'geo-restoration/initial-pseudobulk-rejection.log','geo-restoration/initial-pseudobulk-rejection.stdout',
    'geo-restoration/release-build.log','geo-restoration/native-environment.log','geo-restoration/duplicate-cleanup.json',
    'geo-restoration/plan-bounds.log','geo-restoration/plan-bounds/report.json','geo-restoration/plan-bounds/rejection.log',
    'geo-restoration/plan-bounds/over-two-mib.json']
paths += [str(f.relative_to(a.work)) for f in sorted((a.work/'descriptors').rglob('*')) if f.is_file()]
paths += [str(f.relative_to(a.work)) for f in sorted((a.work/'descriptor-integrity-repeat').iterdir()) if f.suffix in ('.log',) or f.name=='report.json']
assert len(paths)==len(set(paths))
entries=[]
for relative in paths:
    source=a.work/relative;assert source.is_file()
    compressed=source.suffix in ('.log','.stdout') or (source.stat().st_size>100000 and source.suffix not in ('.npz','.gz'))
    destination=a.out/(relative+('.gz' if compressed else ''));destination.parent.mkdir(parents=True,exist_ok=True)
    if compressed:
        with source.open('rb') as inp,destination.open('wb') as out:
            with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as stream:shutil.copyfileobj(inp,stream)
        with gzip.open(destination,'rb') as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==sha(source)
    else:shutil.copyfile(source,destination)
    entries.append(dict(path=str(destination.relative_to(a.out)),sha256=sha(destination),bytes=destination.stat().st_size,
        logicalPath=relative,logicalSHA256=sha(source),logicalBytes=source.stat().st_size,compression='gzip' if compressed else None))
external=[]
for relative in ['source/original.h5ad','numivivo','geo-restoration/annotated.h5ad','geo-restoration/numivivo']:
    source=a.work/relative
    external.append(dict(path=relative,sha256=sha(source),bytes=source.stat().st_size,retainedInGit=False))
write(a.out/'manifest.json',dict(schemaVersion=1,status='verified-GEO-restoration-prediction-pending',files=entries,externalFiles=external,
    python=dict(version=platform.python_version(),packages={n:version(n) for n in ['numpy','scipy','anndata','h5py','pandas']}),
    sourceCodeBase='8898ea2e87d8108089b925b7ff9a4a11a5563923',
    sourceEdit='Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift',
    qualification='Restored original GEO cell-to-guide assignments and exact native counts; controls and prediction scores remain unqualified.',
    nativeArtifactRetention='Both deposited-label and restored-label source bundles retained on the Mac mini. Annotation and pseudobulk steps retain their distinct actual executable identities. Large H5ADs and binaries are external.',
    descriptors='Current GO records only; candidate identities match both source preparations. One MANF term cites this study; no temporal-independence or predictive-accuracy claim.'))
for e in entries:assert sha(a.out/e['path'])==e['sha256']
print(json.dumps(dict(status='verified',files=len(entries),storedBytes=sum(e['bytes'] for e in entries)),indent=2))
