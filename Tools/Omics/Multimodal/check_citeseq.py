#!/usr/bin/env python3
"""Qualify native multi-assay product on the complete public PBMC5k CITE-seq matrix."""
import argparse
import copy
import hashlib
import importlib.metadata
import json
import platform
import shutil
import subprocess
import time
from pathlib import Path
import h5py
import mudata
import numpy as np
import scanpy as sc
from scipy import sparse

SOURCE_SHA = '3b290ad9605b96974c9c16e5ae3427e5e5c496a66e55223df135b388b6d61417'

def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write(p, value):
    p.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--source', type=Path, required=True)
p.add_argument('--plan', type=Path, default=Path(__file__).with_name('pbmc5k-plan.json'))
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
assert sha(a.source) == SOURCE_SHA
a.out.mkdir(parents=True, exist_ok=False)
commands = []
def run(label, args, success=True):
    start = time.monotonic()
    r = subprocess.run([str(a.binary), *map(str, args)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (a.out / (label + '.log')).write_text(r.stdout)
    commands.append(dict(label=label, arguments=list(map(str,args)), expectedSuccess=success, exitCode=r.returncode, seconds=time.monotonic()-start))
    write(a.out / 'commands.json', commands)
    assert (r.returncode == 0) == success, (label, r.stdout[-3000:])

bundle = a.out / 'bundle'
run('import', ['multiassay-10x-import', a.source, '--plan', a.plan, '--output', bundle])
run('verify', ['multiassay-verify', bundle])
run('overwrite', ['multiassay-10x-import', a.source, '--plan', a.plan, '--output', bundle], False)
run('repeat', ['multiassay-10x-import', a.source, '--plan', a.plan, '--output', a.out/'repeat'])
for name in ['original.h5','plan.json','dataset.json','dataset.h5mu','receipt.json']:
    assert sha(bundle/name) == sha(a.out/'repeat'/name), name

with h5py.File(a.source) as f:
    m=f['matrix']
    original=sparse.csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr()
    ids=m['features/id'].asstr()[:]; names=m['features/name'].asstr()[:]; types=m['features/feature_type'].asstr()[:]
    barcodes=m['barcodes'].asstr()[:]
scan=sc.read_10x_h5(a.source,gex_only=False)
np.testing.assert_array_equal(scan.obs_names,barcodes)
np.testing.assert_array_equal(scan.var['gene_ids'],ids)
assert (scan.X != original).nnz == 0
mu=mudata.read_h5mu(bundle/'dataset.h5mu')
np.testing.assert_array_equal(mu.obs['barcode'],barcodes)
with (bundle/'dataset.json').open() as f: native=json.load(f)
assert len(native['observations']) == 5247 and len(native['assays']) == 2
assert [o['identity']['barcode'] for o in native['observations']] == barcodes.tolist()
results=[]
for key,feature_type in [('rna','Gene Expression'),('protein','Antibody Capture')]:
    ix=np.flatnonzero(types==feature_type)
    assay=next(x for x in native['assays'] if x['id']==key)
    assert assay['observationIndices']==list(range(5247))
    assert [f['id'] for f in assay['features']]==ids[ix].tolist()
    assert [f['name'] for f in assay['features']]==names[ix].tolist()
    d=assay['matrix']
    matrix=sparse.csr_matrix((np.array(d['counts'],dtype=np.uint64),d['featureIndices'],d['rowOffsets']),shape=(d['cellCount'],d['featureCount']))
    expected=original[:,ix]
    assert (matrix!=expected).nnz==0 and (mu.mod[key].X!=expected).nnz==0
    np.testing.assert_array_equal(mu.mod[key].var_names,ids[ix])
    np.testing.assert_array_equal(mu.mod[key].obs['barcode'],barcodes)
    np.testing.assert_array_equal(np.asarray(mu.obsmap[key]).ravel(),np.arange(1,5248))
    np.testing.assert_array_equal(np.asarray(matrix.sum(axis=1)).ravel(),np.asarray(expected.sum(axis=1)).ravel())
    results.append(dict(assay=key,cells=5247,features=len(ix),nonzeros=matrix.nnz,totalCounts=int(matrix.sum()),
                        sourceCountsExact=True,scanpyCountsExact=True,mudataCountsExact=True))
# Independent MuData write/read retains both modality matrices and alignment.
mu.write(a.out/'mudata-roundtrip.h5mu')
roundtrip=mudata.read_h5mu(a.out/'mudata-roundtrip.h5mu')
for key in mu.mod:
    assert (roundtrip.mod[key].X!=mu.mod[key].X).nnz==0
    np.testing.assert_array_equal(roundtrip.mod[key].obs_names,mu.mod[key].obs_names)

# Structural fixture uses existing typed metadata, but is explicitly synthetic.
fixture=dict(schemaVersion=1,id='partial-fixture',evidence='synthetic',sourceDescription='Structural control only',
    samples=native['samples'],observations=[dict(identity=dict(barcode='b'+str(i),sampleID=native['samples'][0]['id']),kind='spot' if i==2 else 'cell',position=dict(frameID='slide',coordinates=[2.5,4.]) if i==2 else None) for i in range(3)],
    spatialFrames=[dict(id='slide',unit='micrometer',axes=['x','y'],sourceDescription='Synthetic coordinate frame')],assays=[])
for i,(key,kind,rows,offsets,counts_) in enumerate([('rna','rna',[0,1,2],[0,1,1,2],[2**53+1,4]),('protein','antibodyCapture',[2,0],[0,1,1],[7])]):
    fixture['assays'].append(dict(id=key,kind=kind,featureNamespace='fixture-'+key,countUnit='umiCount',genomeAssembly=None,sourceDescription='Synthetic fixture',
        features=[dict(id='shared-id',name=key,interval=None)],observationIndices=rows,matrix=dict(cellCount=len(rows),featureCount=1,rowOffsets=offsets,featureIndices=[0]*len(counts_),counts=counts_)))
write(a.out/'fixture.json',fixture)
run('partial-export',['multiassay-h5mu-write',a.out/'fixture.json','--output',a.out/'partial.h5mu'])
partial=mudata.read_h5mu(a.out/'partial.h5mu')
np.testing.assert_array_equal(np.asarray(partial.obsmap['protein']).ravel(),[2,0,1])
assert int(partial.mod['rna'].X[0,0])==2**53+1
assert int(partial.mod['protein'].X[0,0])==7 and int(partial.mod['protein'].X[1,0])==0
metadata=json.loads(partial.uns['numivivo'])
assert metadata['observations'][2]['position']==dict(frameID='slide',coordinates=[2.5,4.])
assert metadata['spatialFrames'][0]['unit']=='micrometer'
run('export-overwrite',['multiassay-h5mu-write',a.out/'fixture.json','--output',a.out/'partial.h5mu'],False)
for label,mutation in [
    ('duplicate-observation',lambda d:d['assays'][1].update(observationIndices=[0,0])),
    ('out-of-range-observation',lambda d:d['assays'][1].update(observationIndices=[3,0])),
    ('unknown-spatial-frame',lambda d:d['observations'][2]['position'].update(frameID='missing')),
    ('unknown-unit',lambda d:d['assays'][1].update(countUnit='concentration'))]:
    bad=copy.deepcopy(fixture);mutation(bad);write(a.out/(label+'.json'),bad)
    run(label,['multiassay-h5mu-write',a.out/(label+'.json'),'--output',a.out/(label+'.h5mu')],False)
    assert not (a.out/(label+'.h5mu')).exists()
for label,mutation in [
    ('omit-protein',lambda d:d.update(assays=d['assays'][:1])),
    ('wrong-assembly',lambda d:d['assays'][0].update(genomeAssembly='GRCh37')),
    ('wrong-kind',lambda d:d['assays'][1].update(kind='rna')),
    ('unknown-option',lambda d:d.update(ignoreMissing=True))]:
    bad=json.loads(a.plan.read_text());mutation(bad);write(a.out/(label+'.json'),bad)
    run(label,['multiassay-10x-import',a.source,'--plan',a.out/(label+'.json'),'--output',a.out/label],False)
    assert not (a.out/label).exists()
# Rehashed count tampering must fail reconstruction, not merely digest comparison.
tampered=a.out/'tampered';tampered.mkdir()
for name in ['original.h5','plan.json','dataset.h5mu']:
    shutil.copy2(bundle/name,tampered/name)
assert json.dumps(native,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode() == (bundle/'dataset.json').read_bytes()
changed=copy.deepcopy(native);changed['assays'][0]['matrix']['counts'][0]+=1
# Native canonical JSON uses sorted keys and compact separators.
(tampered/'dataset.json').write_text(json.dumps(changed,sort_keys=True,separators=(',',':'),ensure_ascii=False))
receipt=json.loads((bundle/'receipt.json').read_text())
# Fingerprint encoding is verified from the actual receipt, never guessed.
assert set(receipt['dataset'])=={'bytes'} and len(receipt['dataset']['bytes'])==32
receipt['dataset']={'bytes':list(bytes.fromhex(sha(tampered/'dataset.json')))};write(tampered/'receipt.json',receipt)
run('rehashed-dataset-tamper',['multiassay-verify',tampered],False)
write(a.out/'checks.json',dict(status='passed',sourceSHA256=SOURCE_SHA,binarySHA256=sha(a.binary),
    protocol='complete-source-exact-count-and-alignment',assays=results,repeatBytesExact=True,mudataRoundtripExact=True,
    partialAssayMapExact=True,aboveFloatIntegerPrecisionPreserved=True,spatialMetadataPreserved=True,
    commands=len(commands),expectedRejections=sum(not r['expectedSuccess'] for r in commands),
    platform=platform.platform(),packages={n:importlib.metadata.version(n) for n in ['numpy','scipy','scanpy','anndata','mudata','h5py']},
    implementationSHA256=sha(Path(__file__)),planSHA256=sha(a.plan)))
