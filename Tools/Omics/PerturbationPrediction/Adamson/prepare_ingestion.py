#!/usr/bin/env python3
"""Audit all source counts and prepare lossless guide-label ingestion, without fitting.

No control assignments or biological target names are inferred. Missing labels get
an explicit new category in a derived column; source obs and X remain untouched.
SciPy independently aggregates bounded CSC gene blocks into source-guide groups.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import h5py
import numpy as np
from scipy import sparse

SOURCE = 'e70fcd49808cab8d724de8d5a332940911206e1c8ef44cc7b568d048ed795c85'
MISSING = 'numivivo-unassigned-source-guide'
COLUMN = 'numivivo_source_guide'

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def write(path, value):
    path.write_text(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n')

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    assert sha(a.source) == SOURCE
    with h5py.File(a.source, 'r') as f:
        x = f['X']; cells, features = map(int, x.attrs['shape'])
        assert (cells, features) == (65337, 32738)
        assert x.attrs['encoding-type'] == 'csc_matrix'
        categories = f['obs/perturbation/categories'].asstr()[:].tolist()
        source_codes = f['obs/perturbation/codes'][:].astype(np.int64)
        assert MISSING not in categories and COLUMN not in f['obs']
        assert np.all((source_codes >= -1) & (source_codes < len(categories)))
        categories.append(MISSING)
        codes = np.where(source_codes == -1, len(categories)-1, source_codes)
        ids = f['var/ensembl_id'].asstr()[:].tolist()
        names = f['var/gene_symbol'].asstr()[:].tolist()
        barcodes = f['obs/cell_barcode'].asstr()[:].tolist()
        assert len(set(ids)) == features and len(set(barcodes)) == cells
        mito = [i for i, name in enumerate(names) if name.startswith('MT-')]
        labels = sorted(categories); group_index = {v: i for i, v in enumerate(labels)}
        rows = np.array([group_index[categories[i]] for i in codes])
        membership = sparse.csr_matrix((np.ones(cells, dtype=np.int64), (rows, np.arange(cells))), shape=(len(labels), cells))
        counts = np.zeros((len(labels), features), dtype=np.int64)
        totals = np.zeros(cells, dtype=np.int64); detected = np.zeros(cells, dtype=np.int64)
        mitochondrial = np.zeros(cells, dtype=np.int64)
        pointers = x['indptr'][:].astype(np.int64)
        assert pointers[0] == 0 and pointers[-1] == len(x['data']) == len(x['indices'])
        assert np.all(np.diff(pointers) >= 0)
        zero_entries = canonical = 0; maximum = 0
        for start in range(0, features, 64):
            stop = min(features, start+64); lo, hi = pointers[start], pointers[stop]
            values = x['data'][lo:hi]; indices = x['indices'][lo:hi]
            assert np.isfinite(values).all() and np.all(values >= 0)
            assert np.all(values == np.floor(values)) and np.all(values <= 2**24)
            assert np.all((indices >= 0) & (indices < cells))
            maximum = max(maximum, int(values.max(initial=0)))
            zero_entries += int(np.count_nonzero(values == 0))
            block = sparse.csc_matrix((values.astype(np.int64), indices, pointers[start:stop+1]-lo), shape=(cells, stop-start))
            block.sum_duplicates(); block.eliminate_zeros(); canonical += block.nnz
            counts[:, start:stop] = (membership @ block).toarray()
            totals += np.asarray(block.sum(axis=1)).ravel()
            detected += np.asarray((block > 0).sum(axis=1)).ravel()
            selected = [i-start for i in mito if start <= i < stop]
            if selected: mitochondrial += np.asarray(block[:, selected].sum(axis=1)).ravel()
        assert int(counts.sum()) == int(totals.sum()) and int(detected.sum()) == canonical
        np.savez_compressed(a.out/'reference.npz', counts=counts, totalCounts=totals, detectedFeatures=detected,
                            mitochondrialCounts=mitochondrial, sourceCodes=source_codes, derivedCodes=codes)
        write(a.out/'identities.json', dict(groups=labels, features=ids, names=names, barcodes=barcodes))
        plan = dict(schemaVersion=1, source=dict(bytes=list(bytes.fromhex(SOURCE))),
                    provenance='Adamson GSM2406681 complete scPerturb source. Preserve every original guide label and all source cells; source missing categorical codes become an explicit unassigned category only in this added column. No control or target identity inference.',
                    edits=[dict(path='obs/'+COLUMN, mode='add', value=dict(categorical=dict(codes=codes.tolist(), categories=categories, ordered=False)))])
        write(a.out/'annotation-plan.json', plan)
        samples = [dict(id=s, condition=s, biologicalReplicateID='K562-pooled-replication-unresolved',
                        batchID='GSM2406681-10X010-pooled', organism='NCBITaxon:9606') for s in labels]
        mapping = dict(schemaVersion=1, id='adamson2016-upr-all-source-guides', evidence='measured', countUnit='umiCount', matrixPath='X',
                       sourceDescription='Adamson 2016 GSM2406681 UPR CRISPRi screen, complete scPerturb mirror with published MD5 verified. Aggregated by original source guide labels; missing labels retained explicitly. No biological replication or control identity claim.',
                       sampleColumn=COLUMN, featureIDColumn='ensembl_id', groupColumn='cell_line', samples=samples,
                       mitochondrialFeatureIDs=[ids[i] for i in mito])
        write(a.out/'pseudobulk-plan.json', dict(schemaVersion=1, mapping=mapping, contrasts=[]))
        sizes = Counter(categories[i] for i in codes)
        write(a.out/'audit.json', dict(sourceSHA256=SOURCE, cells=cells, features=features, storedEntries=len(x['data']),
              canonicalNonzeros=canonical, explicitZeroEntries=zero_entries, maximumStoredCount=maximum,
              totalCounts=int(totals.sum()), zeroLibraryCells=int(np.count_nonzero(totals == 0)), groups=len(labels),
              groupCellCounts=dict(sorted(sizes.items())), missingGuideCells=int(np.count_nonzero(source_codes == -1)),
              sourceNpertsCounts={str(int(k)): int(v) for k,v in zip(*np.unique(f['obs/nperts'][:], return_counts=True))},
              sourceNpertsIsTargetMultiplicity=False, controlIdentitiesVerified=False, responseFittingPerformed=False,
              aggregation='SciPy int64 sparse membership times CSC blocks of at most 64 genes; no dense cells-by-genes matrix',
              sourceSlots={k:list(f[k]) for k in ['layers','obsm','obsp','varm','varp','uns']}, elapsedSeconds=time.monotonic()-started))
    print((a.out/'audit.json').read_text())

if __name__ == '__main__': main()
