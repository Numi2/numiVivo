#!/usr/bin/env python3
"""Audit all four original Baron human matrices and prepare lossless sparse H5AD.

No cell/gene filtering, rounding, donor pooling or disease relabeling. Mouse
matrices belong to a different organism and are not merged into human counts.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ARCHIVE_SHA = 'aed2d208d47a36658aa0e63629afe5d4144ef465a8e3d9a0f377422b1f1073dc'
SOFT_SHA = 'd630e7e35aa7e98d78a23b5304ecf54d65e5c3bebf2568be386fe6d3d6f09eb9'
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84133/'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--soft', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    assert digest(a.archive) == ARCHIVE_SHA, 'Unqualified count archive'
    assert digest(a.soft) == SOFT_SHA, 'Unqualified GEO metadata release'
    a.out.mkdir(parents=True, exist_ok=False)
    samples = {}
    sample = None
    for line in gzip.decompress(a.soft.read_bytes()).decode().splitlines():
        if line.startswith('^SAMPLE = '):
            sample = line.split(' = ', 1)[1]
            samples[sample] = {}
        if sample and line.startswith(('!Sample_characteristics_ch1 = ', '!Sample_source_name_ch1 = ', '!Sample_data_processing = ')):
            key, value = line.split(' = ', 1)
            samples[sample].setdefault(key, []).append(value)
    values, columns, offsets, obs, features, members, mapping_samples = [], [], [0], [], None, [], []
    with tarfile.open(a.archive) as archive:
        for donor in range(1, 5):
            accession = f'GSM{2230756 + donor}'
            name = f'{accession}_human{donor}_umifm_counts.csv.gz'
            member = archive.getmember(name)
            assert member.isfile()
            encoded = archive.extractfile(member).read()
            characteristics = samples[accession]['!Sample_characteristics_ch1']
            disease = [v.split(': ', 1)[1] for v in characteristics if v.startswith('type 2 diabetes mellitus: ')]
            assert disease == (['Yes'] if donor == 4 else ['No']), 'Disease provenance differs'
            assert any('No normalization is performed.' in v for v in samples[accession]['!Sample_data_processing'])
            condition = 'T2D' if disease == ['Yes'] else 'nonT2D'
            source_name = samples[accession]['!Sample_source_name_ch1']
            assert len(source_name) == 1 and 'donor' in source_name[0].lower()
            sample_id = f'human{donor}'
            mapping_samples.append(dict(id=sample_id, biologicalReplicateID=sample_id, donorID=sample_id,
                                        condition=condition, batchID='unreported', organism='NCBITaxon:9606'))
            start = len(obs)
            with gzip.GzipFile(fileobj=io.BytesIO(encoded)) as stream:
                reader = csv.reader(io.TextIOWrapper(stream))
                header = next(reader)
                assert header[:3] == ['', 'barcode', 'assigned_cluster']
                if features is None:
                    features = header[3:]
                    assert len(features) == len(set(features))
                assert header[3:] == features, 'Feature identity or order differs; never intersect silently'
                for row in reader:
                    assert len(row) == len(header)
                    counts = np.asarray(row[3:], dtype=np.float64)
                    assert np.isfinite(counts).all() and (counts >= 0).all()
                    assert (counts == np.floor(counts)).all() and (counts <= np.iinfo(np.int32).max).all()
                    keep = np.flatnonzero(counts)
                    values.append(counts[keep].astype(np.int32))
                    columns.append(keep.astype(np.int32))
                    offsets.append(offsets[-1] + len(keep))
                    obs.append(dict(cell_id=row[0], source_barcode=row[1], cell_type=row[2],
                                    native_sample=sample_id, donor=sample_id, condition=condition))
            members.append(dict(name=name, sha256=hashlib.sha256(encoded).hexdigest(), cells=len(obs)-start,
                                accession=accession, condition=condition))
    observations = pd.DataFrame(obs).set_index('cell_id')
    assert observations.index.is_unique and len(observations) == 8569 and len(features) == 20125
    counts = sparse.csr_matrix((np.concatenate(values), np.concatenate(columns), np.asarray(offsets, dtype=np.int64)),
                               shape=(len(obs), len(features)))
    counts.check_format(full_check=True)
    assert counts.has_canonical_format
    data = ad.AnnData(X=counts, obs=observations, var=pd.DataFrame(index=pd.Index(features, name='gene')))
    data.uns['numivivo_source'] = dict(archiveSHA256=ARCHIVE_SHA, metadataSHA256=SOFT_SHA,
        selection='All four human donor matrices; all source human cells and genes. Separate mouse matrices excluded by organism.',
        cellTypeProvenance='Author assigned_cluster column, retained verbatim; not independent ground truth')
    prepared = a.out/'prepared.h5ad'
    data.write_h5ad(prepared, compression='gzip')
    mapping = dict(schemaVersion=1, id='baron2016-all-human', evidence='measured', countUnit='umiCount', matrixPath='X',
        sourceDescription='Baron 2016 GSE84133 original unnormalized UMI counts, all four human donors and source genes; author cell-type assignments retained. Human4 is T2D, human1-3 nonT2D. Disease is confounded with donor; no disease DE qualification.',
        sampleColumn='native_sample', groupColumn='cell_type', mitochondrialFeatureIDs=[], samples=mapping_samples)
    (a.out/'mapping.json').write_text(json.dumps(mapping, indent=2)+'\n')
    (a.out/'stream-plan.json').write_text(json.dumps(dict(schemaVersion=1, mapping=mapping, contrasts=[]), indent=2)+'\n')
    design = pd.crosstab(observations.cell_type, observations.donor)
    summary = dict(status='count-source-qualified-analysis-incomplete', archiveURL=BASE+'suppl/GSE84133_RAW.tar',
        archiveSHA256=ARCHIVE_SHA, metadataURL=BASE+'soft/GSE84133_family.soft.gz', metadataSHA256=SOFT_SHA,
        preparedSHA256=digest(prepared), preparedBytes=prepared.stat().st_size, cells=counts.shape[0], features=counts.shape[1],
        nonzeros=counts.nnz, totalUMIs=int(counts.sum()), zeroCountCells=int(np.sum(np.asarray(counts.sum(axis=1)).ravel() == 0)),
        members=members, cellTypesByDonor={str(t):{str(d):int(design.loc[t,d]) for d in design.columns} for t in design.index},
        residentNonzeroLimit=5_000_000, residentEligible=counts.nnz <= 5_000_000,
        integrationEligibility='Current native donor integration rejects disconnected donor/condition design; the single T2D donor cannot establish a replicated disease effect. Do not relabel disease as unreported to bypass this gate.',
        remaining='Full-cohort streamed reduction and biologically justified integration evaluation; no analysis or competitiveness claim from preparation alone')
    (a.out/'preparation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
