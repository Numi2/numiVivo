#!/usr/bin/env python3
"""H5MU peak-mapping structural controls only; no real ATAC qualification."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--fixture', type=Path, default=Path(__file__).parent / 'evidence/2026-09-09/structural/valid-bundle')
p.add_argument('--out', type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=False)
dataset = json.loads((a.fixture / 'dataset.json').read_text())
plan = dict(schemaVersion=1, id=dataset['id'], evidence='synthetic', sourceDescription=dataset['sourceDescription'],
            samples=dataset['samples'], sampleColumn='sample', barcodeColumn='barcode',
            defaultObservationKind='cell', observationKindColumn='observation_kind', assays=[])
for assay in dataset['assays']:
    entry = {k: assay[k] for k in ['id', 'kind', 'featureNamespace', 'countUnit', 'genomeAssembly']}
    entry.update(sourceName=assay['id'], matrixPath='X', featureNameColumn='name')
    if assay['kind'] == 'chromatinAccessibility':
        entry['peakIDConvention'] = 'contig:start-end:zero-based-half-open'
    plan['assays'].append(entry)
checks = []
for name, mutate in [
    ('valid', lambda d: None),
    ('missing-convention', lambda d: d['assays'][1].pop('peakIDConvention')),
    ('wrong-convention', lambda d: d['assays'][1].update(peakIDConvention='one-based-inclusive')),
    ('missing-assembly', lambda d: d['assays'][1].pop('genomeAssembly')),
    ('wrong-unit', lambda d: d['assays'][1].update(countUnit='umiCount'))]:
    selected = copy.deepcopy(plan)
    mutate(selected)
    path = a.out / (name + '-plan.json')
    path.write_text(json.dumps(selected, indent=2) + '\n')
    dest = a.out / name
    args = [str(a.binary), 'multiassay-h5mu-import', str(a.fixture / 'dataset.h5mu'),
            '--plan', str(path), '--output', str(dest)]
    r = subprocess.run(args, capture_output=True, text=True)
    (a.out / (name + '.log')).write_text(r.stdout + r.stderr)
    assert (r.returncode == 0) == (name == 'valid'), (name, r.stderr)
    if name == 'valid':
        assert json.loads((dest / 'dataset.json').read_text()) == dataset
        verify = subprocess.run([str(a.binary), 'multiassay-verify', str(dest)], capture_output=True, text=True)
        (a.out / 'verify.log').write_text(verify.stdout + verify.stderr)
        assert verify.returncode == 0, verify.stderr
    else:
        assert not dest.exists()
    checks.append(dict(name=name, arguments=args, exitCode=r.returncode, expectedSuccess=name == 'valid'))
def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
(a.out / 'checks.json').write_text(json.dumps(dict(status='passed', scope='synthetic peak mapping only',
    checks=checks, binarySHA256=sha(a.binary), checkerSHA256=sha(Path(__file__)),
    sourceSHA256=sha(a.fixture / 'dataset.h5mu')), indent=2) + '\n')
