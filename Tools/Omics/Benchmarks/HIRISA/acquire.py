#!/usr/bin/env python3
"""Acquire the complete GEO HIRISA release and freeze its metadata-only design."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306664/'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as src:
        for block in iter(lambda: src.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def parse_samples(path):
    samples = []
    current = None
    for line in gzip.open(path, 'rt'):
        line = line.rstrip('\n')
        if line.startswith('^SAMPLE = '):
            current = {'accession': line.split(' = ', 1)[1]}
            samples.append(current)
        elif current is not None and line.startswith('!Sample_'):
            key, value = line[8:].split(' = ', 1)
            if key == 'characteristics_ch1':
                key, value = value.split(': ', 1)
            current.setdefault(key, []).append(value)
    return samples

def one(sample, key):
    values = sample[key]
    assert len(values) == 1, (sample['accession'], key, values)
    return values[0]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    args = parser.parse_args()
    root = args.root
    root.mkdir(parents=True, exist_ok=True)
    raw = root / 'raw'
    raw.mkdir(exist_ok=True)
    protocol_hash = digest(Path(__file__).with_name('PROTOCOL.md'))
    for name, relative in [('series.soft.gz', 'soft/GSE306664_family.soft.gz'),
                           ('filelist.txt', 'suppl/filelist.txt')]:
        path = root / name
        if not path.exists():
            with urllib.request.urlopen(BASE + relative, timeout=60) as response:
                path.write_bytes(response.read())
    sizes = {}
    for line in (root / 'filelist.txt').read_text().splitlines():
        fields = line.split('\t')
        if fields[0] == 'File' and fields[1].endswith('.h5'):
            sizes[fields[1]] = int(fields[3])
    samples = parse_samples(root / 'series.soft.gz')
    assert len(samples) == len(sizes) == 131
    assert len({s['accession'] for s in samples}) == 131
    donors = sorted({one(s, 'subject id') for s in samples})
    assert donors == ['2616BW', '3283BW', '3491BW', '3955BW', '6811BW']
    enriched = [s for s in samples if one(s, 'cell type') != 'PBMC']
    assert len(enriched) == 105
    comparisons = []
    unmatched = []
    for population in ['Bcell', 'Monocyte', 'NK', 'Tcell']:
        for condition in ['IFNa', 'IFNb', 'IFNg', 'IFN-L1']:
            pairs = []
            treated = [s for s in enriched if one(s, 'cell type') == population
                       and one(s, 'treatment') == condition]
            assert len(treated) == 5
            for sample in treated:
                controls = [c for c in enriched if one(c, 'treatment') == 'none'
                            and all(one(c, k) == one(sample, k)
                                    for k in ['subject id', 'cell type', 'batch id', 'pool id'])]
                assert len(controls) <= 1
                if not controls:
                    unmatched.append(sample['accession'])
                    continue
                pairs.append({'donor': one(sample, 'subject id'),
                              'treated': sample['accession'], 'control': controls[0]['accession'],
                              'batch': one(sample, 'batch id'), 'pool': one(sample, 'pool id')})
            comparisons.append({'population': population, 'treatment': condition, 'pairs': pairs})
    assert unmatched == ['GSM9205558']
    assert Counter(len(c['pairs']) for c in comparisons) == {5: 15, 4: 1}
    design = {'protocolSHA256': protocol_hash, 'samples': samples,
              'sourceMetadataSHA256': digest(root / 'series.soft.gz'),
              'fileListingSHA256': digest(root / 'filelist.txt'),
              'comparisons': comparisons, 'unmatchedTreatedSamples': unmatched,
              'assignmentUsesExpression': False}
    design_path = root / 'design.json'
    encoded = json.dumps(design, sort_keys=True, indent=2) + '\n'
    if design_path.exists():
        assert design_path.read_text() == encoded, 'Frozen metadata design changed'
    else:
        design_path.write_text(encoded)
    existing = {}
    receipt_path = root / 'source-receipt.json'
    if receipt_path.exists():
        existing = {s['accession']: s for s in json.loads(receipt_path.read_text())['files']}

    def fetch(sample):
        url = one(sample, 'supplementary_file_1').replace('ftp://', 'https://', 1)
        name = url.rsplit('/', 1)[1]
        assert name in sizes
        destination = raw / name
        initial = root / name
        if initial.exists() and not destination.exists():
            initial.rename(destination)
        if not destination.exists():
            temporary = destination.with_suffix('.partial')
            # Never silently replace a failed partial transfer.
            assert not temporary.exists(), str(temporary)
            with urllib.request.urlopen(url, timeout=90) as response, temporary.open('xb') as out:
                for block in iter(lambda: response.read(1024 * 1024), b''):
                    out.write(block)
            assert temporary.stat().st_size == sizes[name], name
            temporary.rename(destination)
        assert destination.stat().st_size == sizes[name], name
        sha = digest(destination)
        if sample['accession'] in existing:
            assert existing[sample['accession']]['sha256'] == sha
        entry = {'accession': sample['accession'], 'path': str(destination.relative_to(root)),
                 'url': url, 'bytes': sizes[name], 'sha256': sha}
        print(json.dumps(entry), flush=True)
        return entry

    with ThreadPoolExecutor(max_workers=4) as pool:
        files = list(pool.map(fetch, samples))
    assert {Path(f['path']).name for f in files} == set(sizes)
    receipt = {'protocolSHA256': protocol_hash, 'designSHA256': digest(design_path),
               'files': files, 'totalBytes': sum(f['bytes'] for f in files), 'complete': True}
    receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'files': len(files), 'totalBytes': receipt['totalBytes'],
                      'pairedComparisons': len(comparisons), 'unmatched': unmatched}), flush=True)

if __name__ == '__main__':
    main()
