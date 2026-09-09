#!/usr/bin/env python3
"""Prepare public HCC1395 chr6/chr17 RNA, not synthetic or mocked observations.

This is a bounded characterization experiment, not clinical qualification.
The complete published chr6/chr17 RNA archive is consumed and hashed. Genomic
coordinates are preserved. Only the CUL9/CAPN11 published variant loci are used
for splice characterization; a separate real missense SNV exercises pVACseq.
"""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import shutil
import tarfile
import urllib.request
import zipfile

COURSE = 'e8219708889092839d60ea13af31b002812d47f5'
COURSE_SHA = '5454d435af43804bee111dc085c599a3e223445cedc1c52768c6e98e981e37eb'
RNA_URL = 'https://genomedata.org/hcc1395/fastqs/chr6_and_chr17/RNAseq_Tumor.tar'
REF_URL = 'https://genomedata.org/hcc1395/references/genome/chr6_and_chr17/ref_genome.tar'
INDEX_URL = 'https://genomedata.org/hcc1395/references/transcriptome/chr6_and_chr17/ref_genome_hisat_index.tar.gz'
GTF_URL = 'https://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz'
SPLICE_POSITIONS = {('chr6', 43187096), ('chr6', 44181251)}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda: f.read(8 * 1024 * 1024), b''): h.update(data)
    return h.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


class MeasuredReader:
    def __init__(self, response, maximum: int):
        self.response = response
        self.maximum = maximum
        self.size = 0
        self.sha = hashlib.sha256()
    def read(self, size=-1):
        if size < 0: size = 8 * 1024 * 1024
        data = self.response.read(size)
        self.size += len(data)
        if self.size > self.maximum: raise ValueError('Public source exceeds declared download bound')
        self.sha.update(data)
        return data


def download(url: str, path: Path, maximum: int) -> dict:
    request = urllib.request.Request(url, headers={'User-Agent': 'NumiVivo-HCC1395-qualification'})
    with urllib.request.urlopen(request, timeout=120) as r, path.open('xb') as out:
        reader = MeasuredReader(r, maximum)
        for block in iter(lambda: reader.read(8 * 1024 * 1024), b''): out.write(block)
        record = {'url': url, 'sha256': reader.sha.hexdigest(), 'bytes': reader.size,
                  'etag': r.headers.get('ETag'), 'lastModified': r.headers.get('Last-Modified')}
    print('DOWNLOADED', json.dumps(record), flush=True)
    return record


def extract_regular(archive: Path, output: Path) -> list[Path]:
    output.mkdir()
    files = []
    total = 0
    with tarfile.open(archive, 'r:*') as tar:
        for member in tar:
            if member.isdir(): continue
            total += member.size
            if not member.isfile() or member.size > 2 * 1024**3 or total > 3 * 1024**3:
                raise ValueError('Unexpected archive member or size')
            target = output / Path(member.name).name
            if target.exists(): raise ValueError('Ambiguous flattened archive members')
            with tar.extractfile(member) as src, target.open('xb') as dest:
                shutil.copyfileobj(src, dest, 8 * 1024 * 1024)
            files.append(target)
    return files


def rna_archive(output: Path) -> tuple[dict, list[Path]]:
    output.mkdir()
    files = []
    request = urllib.request.Request(RNA_URL, headers={'User-Agent': 'NumiVivo-HCC1395-qualification'})
    with urllib.request.urlopen(request, timeout=120) as response:
        reader = MeasuredReader(response, 8 * 1024**3)
        expected = int(response.headers.get('Content-Length', '0'))
        with tarfile.open(fileobj=reader, mode='r|', bufsize=8 * 1024 * 1024) as tar:
            for member in tar:
                if member.isdir(): continue
                name = Path(member.name).name
                if name.startswith('._') or name == '.DS_Store': continue
                if not member.isfile() or member.size > 8 * 1024**3:
                    raise ValueError('Unexpected RNA archive member')
                if not re.search(r'\.(fastq|fq)(\.gz)?$', name):
                    print('UNUSED_RNA_ARCHIVE_MEMBER', name, member.size, flush=True)
                    continue
                target = output / (name if name.endswith('.gz') else name + '.gz')
                if target.exists(): raise ValueError('Duplicate RNA file name')
                with tar.extractfile(member) as src, target.open('xb') as dest:
                    if name.endswith('.gz'):
                        shutil.copyfileobj(src, dest, 8 * 1024 * 1024)
                    else:
                        with gzip.GzipFile(filename='', mode='wb', fileobj=dest, compresslevel=1, mtime=0) as gz:
                            shutil.copyfileobj(src, gz, 8 * 1024 * 1024)
                files.append(target)
                print('RNA_MEMBER', json.dumps({'name': name, 'originalBytes': member.size,
                      'savedBytes': target.stat().st_size, 'sha256': digest(target)}), flush=True)
        # Include padding/trailing bytes in the original archive identity.
        for _ in iter(lambda: reader.read(8 * 1024 * 1024), b''): pass
        if expected and reader.size != expected: raise ValueError('RNA archive download was incomplete')
        source = {'url': RNA_URL, 'bytes': reader.size, 'sha256': reader.sha.hexdigest(),
                  'etag': response.headers.get('ETag'), 'scope': 'all supplied chr6/chr17 RNA reads, no downsampling'}
    if not files: raise ValueError('No real RNA FASTQ files were found')
    return source, files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    evidence = root / 'evidence'; evidence.mkdir()
    sources = {}
    course = root / 'course.zip'
    sources['course'] = download(f'https://raw.githubusercontent.com/griffithlab/pVACtools_Intro_Course/{COURSE}/HCC1395_inputs.zip', course, 40 * 1024**2)
    if sources['course']['sha256'] != COURSE_SHA: raise ValueError('Pinned HCC1395 course archive changed')
    with zipfile.ZipFile(course) as z:
        prefix = 'HCC1395_inputs/'
        vcf = gzip.decompress(z.read(prefix + 'annotated.expression.vcf.gz')).decode()
        (root / 'source.vcf').write_text(vcf)
        (evidence / 'optitype_normal_result.tsv').write_bytes(z.read(prefix + 'optitype_normal_result.tsv'))
        (evidence / 'reference-junctions.tsv').write_bytes(z.read(prefix + 'HCC1395.splice_junctions.tsv'))
    header = [l for l in vcf.splitlines() if l.startswith('#')]
    records = [l for l in vcf.splitlines() if l and not l.startswith('#')]
    splice = [l for l in records if (l.split('\t')[0], int(l.split('\t')[1])) in SPLICE_POSITIONS]
    if len(splice) != 2: raise ValueError('The two predeclared public splice SNVs were not found exactly once')
    missense = next((l for l in records if l.startswith('chr6\t') and 'missense_variant' in l
                     and len(l.split('\t')[3]) == len(l.split('\t')[4]) == 1 and l.split('\t')[6] == 'PASS'), None)
    if missense is None: raise ValueError('No actual chr6 missense SNV available for parent report')
    (root / 'splice.vcf').write_text('\n'.join(header + splice) + '\n')
    (root / 'parent.vcf').write_text('\n'.join(header + [missense]) + '\n')
    selected = [{'chromosome': l.split('\t')[0], 'position1': int(l.split('\t')[1]),
                 'reference': l.split('\t')[3], 'alternate': l.split('\t')[4]} for l in splice + [missense]]
    write_json(evidence / 'selected-variants.json', selected)
    print('PREDECLARED_LOCI', json.dumps(selected), flush=True)
    reference_tar = root / 'reference.tar'
    sources['referenceArchive'] = download(REF_URL, reference_tar, 300 * 1024**2)
    ref_files = extract_regular(reference_tar, root / 'reference-source')
    refs = [p for p in ref_files if re.search(r'\.(fa|fasta)(\.gz)?$', p.name)]
    if len(refs) != 1: raise ValueError('Expected one original chromosome-reference FASTA')
    opener = gzip.open if refs[0].name.endswith('.gz') else open
    with opener(refs[0], 'rb') as src, (root / 'reference.fa').open('xb') as dst:
        shutil.copyfileobj(src, dst, 8 * 1024 * 1024)
    index_tar = root / 'hisat-index.tar.gz'
    sources['hisatIndex'] = download(INDEX_URL, index_tar, 500 * 1024**2)
    index_files = extract_regular(index_tar, root / 'hisat-index')
    first = [p for p in index_files if p.name.endswith('.1.ht2') or p.name.endswith('.1.ht2l')]
    if len(first) != 1: raise ValueError('Expected exactly one HISAT2 index')
    index_prefix = re.sub(r'\.1\.ht2l?$', '', str(first[0]))
    gtf_gz = root / 'ensembl105.gtf.gz'
    sources['ensembl105'] = download(GTF_URL, gtf_gz, 150 * 1024**2)
    count = 0
    with gzip.open(gtf_gz, 'rt') as src, (root / 'annotation.gtf').open('w') as dst:
        for line in src:
            if line.startswith('#'): dst.write(line); continue
            chrom, rest = line.split('\t', 1)
            if chrom in ('6', '17'):
                dst.write('chr' + chrom + '\t' + rest); count += 1
    if count < 1000: raise ValueError('Reference annotation is unexpectedly empty')
    sources['annotationTransform'] = {'operation': 'retain chromosomes 6 and 17; explicit chr-prefix rename only',
                                     'coordinatesChanged': False, 'rows': count, 'sha256': digest(root / 'annotation.gtf')}
    sources['rnaArchive'], fastqs = rna_archive(root / 'rna')
    pairs = {}
    for path in fastqs:
        match = re.fullmatch(r'(.+?)(?:_R|_)([12])(?:_001)?\.(?:fastq|fq)\.gz', path.name)
        if match is None: raise ValueError('Unrecognized mate filename: ' + path.name)
        stem, mate = match.groups()
        if mate in pairs.setdefault(stem, {}): raise ValueError('Duplicate read mate')
        pairs[stem][mate] = str(path)
    if not pairs or any(set(p) != {'1', '2'} for p in pairs.values()): raise ValueError('Incomplete RNA read pairs')
    plan = {'referenceFASTA': str(root / 'reference.fa'), 'annotationGTF': str(root / 'annotation.gtf'),
            'hisatIndex': index_prefix, 'read1': [pairs[k]['1'] for k in sorted(pairs)],
            'read2': [pairs[k]['2'] for k in sorted(pairs)], 'strand': 'RF',
            'scope': 'published chr6/chr17 RNA subset, not a whole-genome clinical study',
            'spliceCharacterization': {'junctionMinimumReads': 1, 'variantDistance': 2000,
                                       'epitopeLengths': [9], 'alleles': ['HLA-A*29:02'],
                                       'purpose': 'software/reference characterization; not candidate qualification'}}
    write_json(root / 'alignment-plan.json', plan)
    write_json(evidence / 'source-manifest.json', sources)
    print('REAL_INPUT_PREPARATION_COMPLETE', json.dumps({'referenceSHA256': digest(root / 'reference.fa'),
          'annotationSHA256': digest(root / 'annotation.gtf'), 'rnaArchiveSHA256': sources['rnaArchive']['sha256'],
          'matePairs': len(pairs)}), flush=True)


if __name__ == '__main__': main()
