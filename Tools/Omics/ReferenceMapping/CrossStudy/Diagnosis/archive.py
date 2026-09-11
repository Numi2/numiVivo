#!/usr/bin/env python3
"""Retain complete diagnostics; restore native programs using exact prior inputs."""
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def pack(study, out):
    out.mkdir(parents=True, exist_ok=False)
    records, text_files, recipe, dependencies = [], [], [], []

    def store(name, raw, bundled=False):
        encoded = gzip.compress(raw, compresslevel=9, mtime=0)
        stored = name.replace('/', '--') + '.gz'
        (out/stored).write_bytes(encoded)
        records.append(dict(sourcePath=name, sourceBytes=len(raw), sourceSHA256=digest(raw),
                            storedPath=stored, storedBytes=len(encoded), storedSHA256=digest(encoded),
                            gzipEncoded=True, bundledSourceFiles=bundled))

    def text_record(path, label):
        raw = path.read_bytes()
        return dict(path=label, bytes=len(raw), SHA256=digest(raw), rawUTF8=raw.decode())

    selected = [p for p in study.iterdir() if p.is_file()]
    selected += list((study/'programs').rglob('*')) + list((study/'checks').rglob('*'))
    for path in sorted(selected):
        if not path.is_file():
            continue
        assert not path.is_symlink()
        name = str(path.relative_to(study))
        if path.name == 'original.h5ad':
            cohort = path.parent.name
            source_name = 'kang.h5ad' if cohort == 'kang' else 'ding-query.h5ad'
            raw = path.read_bytes()
            dependencies.append(dict(path=name, priorInputName=source_name, bytes=len(raw), SHA256=digest(raw)))
        elif path.suffix in {'.bin', '.npz'}:
            store(name, path.read_bytes())
        elif path.suffix in {'.json', '.log', '.py', '.md'}:
            text_files.append(text_record(path, 'study/'+name))
        else:
            raise ValueError('Unclassified evidence: '+name)
    for path in sorted(Path(__file__).resolve().parent.iterdir()):
        if path.is_file():
            recipe.append(text_record(path, 'recipe/'+path.name))
    for name, files in [('execution', text_files), ('recipe', recipe)]:
        store(name+'.json', (json.dumps(dict(schemaVersion=1, sourceFiles=files), sort_keys=True,
                                       separators=(',', ':'))+'\n').encode(), True)
    assert len(dependencies) == 2
    manifest = dict(schemaVersion=1, records=records, externalNativeSources=dependencies,
                    sourcePolicy='Exact complete prior H5AD inputs required for native replay; no source counts discarded.')
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(dict(status='archived', records=len(records), bundledFiles=len(text_files)+len(recipe),
                          storedBytes=sum(r['storedBytes'] for r in records))))


def restore(archive, prior_inputs, out):
    verifier = Path(__file__).resolve().parents[3]/'PerturbationPrediction/Replogle2020/verify_target_archive.py'
    subprocess.run([sys.executable, str(verifier), str(archive)], check=True)
    manifest = json.loads((archive/'manifest.json').read_text())
    sources = []
    for item in manifest['externalNativeSources']:
        source = prior_inputs/item['priorInputName']
        assert source.resolve().is_relative_to(prior_inputs.resolve())
        assert source.stat().st_size == item['bytes'] and digest(source.read_bytes()) == item['SHA256']
        sources.append((source, item))
    out.mkdir(parents=True, exist_ok=False)

    def emit(name, raw):
        if not name.startswith('programs/'):
            return
        path = out/name.removeprefix('programs/')
        assert path.resolve().is_relative_to(out.resolve()) and not path.exists()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)

    for item in manifest['records']:
        raw = gzip.decompress((archive/item['storedPath']).read_bytes())
        if item.get('bundledSourceFiles'):
            for source in json.loads(raw)['sourceFiles']:
                if source['path'].startswith('study/programs/'):
                    emit(source['path'].removeprefix('study/'), source['rawUTF8'].encode())
        else:
            emit(item['sourcePath'], raw)
    for source, item in sources:
        path = out/item['path'].removeprefix('programs/')
        assert path.resolve().is_relative_to(out.resolve()) and not path.exists()
        cloned = subprocess.run(['/bin/cp', '-c', str(source), str(path)], capture_output=True).returncode == 0
        if not cloned:
            assert not path.exists()
            shutil.copyfile(source, path)
        assert path.stat().st_size == item['bytes'] and digest(path.read_bytes()) == item['SHA256']
    print(json.dumps(dict(status='restored', nativePrograms=2, output=str(out))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    p = sub.add_parser('pack')
    p.add_argument('--study', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    r = sub.add_parser('restore')
    r.add_argument('--archive', type=Path, required=True)
    r.add_argument('--prior-inputs', type=Path, required=True)
    r.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'pack':
        pack(args.study, args.out)
    else:
        restore(args.archive, args.prior_inputs, args.out)
