"""Restore the dated complete Parse count evidence and check it without network I/O.

Optional APFS reuse copies only archive-verified bytes into separate inodes.
It never hard-links a live study, overwrites a destination, or reruns source I/O.
"""
import argparse
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import uuid


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--reuse-from', type=Path,
                        help='Optional existing study for verified APFS copy-on-write copies')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[4]
    archive_module = repo / 'Tools/Omics/PerturbationPrediction/Duration/archive.py'
    sys.path.insert(0, str(archive_module.parent))
    spec = importlib.util.spec_from_file_location('parse_evidence_archive', archive_module)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    evidence = args.evidence.resolve()
    outer = repo / 'Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py'
    subprocess.run([sys.executable, str(outer), str(evidence)], check=True)
    summary = json.loads((evidence / 'summary.json').read_text())
    assert summary['status'] == 'passed' and summary['donors'] == 12
    assert summary['allDonorNativeVerificationsPassed'] and not summary['predictionFitOrScoring']
    dependencies = json.loads((evidence / 'dependencies.json').read_text())['required']
    archives = []
    for dependency in dependencies:
        directory = repo / dependency['path']
        assert sha(directory / 'manifest.json') == dependency['manifestSHA256']
        subprocess.run([sys.executable, str(outer), str(directory)], check=True)
        archive = directory / dependency['archive']
        assert sha(archive) == dependency['archiveSHA256']
        members = module.verify(archive)
        if dependency['archiveSHA256'] == 'f08c78ca7a74ab01eab6c3f7748615f45e217cbf95446ba35ef830eaf9c75d02':
            selected = members
        elif dependency['archiveSHA256'] == 'fa98901c2f0ce867620f7cd37494f74fca4c238dd32159c4a7751e55cbe5abe4':
            selected = {'build-pooled/numivivo-omics': members['runtime/corrected/numivivo-omics']}
        else:
            raise ValueError('Unknown dependency for this dated qualification')
        archives.append((archive, sha(archive), selected))
    assert len(archives) == 2 and len({digest for _, digest, _ in archives}) == 2
    archive = evidence / 'results.tar.gz'
    archives.append((archive, sha(archive), module.verify(archive)))
    expected = {}
    for _, _, members in archives:
        for name, identity in members.items():
            assert name not in expected or expected[name] == identity, ('Conflicting archives', name)
            expected[name] = identity
    destination = args.out.absolute()
    assert not destination.exists() and not destination.is_symlink()
    staging = destination.with_name(destination.name + '.incomplete-' + uuid.uuid4().hex)
    staging.mkdir(parents=True, exist_ok=False)
    reused = extracted = 0
    clonefile = getattr(ctypes.CDLL(None, use_errno=True), 'clonefile', None)
    if clonefile is not None:
        clonefile.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
        clonefile.restype = ctypes.c_int
    for archive, _, members in archives:
        by_hash = {}
        for name, identity in members.items():
            by_hash.setdefault(identity['SHA256'], []).append((name, identity))
        # Consume each compressed archive sequentially. Random per-path seeks
        # would repeatedly decompress the archive for thousands of small files.
        with tarfile.open(archive, 'r|gz') as container:
            for item in container:
                if not item.name.startswith('objects/'):
                    continue
                digest = item.name.split('/')[1]
                previous = None
                for name, identity in by_hash.get(digest, []):
                    path = staging / name
                    if path.exists():
                        assert path.is_file() and not path.is_symlink() and sha(path) == digest
                        previous = path
                        continue
                    path.parent.mkdir(parents=True, exist_ok=True)
                    candidate = previous or (args.reuse_from / name if args.reuse_from else None)
                    cloned = False
                    if clonefile is not None and candidate and candidate.is_file() and not candidate.is_symlink() and sha(candidate) == digest:
                        cloned = clonefile(os.fsencode(candidate), os.fsencode(path), 0) == 0
                        if not cloned and path.exists():
                            path.unlink()
                    if cloned:
                        assert path.stat().st_ino != candidate.stat().st_ino
                        reused += 1
                    else:
                        assert shutil.disk_usage(staging).free > identity['bytes'] + 200000000
                        source = previous.open('rb') if previous else container.extractfile(item)
                        with source, path.open('xb') as output:
                            shutil.copyfileobj(source, output, 1048576)
                        extracted += 1
                    assert path.stat().st_size == identity['bytes'] and sha(path) == digest
                    path.chmod(identity['mode'])
                    previous = path
    recipe = staging / 'donor-recipe-v2'
    freeze = json.loads((recipe / 'source-freeze.json').read_text())
    for item in freeze['sourceFiles']:
        assert sha(recipe / item['path']) == item['SHA256']
    assert sha(staging / 'build-pooled/numivivo-omics') == freeze['binarySHA256']
    # The frozen checker compares saved complete arrays and receipts. It does
    # not fetch count ranges, execute the native binary, fit, or score a model.
    environment = dict(os.environ, NUMIVIVO_PARSE_STUDY=str(staging),
                       OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    offline_launcher = (
        "import sys,runpy; "
        "sys.addaudithook(lambda event,args: "
        "(_ for _ in ()).throw(RuntimeError('Offline checker attempted '+event)) "
        "if event in {'socket.connect','socket.getaddrinfo','subprocess.Popen','os.system','os.exec','os.posix_spawn'} else None); "
        "sys.path.insert(0,sys.argv[1]); runpy.run_path(sys.argv[1]+'/check_donors.py',run_name='__main__')"
    )
    with (staging / 'restored-offline-check.log').open('xb') as log:
        subprocess.run([sys.executable, '-c', offline_launcher, str(recipe)], env=environment,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    checked = json.loads((staging / 'donor-counts-verification.json').read_text())
    assert checked['status'] == 'passed' and not checked['predictionFitOrScoring']
    # check_donors regenerates two deterministic products; check them too.
    assert all(sha(staging / name) == identity['SHA256'] for name, identity in expected.items())
    assert all(sha(archive) == digest for archive, digest, _ in archives)
    result = dict(status='restored-and-offline-checked',files=len(expected),
                  copyOnWriteFiles=reused, extractedFiles=extracted,
                  cells=checked['cells'], donors=checked['donors'],
                  nativeReplayExecuted=False, sourceRangesFetched=False,
                  offlineCheckerNetworkAndChildProcessesForbidden=True,
                  predictionFitOrScoring=False,
                  archives=[dict(path=str(path), SHA256=digest) for path, digest, _ in archives])
    (staging / 'restoration.json').write_text(json.dumps(result, indent=2) + '\n')
    assert not destination.exists() and not destination.is_symlink()
    os.rename(staging, destination)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
