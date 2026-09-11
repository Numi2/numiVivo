"""Execute frozen native predictions; deduplicate verified artifacts before scoring."""
from pathlib import Path
import gzip, hashlib, json, os, shutil, subprocess, time

root = Path.home() / 'numivivo-gse181897-20260911'
src = root / 'prediction-inputs'
out = root / 'prediction-execution'
assert out.is_dir(), 'Resume requires the preserved failed run'

binary = Path.home() / 'numivivo-donor-intervals-20260911/runtime/numivivo-omics'
library = Path.home() / 'numivivo-multiassay-hdf5-20260909/libhdf5.dylib'

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while b := f.read(1048576): h.update(b)
    return h.hexdigest()

def write(p, d):
    temp = p.with_suffix(p.suffix + '.tmp')
    temp.write_text(json.dumps(d, indent=2, sort_keys=True, allow_nan=False) + '\n')
    os.replace(temp, p)

freeze = json.loads((src / 'input-freeze.json').read_bytes())
for name, digest in freeze['files'].items(): assert sha(src / name) == digest, name
assert freeze['primaryConditionMeaningQualified'] and not freeze['fitStarted']
assert sha(binary) == '0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b'
env = dict(os.environ, NUMIVIVO_HDF5_LIBRARY=str(library))
commands = json.loads((out / 'commands.json').read_bytes())
old_index = json.loads((out / 'retention.json').read_bytes())
bundles, objects = old_index['bundles'], old_index['objects']
old_state = json.loads((out / 'pipeline.json').read_bytes())
assert old_state['status'] == 'failed-native-command' and old_state['activeNativePID'] is None
assert subprocess.run(['ps', '-p', str(old_state['pid']), '-o', 'command='], capture_output=True).returncode == 1
assert old_state['inputFreezeSHA256'] == sha(src / 'input-freeze.json')
for digest, item in objects.items():
    assert sha(out / 'objects' / (digest + '.gz')) == item['compressedSHA256']
shutil.copyfile(out / 'pipeline.json', out / 'failed-pipeline.json')
shutil.copyfile(root / 'run_prediction.py', out / 'initial-runner.py')
shutil.copyfile(__file__, out / 'resume-runner.py')
state = dict(status='running-native-prediction', pid=os.getpid(), startedUnix=time.time(),
    binary=str(binary), binarySHA256=sha(binary), hdf5LibrarySHA256=sha(library),
    runnerSHA256=sha(__file__), inputFreezeSHA256=sha(src / 'input-freeze.json'),
    scoringStarted=False, completedOriginChunks=old_state['completedOriginChunks'],
    resumedFromRunnerSHA256=old_state['runnerSHA256'], priorStartedUnix=old_state['startedUnix'])
write(out / 'pipeline.json', state)

def run(args, name):
    start = time.time()
    with (out / 'logs' / (name + '.log')).open('xb') as log:
        p = subprocess.Popen([str(binary), *map(str, args)], env=env, stdout=log, stderr=subprocess.STDOUT)
        state['activeNativePID'] = p.pid; state['activeCommand'] = name
        write(out / 'pipeline.json', state)
        code = p.wait()
    commands.append(dict(name=name, arguments=list(map(str, args)), nativePID=p.pid,
        exitCode=code, seconds=time.time()-start))
    write(out / 'commands.json', commands)
    state['activeNativePID'] = None
    if code:
        state['status'] = 'failed-native-command'; write(out / 'pipeline.json', state)
        raise RuntimeError(name + ': ' + (out / 'logs' / (name + '.log')).read_text()[-2000:])

def retain(folder, key, remove=True):
    assert folder.parent == out and folder.is_dir() and not folder.is_symlink()
    entries = list(folder.rglob('*'))
    assert not any(p.is_symlink() for p in entries)
    files = sorted(p for p in entries if p.is_file())
    assert files and (folder / 'receipt.json').is_file()
    members = {}
    for p in files:
        digest = sha(p)
        path = out / 'objects' / (digest + '.gz')
        if digest not in objects:
            assert not path.exists()
            with p.open('rb') as reader, path.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', mtime=0, filename='', compresslevel=6) as writer:
                shutil.copyfileobj(reader, writer, 1048576)
            h = hashlib.sha256(); size = 0
            with gzip.open(path, 'rb') as f:
                while b := f.read(1048576): h.update(b); size += len(b)
            assert h.hexdigest() == digest and size == p.stat().st_size
            objects[digest] = dict(bytes=size, compressedBytes=path.stat().st_size, compressedSHA256=sha(path))
        assert objects[digest]['bytes'] == p.stat().st_size
        members[str(p.relative_to(folder))] = dict(SHA256=digest, bytes=p.stat().st_size)
    bundles[key] = members
    write(out / 'retention.json', dict(bundles=bundles, objects=objects))
    if remove:
        handles = subprocess.run(['/usr/sbin/lsof', '+D', str(folder)], capture_output=True, text=True)
        assert handles.returncode == 1 and not handles.stdout.strip() and not handles.stderr.strip()
        # Exact invocation-owned, terminal native bundle only; every file now has
        # a fully decoded and hash-verified retained object and logical mapping.
        shutil.rmtree(folder)

chunks = json.loads((src / 'chunks.json').read_bytes())
for origin in ['Kang', 'HIRISA']:
    model = out / ('model-' + origin)
    if not model.exists():
        run(['singlecell-perturbation-fit', src / origin / 'training.h5ad', '--plan', src / origin / 'training.json', '--output', model], origin + '-fit')
        run(['singlecell-perturbation-verify', model], origin + '-model-verify')
    else:
        assert origin == 'Kang' and 'model-Kang' in bundles
        for rel, item in bundles['model-Kang'].items(): assert sha(model / rel) == item['SHA256']
    retain(model, 'model-' + origin, remove=False)
    for chunk in chunks:
        key = origin + '-' + chunk['id']
        if 'prediction-' + key in bundles: continue
        dest = out / ('prediction-' + key)
        args = ['singlecell-perturbation-predict', src / 'queries' / (chunk['id'] + '.h5ad'),
                '--plan', src / 'queries' / (chunk['id'] + '.json'), '--reference', model, '--output', dest]
        if dest.exists():
            assert key == 'Kang-14'
            assert any(c['name'] == key + '-predict' and c['exitCode'] == 0 for c in commands)
            run(['singlecell-perturbation-prediction-verify', dest], key + '-verify-after-disk-recovery')
        else:
            run(args, key + '-predict')
            run(['singlecell-perturbation-prediction-verify', dest], key + '-verify')
        retain(dest, 'prediction-' + key)
        state['completedOriginChunks'] += 1
        state['retainedCompressedBytes'] = sum(x['compressedBytes'] for x in objects.values())
        write(out / 'pipeline.json', state)
        print(json.dumps({k:state[k] for k in ['completedOriginChunks','retainedCompressedBytes']}), flush=True)
    retain(model, 'model-' + origin)

assert state['completedOriginChunks'] == 32
files = {str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'pipeline.json'}
write(out / 'prediction-freeze.json', dict(schemaVersion=1, inputFreezeSHA256=sha(src / 'input-freeze.json'),
    files=files, donorPredictions=124, originChunks=32, nativeCommands=len(commands), scoringStarted=False,
    completedUnix=time.time(), retainedLogicalFiles=sum(len(x) for x in bundles.values()), uniqueObjects=len(objects)))
state.update(status='completed-native-predictions-and-reconstructions', activeCommand=None,
    finishedUnix=time.time(), predictionFreezeSHA256=sha(out / 'prediction-freeze.json'))
write(out / 'pipeline.json', state)
print(json.dumps(state), flush=True)
