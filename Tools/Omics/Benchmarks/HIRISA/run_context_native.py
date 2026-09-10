#!/usr/bin/env python3
"""Run each complete frozen lineage aggregation, independent check and replay once."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('study', 'repo', 'runtime', 'python', 'hdf5'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    root = args.study/'context-transfer'
    frozen = json.loads((root/'freeze.json').read_text())
    reference = json.loads((root/'reference/check.json').read_text())
    runtime = json.loads((args.runtime/'environment.json').read_text())
    binary = args.runtime/'numivivo'
    assert reference['status'] == 'passed' and reference['freezeSHA256'] == sha(root/'freeze.json')
    assert sha(binary) == runtime['binarySHA256'] == '70a5bdcb258f35f947811bb7b5f2ddc777dfef19ce26880a74de818c8586ac95'
    assert all(sha(args.repo/name) == digest for name, digest in runtime['sourceFiles'].items())
    for name, identity in frozen['files'].items():
        assert sha(root/name) == identity['sha256']
    assert sha(args.study/'native-release/original.h5ad') == frozen['sourceSHA256']
    assert sha(args.hdf5) == 'a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e'
    checker = Path(__file__).with_name('check_context_native.py')
    write(root/'native-execution-freeze.json', dict(startedAtUnix=time.time(), coordinatorPID=os.getpid(),
          binarySHA256=sha(binary), runtimeEnvironmentSHA256=sha(args.runtime/'environment.json'),
          sourceSHA256=frozen['sourceSHA256'], frozenMetadataSHA256=sha(root/'freeze.json'),
          referenceCheckSHA256=sha(root/'reference/check.json'), runnerSHA256=sha(Path(__file__)),
          checkerSHA256=sha(checker), minimumAvailableBytes=3*2**30,
          psBefore=subprocess.check_output(['ps', '-axo', 'pid,etime,pcpu,rss,comm'], text=True),
          scope='Count aggregation and replay only; no response-model fit or scoring'))
    env = {**os.environ, 'NUMIVIVO_HDF5_LIBRARY': str(args.hdf5),
           'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}
    commands = []
    for lineage in ('B', 'Mono', 'NK', 'CD4-T', 'CD8-T', 'other-T'):
        bundle = root/(lineage+'-native')
        phases = [
            ('publish', ['/usr/bin/time', '-l', str(binary), 'singlecell-h5ad-pseudobulk',
                         str(args.study/'native-release/original.h5ad'), '--plan', str(root/(lineage+'-plan.json')), '--output', str(bundle)]),
            ('independent', [str(args.python), str(checker), '--bundle', str(bundle), '--frozen', str(root),
                             '--reference', str(root/'reference'), '--lineage', lineage, '--out', str(root/(lineage+'-independent.json'))]),
            ('verify', ['/usr/bin/time', '-l', str(binary), 'singlecell-h5ad-pseudobulk-verify', str(bundle)])]
        for phase, command in phases:
            name = lineage+'-'+phase
            free = shutil.disk_usage(root).free
            assert free >= 3*2**30, (name, free, 'Headroom; retain existing work, do not relaunch completed phases')
            started = time.time()
            write(root/(name+'-start.json'), dict(command=command, startedAtUnix=started, freeBytes=free))
            with (root/(name+'.log')).open('x') as stream:
                process = subprocess.Popen(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
                write(root/(name+'-process.json'), dict(pid=process.pid, coordinatorPID=os.getpid(), command=command))
                return_code = process.wait()
            status = dict(phase=name, returnCode=return_code, seconds=time.time()-started, freeBytesAfter=shutil.disk_usage(root).free)
            write(root/(name+'-status.json'), status)
            commands.append(status)
            print(json.dumps(status), flush=True)
            assert return_code == 0, name
        write(root/(lineage+'-complete.json'), dict(status='passed', lineage=lineage,
              reportSHA256=sha(bundle/'report.json'), receiptSHA256=sha(bundle/'receipt.json'),
              independentCheckSHA256=sha(root/(lineage+'-independent.json')), nativeReplayPassed=True))
    write(root/'native-complete.json', dict(status='passed', commands=commands,
          selectedCells=frozen['selectedCells'], originalSourceCells=frozen['sourceCells'],
          allSixLineagesAndReplaysPassed=True, responseFittingPerformed=False))


if __name__ == '__main__':
    main()
