#!/usr/bin/env python3
"""Wait on the existing count job, run frozen native folds, and seal all outputs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import shutil
import subprocess
import time

from prepare_context_predictions import sha, read, write

LINEAGES = ('B', 'Mono', 'NK', 'CD4-T', 'CD8-T', 'other-T')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('root', 'repo', 'runtime', 'python', 'hdf5'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--count-coordinator', type=int, required=True)
    args = parser.parse_args()
    root = args.root
    runtime = read(args.runtime/'environment.json')
    binary = args.runtime/'numivivo'
    assert sha(binary) == runtime['binarySHA256'] == '70a5bdcb258f35f947811bb7b5f2ddc777dfef19ce26880a74de818c8586ac95'
    assert all(sha(args.repo/name) == digest for name, digest in runtime['sourceFiles'].items())
    write(root/'prediction-coordinator-start.json', dict(pid=os.getpid(), countCoordinatorPID=args.count_coordinator,
          startedAtUnix=time.time(), runnerSHA256=sha(Path(__file__)), binarySHA256=sha(binary),
          preparerSHA256=sha(Path(__file__).with_name('prepare_context_predictions.py')),
          status='waiting-for-existing-count-job; never restart it'))
    while not (root/'native-complete.json').exists():
        try:
            os.kill(args.count_coordinator, 0)
        except ProcessLookupError:
            assert (root/'native-complete.json').exists(), 'Original count job ended without complete qualification; inspect its retained status, do not restart blindly'
        time.sleep(10)
    count_result = read(root/'native-complete.json')
    assert count_result['status'] == 'passed' and count_result['allSixLineagesAndReplaysPassed']
    for lineage in LINEAGES:
        inputs = root/'prediction-inputs'/lineage
        if not inputs.exists():
            command = [str(args.python), str(Path(__file__).with_name('prepare_context_predictions.py')),
                       '--root', str(root), '--lineage', lineage, '--out', str(inputs)]
            with (root/(lineage+'-prediction-prepare.log')).open('x') as stream:
                subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True)
        transport = read(inputs/'transport.json')
        assert transport['status'] == 'passed' and transport['scriptSHA256'] == sha(Path(__file__).with_name('prepare_context_predictions.py'))
        assert transport['nativeCompleteSHA256'] == sha(root/(lineage+'-complete.json'))
        assert all(sha(inputs/name) == digest for name, digest in transport['files'].items())
        assert not (root/(lineage+'-prediction')).exists()
    write(root/'prediction-execution-freeze.json', dict(startedAtUnix=time.time(), coordinatorPID=os.getpid(),
          nativeCompleteSHA256=sha(root/'native-complete.json'), frozenMetadataSHA256=sha(root/'freeze.json'),
          transportHashes={s: sha(root/'prediction-inputs'/s/'transport.json') for s in LINEAGES},
          binarySHA256=sha(binary), hdf5SHA256=sha(args.hdf5), maximumConcurrentLineages=2,
          minimumAvailableBytesBeforePhase=3*2**30, targetsOpenedForScoring=False,
          psBefore=subprocess.check_output(['ps', '-axo', 'pid,etime,pcpu,rss,comm'], text=True)))
    env = {**os.environ, 'NUMIVIVO_HDF5_LIBRARY': str(args.hdf5),
           'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'}

    def run(lineage):
        bundle = root/(lineage+'-prediction')
        inputs = root/'prediction-inputs'/lineage
        statuses = []
        for phase, native_args in [
            ('publish', ['singlecell-perturbation-batch', str(root/(lineage+'-native')), '--plan', str(inputs/'plan.json'), '--output', str(bundle)]),
            ('verify', ['singlecell-perturbation-batch-verify', str(bundle)])]:
            name = lineage+'-prediction-'+phase
            free = shutil.disk_usage(root).free
            assert free >= 3*2**30, (name, free, 'Headroom; retain outputs for a targeted continuation')
            command = ['/usr/bin/time', '-l', str(binary), *native_args]
            started = time.time()
            write(root/(name+'-start.json'), dict(command=command, startedAtUnix=started, freeBytes=free))
            with (root/(name+'.log')).open('x') as stream:
                process = subprocess.Popen(command, env=env, stdout=stream, stderr=subprocess.STDOUT)
                write(root/(name+'-process.json'), dict(pid=process.pid, coordinatorPID=os.getpid(), command=command))
                code = process.wait()
            status = dict(phase=name, returnCode=code, seconds=time.time()-started, freeBytesAfter=shutil.disk_usage(root).free)
            write(root/(name+'-status.json'), status)
            statuses.append(status)
            print(status, flush=True)
            assert code == 0, name
            if phase == 'publish':
                receipt = read(bundle/'receipt.json')
                checks = read(inputs/'projection-checks.json')
                assert len(receipt['folds']) == len(checks) == 20
                for actual, expected in zip(receipt['folds'], checks):
                    assert actual['id'] == expected['id']
                    if actual['status'] == 'completed':
                        assert actual['trainingAggregate']['bytes'] == list(bytes.fromhex(expected['trainingAggregateSHA256']))
                        assert actual['queryAggregate']['bytes'] == list(bytes.fromhex(expected['queryAggregateSHA256']))
                paths = [bundle/'receipt.json', bundle/'plan.json', *sorted((bundle/'folds').rglob('*.json'))]
                write(root/(lineage+'-prediction-freeze.json'), dict(createdAtUnix=time.time(), folds=20,
                      completed=sum(s['status'] == 'completed' for s in receipt['folds']), targetsOpenedForScoring=False,
                      files={str(p.relative_to(root)): sha(p) for p in paths}))
        frozen = read(root/(lineage+'-prediction-freeze.json'))
        assert all(sha(root/name) == digest for name, digest in frozen['files'].items())
        write(root/(lineage+'-prediction-complete.json'), dict(status='passed', commands=statuses,
              nativeReplayPassed=True, freezeSHA256=sha(root/(lineage+'-prediction-freeze.json'))))
        return frozen

    with ThreadPoolExecutor(max_workers=2) as executor:
        outputs = list(executor.map(run, LINEAGES))
    files = {}
    for result in outputs:
        assert set(files).isdisjoint(result['files'])
        files.update(result['files'])
    write(root/'prediction-output-freeze.json', dict(createdAtUnix=time.time(), folds=120,
          completed=sum(r['completed'] for r in outputs), frozenMetadataSHA256=sha(root/'freeze.json'),
          targetsOpenedForScoring=False, allNativeReplaysPassed=True, files=files))
    write(root/'prediction-complete.json', dict(status='passed', folds=120,
          completed=sum(r['completed'] for r in outputs), allNativeReplaysPassed=True,
          outputFreezeSHA256=sha(root/'prediction-output-freeze.json'), scoringPerformed=False))


if __name__ == '__main__':
    main()
