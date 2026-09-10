#!/usr/bin/env python3
"""Run every frozen experimental request through the existing native cohort owner."""
import argparse
import gzip
import json
import os
from pathlib import Path
import platform
import subprocess
import time
from acquire import digest


def contains(actual, expected):
    return all(k in actual and contains(actual[k], v) for k, v in expected.items()) if isinstance(expected, dict) else actual == expected

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs', type=Path, required=True)
    p.add_argument('--binary', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args();root = a.inputs;binary = a.binary.resolve();out = a.output
    source = json.loads((root / 'receipt.json').read_text())
    assert source['cases'] == 48 and source['cohorts'] == 16
    assert source['protocolSHA256'] == digest(Path(__file__).with_name('PROTOCOL.md'))
    assert source['executionProtocolSHA256'] == digest(Path(__file__).with_name('INFERENCE_EXECUTION.md'))
    assert digest(root / 'requests/manifest.json') == source['requestManifestSHA256']
    binary_hash = digest(binary)
    scope = []
    for line in binary.with_name('sources.sha256').read_text().splitlines():
        sha, path = line.split(maxsplit=1)
        assert digest(Path(path)) == sha, ('native owner source differs from build scope', path)
        scope.append({'path': path, 'sha256': sha})
    out.mkdir(exist_ok=True)
    lock = out / 'driver.lock'
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, str(os.getpid()).encode());os.close(fd)
    environment = {'binarySHA256': binary_hash, 'sources': scope, 'machine': platform.machine(),
                   'macOS': platform.mac_ver(), 'driverSHA256': digest(Path(__file__))}
    (out / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    try:
        manifest = json.loads((root / 'requests/manifest.json').read_text())
        for case in manifest:
            request = root / 'requests' / case['path'];report = root / case['sourceReport']
            assert digest(request) == case['requestSHA256'] and digest(report) == case['sourceReportSHA256']
            identity = {'binarySHA256': binary_hash, 'sourceReportSHA256': case['sourceReportSHA256'],
                        'protocolSHA256': source['protocolSHA256'], 'executionProtocolSHA256': source['executionProtocolSHA256'],
                        'requestSHA256': case['requestSHA256']}
            name = case['cohort'] + '-' + case['method'];directory = out / name
            if (directory / 'receipt.json').exists():
                saved = json.loads((directory / 'receipt.json').read_text())
                assert all(saved[k] == v for k, v in identity.items())
                assert saved['outputSHA256'] == digest(directory / 'output.json.gz')
                print(name, 'retained terminal checkpoint', flush=True);continue
            directory.mkdir(exist_ok=False)
            command = ['/usr/bin/time', '-l', str(binary), str(report), str(request)]
            start = time.time()
            with (directory / 'output.json').open('wb') as stdout, (directory / 'stderr.log').open('wb') as stderr:
                completed = subprocess.run(command, stdout=stdout, stderr=stderr)
            payload = (directory / 'output.json').read_bytes()
            (directory / 'output.json.gz').write_bytes(gzip.compress(payload, mtime=0))
            error = None
            try:
                result = json.loads(payload)
                assert result['sourceReportSHA256'] == identity['sourceReportSHA256']
                assert contains(result['request'], json.loads(request.read_text()))
                error = result.get('error')
            except Exception as failure:
                error = 'Invalid or missing harness output: ' + str(failure)
            record = {**case, **identity, 'command': command, 'returnCode': completed.returncode,
                      'seconds': time.time() - start, 'error': error,
                      'outputSHA256': digest(directory / 'output.json.gz')}
            (directory / 'receipt.json').write_text(json.dumps(record, sort_keys=True, indent=2) + '\n')
            (directory / 'output.json').unlink()
            print(name, 'returnCode', completed.returncode, 'seconds', round(record['seconds'], 2), 'error', error, flush=True)
        (out / 'driver-complete.json').write_text(json.dumps({'attemptedCases': len(manifest), 'allCasesTerminal': True}) + '\n')
    finally:
        lock.unlink()

if __name__ == '__main__':
    main()
