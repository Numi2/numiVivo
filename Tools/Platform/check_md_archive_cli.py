#!/usr/bin/env python3
"""Verify retained long-archive fixtures through real CLI processes on macOS.

Generate the receipt with NUMIVIVO_LONG_ARCHIVE_TESTS=1 and
NUMIVIVO_TEST_ARTIFACTS set when running the named Swift archive campaign.
Each command runs in a fresh process. RSS excludes fixture generation; elapsed
times are observations without a cache-control or throughput qualification.
"""
import argparse
import hashlib
import json
import os
import pathlib
import platform
import re
import subprocess


def digest(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=pathlib.Path, required=True)
    parser.add_argument('--receipt', type=pathlib.Path, required=True)
    parser.add_argument('--out', type=pathlib.Path, required=True)
    parser.add_argument('--baseline-binary', type=pathlib.Path)
    args = parser.parse_args()
    if platform.system() != 'Darwin':
        raise RuntimeError('This measurement requires macOS /usr/bin/time -l; RSS is recorded in bytes')
    receipt_path, out = args.receipt.resolve(), args.out.resolve()
    receipt = json.loads(receipt_path.read_text())
    if (receipt.get('schema') != 'numivivo.org/test-evidence/md-trajectory-archive/v1'
            or receipt.get('status') != 'success' or receipt.get('retained') is not True):
        raise RuntimeError('Expected a successful retained Swift long-archive campaign receipt')
    root = pathlib.Path(receipt['rootPath']).resolve(strict=True)
    if out == root or root in out.parents:
        raise RuntimeError('Keep measurement outputs outside the immutable archive store')
    if out.exists() and any(out.iterdir()):
        raise RuntimeError('Use an empty output directory; existing measurements are not replaced')
    out.mkdir(parents=True, exist_ok=True)
    prefixes = {}
    for key, count in [('benchmark10000', 10000), ('original100001', 100001), ('extended100002', 100002)]:
        prefix = receipt[key]
        if (prefix['chunkCount'] != count or prefix['frameCount'] != count
                or not re.fullmatch(r'[0-9a-f]{64}', prefix['manifestFingerprint'])):
            raise RuntimeError('Unexpected retained fixture shape: ' + key)
        prefixes[key] = prefix
    binaries = {'candidate': args.binary.resolve(strict=True)}
    if args.baseline_binary:
        binaries['baseline'] = args.baseline_binary.resolve(strict=True)
    report = {
        'schema': 'numivivo.org/test-evidence/md-archive-cli/v1',
        'receiptSHA256': digest(receipt_path), 'archiveRoot': str(root),
        'platform': platform.platform(), 'rssUnit': 'bytes',
        'method': 'Fresh serial CLI processes; /usr/bin/time -l; no controlled filesystem-cache policy',
        'binaries': {key: {'path': str(path), 'sha256': digest(path)} for key, path in binaries.items()},
        'runs': [], 'passed': False,
    }

    def save():
        (out / 'checks.json').write_text(json.dumps(report, indent=2) + '\n')

    def run(label, binary, prefix_key, verify=False, maximum=None, expected=0):
        prefix = prefixes[prefix_key]
        command = ['/usr/bin/time', '-l', str(binaries[binary]), 'md-trajectory-inspect',
                   '--store', str(root), '--manifest', prefix['manifestFingerprint']]
        if verify:
            command.append('--verify')
        if maximum is not None:
            command += ['--maximum-chunks', str(maximum)]
        stdout_path, stderr_path = out / (label + '.json'), out / (label + '.stderr')
        item = {'label': label, 'binary': binary, 'command': command, 'expectedExit': expected, 'passed': False}
        report['runs'].append(item)
        save()
        with stdout_path.open('x') as stdout, stderr_path.open('x') as stderr:
            process = subprocess.run(command, stdout=stdout, stderr=stderr, timeout=600,
                                     env=dict(os.environ, LC_ALL='C'))
        item['exitCode'] = process.returncode
        stderr = stderr_path.read_text()
        rss = re.search(r'^\s*(\d+)\s+maximum resident set size\s*$', stderr, re.MULTILINE)
        timing = re.search(r'([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys', stderr)
        footprint = re.search(r'^\s*(\d+)\s+peak memory footprint\s*$', stderr, re.MULTILINE)
        if rss:
            item['maximumResidentSetBytes'] = int(rss.group(1))
        if timing:
            item.update(zip(['elapsedSeconds', 'userSeconds', 'systemSeconds'], map(float, timing.groups())))
        if footprint:
            item['peakMemoryFootprintBytes'] = int(footprint.group(1))
        save()
        if process.returncode != expected or not rss or not timing or (expected == 0 and int(rss.group(1)) <= 0):
            raise RuntimeError(label + ': unexpected exit or missing process measurements; inspect retained logs')
        if expected == 0:
            value = json.loads(stdout_path.read_text())
            if (value['indexedChunks'] != prefix['chunkCount']
                    or value['manifest']['chunkCount'] != prefix['chunkCount']
                    or value['manifest']['frameCount'] != prefix['frameCount']
                    or value['allPayloadsVerified'] is not verify):
                raise RuntimeError(label + ': CLI returned an incorrect archive summary or verification scope')
        elif stdout_path.read_text().strip() or 'limit' not in stderr:
            raise RuntimeError(label + ': expected an explicit work-limit rejection without a success report')
        item['passed'] = True
        save()

    run('candidate-index-10000', 'candidate', 'benchmark10000')
    run('candidate-index-100001', 'candidate', 'original100001')
    run('candidate-full-10000', 'candidate', 'benchmark10000', verify=True)
    run('candidate-full-100001', 'candidate', 'original100001', verify=True)
    run('candidate-extended-100002', 'candidate', 'extended100002', verify=True)
    run('candidate-explicit-limit', 'candidate', 'original100001', maximum=100000, expected=65)
    if 'baseline' in binaries:
        run('baseline-full-10000', 'baseline', 'benchmark10000', verify=True, maximum=10000)
        run('baseline-full-100001', 'baseline', 'original100001', verify=True, maximum=100001)
        run('baseline-default-limit', 'baseline', 'original100001', expected=65)
    report['passed'] = True
    save()
    print('PASS', len(report['runs']), 'real archive CLI invocations; process measurements retained in', out)


if __name__ == '__main__':
    main()
