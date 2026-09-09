#!/usr/bin/env python3
"""Additional integrity checks on a completed native composition qualification."""
import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import combinations as reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('binary', 'qualified', 'pseudobulk', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    checks = json.loads((args.qualified / 'checks.json').read_text())
    assert checks['status'] == 'passed' and checks['binarySHA256'] == reference.sha(args.binary)
    assert reference.sha(args.qualified / 'prediction/report.json') == checks['nativeReportSHA256']
    args.out.mkdir(parents=True, exist_ok=False)
    results = []

    def reject(name, command):
        with (args.out / (name + '.log')).open('w') as log:
            result = subprocess.run([str(args.binary.resolve()), *map(str, command)], stdout=log, stderr=subprocess.STDOUT)
        assert result.returncode == 65, (name, result.returncode)
        results.append(name)

    training = args.qualified / 'training.json'
    before = reference.sha(training)
    reject('prepare-overwrite', ['singlecell-composition-prepare', args.pseudobulk, '--plan', args.qualified / 'selection.json', '--output', training])
    assert reference.sha(training) == before
    target = args.out / 'tampered-prediction'
    shutil.copytree(args.qualified / 'prediction', target)
    data = (target / 'report.json').read_bytes()
    # Preserve all other canonical bytes and change one predicted value.
    changed, count = re.subn(rb'("expression":\[\[)([^,\]]+)', lambda m: m[1] + str(float(m[2]) + 1).encode(), data, count=1)
    assert count == 1 and changed != data
    (target / 'report.json').write_bytes(changed)
    receipt = json.loads((target / 'receipt.json').read_text())
    receipt['result'] = dict(bytes=list(hashlib.sha256(changed).digest()))
    reference.write(target / 'receipt.json', receipt)
    reject('rehashed-prediction', ['singlecell-composition-prediction-verify', target])
    reference.write(args.out / 'checks.json', dict(status='passed', checks=results,
                                                   binarySHA256=reference.sha(args.binary),
                                                   qualifiedChecksSHA256=reference.sha(args.qualified / 'checks.json')))
    print(json.dumps(dict(status='passed', checks=results)))


if __name__ == '__main__':
    main()
