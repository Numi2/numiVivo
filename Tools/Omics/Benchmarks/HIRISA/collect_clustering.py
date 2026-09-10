#!/usr/bin/env python3
"""Collect a native cluster result while reusing exactly equal local graph bytes.

Only the result and new receipts are transferred. Every reused remote payload
is hashed before and after collection. Native replay remains a separate gate.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from reference_clustering import GRAPH_FILES, PCA_FILES

CLUSTER_FILES = ['execution.json', 'plan.json', 'receipt.json', 'result.json']
REUSED = ['input/'+name for name in GRAPH_FILES]
REUSED += ['input/input/'+name for name in PCA_FILES]
COPIED = CLUSTER_FILES+['input/receipt.json', 'input/input/receipt.json']


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--host', required=True)
    parser.add_argument('--remote-bundle', required=True)
    for name in ['reference-graph', 'protocol', 'graph-check', 'reference', 'out']:
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    names = COPIED+REUSED

    def remote_snapshot():
        # No third-party package is needed on the source host. Before/after stat
        # signatures reject a file that changes while its digest is collected.
        code = '''import hashlib,json,os,stat
from pathlib import Path
root=Path(ROOT);names=NAMES;result={}
for name in names:
 p=root/name
 assert p.resolve().is_relative_to(root.resolve()) and not p.is_symlink(),name
 before=p.stat();assert stat.S_ISREG(before.st_mode),name
 h=hashlib.sha256()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 after=p.stat()
 def signature(s):return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)
 assert signature(before)==signature(after),('Changed while hashing',name)
 result[name]={'bytes':after.st_size,'SHA256':h.hexdigest()}
print(json.dumps(result,sort_keys=True))
'''.replace('ROOT', repr(args.remote_bundle)).replace('NAMES', repr(names))
        return json.loads(subprocess.check_output(
            ['ssh', '-4', args.host, 'python3 -c '+shlex.quote(code)], text=True))

    before = remote_snapshot()
    assert set(before) == set(names)
    write(args.out/'remote-before.json', before)
    # Every parent payload is from the original qualified local graph. New
    # runtime-bound receipts are copied, never synthesized from the old ones.
    sources = {}
    for name in REUSED:
        source = args.reference_graph/name.removeprefix('input/')
        assert source.is_file() and not source.is_symlink(), name
        assert dict(bytes=source.stat().st_size, SHA256=sha(source)) == before[name], name
        sources[name] = source
    bundle = args.out/'bundle'
    bundle.mkdir()
    for name, source in sources.items():
        destination = bundle/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(source, destination)
    filelist = args.out/'transfer-files.txt'
    filelist.write_text('\n'.join(COPIED)+'\n')
    command = ['rsync', '-r', '-e', 'ssh -4', '--files-from='+str(filelist),
               args.host+':'+args.remote_bundle.rstrip('/')+'/', str(bundle)+'/']
    subprocess.run(command, check=True)
    assert {str(p.relative_to(bundle)) for p in bundle.rglob('*') if p.is_file()} == set(names)
    for name in names:
        source = bundle/name
        assert dict(bytes=source.stat().st_size, SHA256=sha(source)) == before[name], name
    after = remote_snapshot()
    write(args.out/'remote-after.json', after)
    assert before == after, 'Remote bundle changed during collection'
    collection = dict(status='collected', remoteHost=args.host, remoteBundle=args.remote_bundle,
                      remotePayloads=before, copied= COPIED, reused=REUSED,
                      copiedLogicalBytes=sum(before[name]['bytes'] for name in COPIED),
                      reusedLogicalBytes=sum(before[name]['bytes'] for name in REUSED),
                      seconds=time.time()-started,
                      scope='All collected and reused bytes match two remote snapshots. Exact original receipts are preserved. Native replay and independent numerical acceptance are separate gates.')
    write(args.out/'collection.json', collection)
    checker = Path(__file__).with_name('reference_clustering.py')
    command = [sys.executable, str(checker), 'check', '--bundle', str(bundle),
               '--reference-graph', str(args.reference_graph), '--protocol', str(args.protocol),
               '--graph-check', str(args.graph_check), '--reference', str(args.reference),
               '--out', str(args.out/'independent-check')]
    with (args.out/'independent-check.log').open('x') as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    write(args.out/'independent-check-status.json', dict(returnCode=result.returncode, command=command,
                                                        checkerSHA256=sha(checker)))
    assert result.returncode == 0, 'Independent gate failed; collection and failure log retained'
    checks = json.loads((args.out/'independent-check/checks.json').read_text())
    assert checks['status'] == 'passed'
    complete = dict(status='passed', collectionSHA256=sha(args.out/'collection.json'),
                    independentCheckSHA256=sha(args.out/'independent-check/checks.json'),
                    cells=checks['cells'], nativeReplay='separate-unasserted-gate',
                    collectorSHA256=sha(Path(__file__)))
    write(args.out/'complete.json', complete)
    print(json.dumps({**complete, 'copiedLogicalBytes': collection['copiedLogicalBytes'],
                      'reusedLogicalBytes': collection['reusedLogicalBytes']}))


if __name__ == '__main__':
    main()
