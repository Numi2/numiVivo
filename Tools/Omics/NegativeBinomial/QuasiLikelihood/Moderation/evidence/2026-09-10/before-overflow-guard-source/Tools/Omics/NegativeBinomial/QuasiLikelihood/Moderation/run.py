#!/usr/bin/env python3
"""Qualify all frozen native QL inputs with checkpointed, remote-spooled output."""
import argparse, gzip, hashlib, json, os, shlex, subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs', type=Path, required=True)
parser.add_argument('--out', type=Path, required=True)
parser.add_argument('--remote-out', required=True)
parser.add_argument('--binary', required=True)
parser.add_argument('--limit', type=int)
args = parser.parse_args()
source = Path(__file__).parent
sha = lambda b: hashlib.sha256(b).hexdigest()
ssh = lambda command, **kw: subprocess.run(['ssh', '-4', 'macmini', command], **kw)
binary_sha = ssh('shasum -a 256 '+shlex.quote(args.binary), check=True, capture_output=True).stdout.decode().split()[0]
protocol_sha = sha((source/'PROTOCOL.md').read_bytes())
checker_sha = sha((source/'reference.R').read_bytes())
runs = json.loads((args.inputs/'complete.json').read_text())
assert len(runs) == 58
args.out.mkdir(parents=True, exist_ok=True)
receipts = []
for run in runs[:args.limit]:
    relative = Path(run['case'])/run['method']
    dest = args.out/relative
    raw = (args.inputs/relative/'input.json.gz').read_bytes()
    assert sha(raw) == run['inputSHA256']
    identity = dict(inputSHA256=sha(raw), binarySHA256=binary_sha,
                    protocolSHA256=protocol_sha, checkerSHA256=checker_sha)
    receipt_file = dest/'receipt.json'
    if receipt_file.exists():
        receipt = json.loads(receipt_file.read_text())
        assert all(receipt[k] == v for k, v in identity.items())
        assert sha((dest/'native.json.gz').read_bytes()) == receipt['nativeSHA256']
        assert sha((dest/'reference.json.gz').read_bytes()) == receipt['referenceSHA256']
        receipts.append(receipt)
        print('checkpoint',run['case'],run['method'],receipt['status'],flush=True)
        continue
    assert not dest.exists(), 'retain incomplete attempt; inspect remote state before retrying'
    dest.mkdir(parents=True)
    (dest/'input.json.gz').write_bytes(raw)
    remote = args.remote_out+'/'+relative.as_posix()
    q = shlex.quote
    upload = ssh('mkdir -p '+q(str(Path(remote).parent))+' && mkdir '+q(remote)+' && cat > '+q(remote+'/input.json.gz'),
                 input=raw, capture_output=True)
    (dest/'upload.log').write_bytes(upload.stderr)
    assert upload.returncode == 0
    command = ('set -o pipefail; gzip -dc '+q(remote+'/input.json.gz')+
        ' | /usr/bin/time -l '+q(args.binary)+' > '+q(remote+'/native.json')+' 2> '+q(remote+'/native.log')+
        ' && gzip -n '+q(remote+'/native.json'))
    native_run = ssh(command,capture_output=True)
    (dest/'transport.log').write_bytes(native_run.stdout+native_run.stderr)
    assert native_run.returncode == 0, 'native process failed; retained remote output'
    for name in ['native.json.gz','native.log']:
        raw_output = ssh('cat '+q(remote+'/'+name),check=True,capture_output=True).stdout
        (dest/name).write_bytes(raw_output)
    native = json.loads(gzip.decompress((dest/'native.json.gz').read_bytes()))
    if 'fit' not in native:
        print('native failure',relative,native,flush=True)
        raise SystemExit(1)
    with (dest/'reference.log').open('wb') as log:
        reference_run = subprocess.run(['/opt/homebrew/bin/Rscript',str(source/'reference.R'),
            str(dest/'input.json.gz'),str(dest/'native.json.gz'),str(dest/'reference.json')],
            stdout=log,stderr=log,env={**os.environ,'R_LIBS_USER':'/Users/home/numivivo-r-library-20260909'})
    if not (dest/'reference.json').exists():
        raise RuntimeError('reference did not produce output; inspect retained log')
    reference = json.loads((dest/'reference.json').read_text())
    reference_bytes = gzip.compress((dest/'reference.json').read_bytes(),mtime=0)
    (dest/'reference.json.gz').write_bytes(reference_bytes)
    (dest/'reference.json').unlink()
    receipt = dict(**identity,case=run['case'],method=run['method'],features=run['genes'],
        nativeSHA256=sha((dest/'native.json.gz').read_bytes()),referenceSHA256=sha(reference_bytes),
        seconds=native['seconds'],profiles=len(native['fit']['profiles']),
        status=reference['status'],referenceReturnCode=reference_run.returncode)
    receipt_file.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    receipts.append(receipt)
    print(run['case'],run['method'],receipt['status'],round(native['seconds'],3),'seconds',flush=True)
(args.out/'complete.json').write_text(json.dumps(receipts,sort_keys=True,indent=2)+'\n')
raise SystemExit(0 if all(r['status']=='passed' for r in receipts) else 1)
