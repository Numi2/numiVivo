#!/usr/bin/env python3
"""Archive exact QL inference inputs, native outputs, and retained failures."""
import argparse,gzip,hashlib,json
from pathlib import Path
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--root',type=Path)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--verify',action='store_true')
args=parser.parse_args();sha=lambda raw:hashlib.sha256(raw).hexdigest()
if args.verify:
    manifest=json.loads((args.out/'manifest.json').read_text())
    for entry in manifest['entries']:
        p=args.out/entry['path']; assert not p.is_symlink()
        raw=p.read_bytes(); assert sha(raw)==entry['sha256'],entry['path']
        if 'logicalSHA256' in entry:assert sha(gzip.decompress(raw))==entry['logicalSHA256'],entry['path']
    assert {str(p.relative_to(args.out)) for p in args.out.rglob('*') if p.is_file()}=={e['path'] for e in manifest['entries']}|{'manifest.json'}
    print(json.dumps(dict(status='all-archive-members-verified',files=len(manifest['entries']))))
    raise SystemExit(0)
assert args.root and not args.out.exists()
summary=json.loads((args.root/'summary.json').read_text());assert summary['status']=='passed' and summary['arms']==58
args.out.mkdir(parents=True);entries=[]
for path in sorted(args.root.rglob('*')):
    if not path.is_file():continue
    assert not path.is_symlink()
    relative=path.relative_to(args.root);raw=path.read_bytes()
    compress=path.suffix in ['.log','.txt','.swift','.R'] or (path.suffix=='.json' and len(raw)>32768)
    data=gzip.compress(raw,mtime=0) if compress else raw
    target=args.out/(str(relative)+'.gz' if compress else str(relative))
    target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    entry=dict(path=str(target.relative_to(args.out)),bytes=len(data),sha256=sha(data))
    if compress or path.suffix=='.gz':entry['logicalSHA256']=sha(raw if compress else gzip.decompress(raw))
    entries.append(entry)
execution=json.loads((args.root/'execution.json').read_text())
(args.out/'manifest.json').write_text(json.dumps(dict(entries=entries,external=execution['external'],
    protocolSHA256=execution['protocolSHA256'],qualification=summary['scope']),sort_keys=True,indent=2)+'\n')
print(json.dumps(dict(files=len(entries),storedBytes=sum(e['bytes'] for e in entries))))
