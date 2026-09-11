"""Restore one complete native bundle from the retained shared object store."""
import argparse, gzip, hashlib, json
from pathlib import Path, PurePosixPath

parser=argparse.ArgumentParser()
parser.add_argument('--native',type=Path,required=True)
parser.add_argument('--bundle',required=True)
parser.add_argument('--out',type=Path,required=True)
args=parser.parse_args()
assert not args.out.exists()
freeze=json.loads((args.native/'prediction-freeze.json').read_bytes())
raw=(args.native/'retention.json').read_bytes()
assert hashlib.sha256(raw).hexdigest()==freeze['files']['retention.json']
index=json.loads(raw); members=index['bundles'][args.bundle]
for name,item in members.items():
    p=PurePosixPath(name)
    assert not p.is_absolute() and '..' not in p.parts and str(p)==name
    assert item['SHA256'] in index['objects'] and item['bytes']==index['objects'][item['SHA256']]['bytes']
args.out.mkdir(parents=True)
for name,item in members.items():
    source=args.native/'objects'/(item['SHA256']+'.gz')
    expected=index['objects'][item['SHA256']]
    compressed=hashlib.sha256(source.read_bytes()).hexdigest()
    assert compressed==expected['compressedSHA256']==freeze['files']['objects/'+source.name]
    dest=args.out/name;dest.parent.mkdir(parents=True,exist_ok=True)
    h=hashlib.sha256();size=0
    with gzip.open(source,'rb') as reader,dest.open('xb') as writer:
        while block:=reader.read(1048576): writer.write(block);h.update(block);size+=len(block)
    assert h.hexdigest()==item['SHA256'] and size==item['bytes']
assert {str(p.relative_to(args.out)) for p in args.out.rglob('*') if p.is_file()}==set(members)
print(json.dumps(dict(status='restored-complete-native-bundle',bundle=args.bundle,files=len(members),
    bytes=sum(x['bytes'] for x in members.values()),native=str(args.native),output=str(args.out))))
