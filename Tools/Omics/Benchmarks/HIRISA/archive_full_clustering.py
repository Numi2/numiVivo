#!/usr/bin/env python3
"""Preserve full native clustering/replay and original-reference comparison."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 r=a.study/'optimized-clustering';c=read(r/'collection-complete.json');n=read(r/'native-complete.json');collection=read(r/'collected/collection.json');check=read(r/'collected/independent-check/checks.json')
 assert c['status']=='passed' and c['nativeReplayPassed'] and c['independentAllCellPartitionCheckPassed']
 assert c['nativeCompletionSHA256']==sha(r/'native-complete.json') and c['independentCheckSHA256']==sha(r/'collected/independent-check/checks.json')
 assert len(n['commands'])==4 and all(x['returnCode']==0 for x in n['commands']) and n['binarySHA256']==c['nativeRuntimeSHA256']
 assert check['status']=='passed' and check['cells']==1612594 and check['disconnectedCommunities']==0
 for name,identity in collection['remotePayloads'].items():
  path=r/'collected/bundle'/name;assert dict(bytes=path.stat().st_size,SHA256=sha(path))==identity,name
 paths=[('run/'+x.name,x) for x in r.iterdir() if x.is_file() and x.suffix in {'.json','.log','.py','.jsonl'}]
 paths += [('collected/'+str(x.relative_to(r/'collected')),x) for x in (r/'collected').rglob('*') if x.is_file() and 'bundle' not in x.relative_to(r/'collected').parts]
 paths += [('bundle/'+name,r/'collected/bundle'/name) for name in collection['copied']]
 paths += [('tools/'+name,Path(__file__).with_name(name)) for name in ['archive_full_clustering.py','reference_clustering.py','collect_clustering.py']]
 # Earlier graph/storage archives supply byte-identical parent payloads; current
 # receipt identities are retained in the newly collected bundle, not relabeled.
 upstream={}
 for folder in ['2026-09-10-graph','2026-09-10-storage-access','2026-09-10-clustering-reference']:
  previous=read(Path(__file__).parent/'evidence'/folder/'manifest.json')
  for x in previous['records']:upstream[x['sourceSHA256']]=dict(archive=folder,sourcePath=x['sourcePath'])
  for name,x in previous.get('fullSources',{}).items():upstream[x['SHA256']]=dict(archive=folder,sourcePath=name)
 restoration={}
 for name in collection['reused']:
  identity=collection['remotePayloads'][name];assert identity['SHA256'] in upstream,name
  restoration[name]=dict(**identity,**upstream[identity['SHA256']])
 assert len({n for n,_ in paths})==len(paths)
 for _,source in paths:assert source.is_file() and not source.is_symlink()
 a.out.mkdir(parents=True,exist_ok=False);records=[]
 for name,source in sorted(paths):
  raw=source.read_bytes();encoded=gzip.compress(raw,compresslevel=6,mtime=0);target=a.out/(name+'.gz');target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(encoded)
  digest=hashlib.sha256(raw).hexdigest();assert sha(source)==digest
  records.append(dict(sourcePath=name,sourceBytes=len(raw),sourceSHA256=digest,storedPath=name+'.gz',storedBytes=len(encoded),storedSHA256=hashlib.sha256(encoded).hexdigest(),gzipEncoded=True))
 manifest=dict(schemaVersion=1,records=records,reusedParentPayloadRestoration=restoration,scope='Complete 1612594-cell current-runtime native clustering, replay and every-label/edge/connectivity reference checks. Thirty communities are graph partitions, not calibrated cell types. Original mapped-reader baseline remains separate and is not asserted complete. Exact original igraph partitions and reused graph/PCA payloads restore from earlier archives; new runtime receipts retained unchanged.')
 (a.out/'manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n');print(json.dumps(dict(status='archived',members=len(records),storedBytes=sum(x['storedBytes'] for x in records),manifestSHA256=sha(a.out/'manifest.json'))))
if __name__=='__main__':main()
