"""Losslessly retain the old native report through its exact membership-only difference."""
import argparse,gzip,hashlib,json
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()

def observation_span(reader):
 # Both reports are already fully parsed and compared by the native comparison.
 # Bound this additional archive transform to the initial design/observation span.
 prefix=bytearray();marker=b'"observations":';start=None;cursor=0;depth=0;quoted=False;escaped=False
 while True:
  block=reader.read(1048576);assert block,'Missing complete observation value';prefix.extend(block)
  assert len(prefix)<=16777216,'Observation span exceeds this dated archive contract'
  if start is None:
   design=prefix.find(b'"design":');position=prefix.find(marker,design) if design>=0 else -1
   if position<0:continue
   start=position+len(marker)
   while start<len(prefix) and prefix[start] in b' \t\r\n':start+=1
   if start==len(prefix):start=None;continue
   assert prefix[start]==ord('[');cursor=start
  while cursor<len(prefix):
   byte=prefix[cursor]
   if quoted:
    if escaped:escaped=False
    elif byte==92:escaped=True
    elif byte==34:quoted=False
   elif byte==34:quoted=True
   elif byte in (91,123):depth+=1
   elif byte in (93,125):
    depth-=1
    if depth==0:
     end=cursor+1;value=bytes(prefix[start:end]);assert len(json.loads(value))==24;return start,end,value
   cursor+=1

def chunks(source,delta):
 descriptor=json.loads((delta/'descriptor.json').read_text());assert sha(source)==descriptor['sourceReportSHA256']
 assert sha(delta/'observations.json.gz')==descriptor['storedObservationSHA256']
 payload=gzip.decompress((delta/'observations.json.gz').read_bytes());assert hashlib.sha256(payload).hexdigest()==descriptor['observationSHA256']
 with source.open('rb') as f:
  start,end,original=observation_span(f)
 assert [start,end]==descriptor['sourceObservationSpan'] and hashlib.sha256(original).hexdigest()==descriptor['sourceObservationSHA256']
 with source.open('rb') as f:
  remaining=start
  while remaining:
   data=f.read(min(remaining,1048576));assert data;remaining-=len(data);yield data
  for offset in range(0,len(payload),1048576):yield payload[offset:offset+1048576]
  f.seek(end)
  for data in iter(lambda:f.read(1048576),b''):yield data

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','verify','restore']);p.add_argument('--study',type=Path,required=True);p.add_argument('--out',type=Path);a=p.parse_args();s=a.study.resolve();source=s/'native-file/report.json';delta=s/'baseline-delta'
 if a.action=='prepare':
  comparison=json.loads((s/'native-comparison.json').read_text());assert comparison['status']=='passed-full-parse-native-aggregate-handoff'
  assert comparison['allNonDesignJSONValuesByteExact'] and comparison['allNonMembershipDesignJSONValuesByteExact']
  assert not delta.exists();delta.mkdir()
  with source.open('rb') as f:start,end,current=observation_span(f)
  with gzip.open(s/'resident-report.json.gz','rb') as f:_,_,payload=observation_span(f)
  (delta/'observations.json.gz').write_bytes(gzip.compress(payload,mtime=0))
  expected=json.loads((s/'resident-execution.json').read_text())
  descriptor=dict(schemaVersion=1,encoding='replace-design-observations-v1',sourceReportSHA256=sha(source),sourceObservationSpan=[start,end],sourceObservationSHA256=hashlib.sha256(current).hexdigest(),observationSHA256=hashlib.sha256(payload).hexdigest(),storedObservationSHA256=sha(delta/'observations.json.gz'),nativeReportSHA256=expected['uncompressedReportSHA256'],nativeReportBytes=expected['uncompressedReportBytes'],originalTransportGzipSHA256=sha(s/'resident-report.json.gz'),scope='Reconstructs the exact original native JSON bytes. Original Python gzip transport encoding is not required or reproduced.')
  (delta/'descriptor.json').write_text(json.dumps(descriptor,indent=2)+'\n')
 descriptor=json.loads((delta/'descriptor.json').read_text());h=hashlib.sha256();size=0
 if a.action=='prepare':
  with gzip.open(s/'resident-report.json.gz','rb') as expected:
   for block in chunks(source,delta):assert expected.read(len(block))==block;h.update(block);size+=len(block)
   assert expected.read(1)==b''
 elif a.action=='verify':
  for block in chunks(source,delta):h.update(block);size+=len(block)
 else:
  assert a.out is not None and not a.out.exists()
  with a.out.open('xb') as raw,gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as output:
   for block in chunks(source,delta):h.update(block);size+=len(block);output.write(block)
 assert h.hexdigest()==descriptor['nativeReportSHA256'] and size==descriptor['nativeReportBytes']
 result=dict(status='passed-exact-native-report-reconstruction',nativeReportBytes=size,nativeReportSHA256=h.hexdigest(),originalBytesComparedDirectly=a.action=='prepare',storedDeltaBytes=sum(p.stat().st_size for p in delta.iterdir() if p.is_file()),originalTransportEncodingRetained=False)
 if a.action=='prepare':(delta/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result))

if __name__=='__main__':main()
