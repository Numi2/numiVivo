from pathlib import Path
import hashlib,json,os,time,urllib.request,zlib
r=Path(__file__).parent
url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE181nnn/GSE181897/suppl/GSE181897_concat.4.raw.h5ad.gz'
p=r/'GSE181897_concat.4.raw.h5ad.gz'; part=p.with_suffix('.gz.partial')
assert not p.exists() and not part.exists()
expected=1011162509; expected_decoded=3063713137
assert os.statvfs(r).f_bavail*os.statvfs(r).f_frsize>expected+700_000_000
start=time.time(); compressed=hashlib.sha256();decoded=hashlib.sha256();n=m=0; d=zlib.decompressobj(16+zlib.MAX_WBITS); last=start
with urllib.request.urlopen(url,timeout=45) as f,part.open('xb') as out:
 assert f.status==200 and int(f.headers['Content-Length'])==expected
 headers=dict(f.headers)
 while b:=f.read(1024*1024):
  compressed.update(b);n+=len(b);x=d.decompress(b);m+=len(x)
  assert m<=expected_decoded and n<=expected
  if os.statvfs(r).f_bavail*os.statvfs(r).f_frsize<600_000_000:raise RuntimeError('Capacity floor reached; retain partial')
  decoded.update(x);out.write(b)
  if time.time()-last>15:print(json.dumps(dict(compressedBytes=n,decodedBytes=m,elapsed=time.time()-start)),flush=True);last=time.time()
 x=d.flush();decoded.update(x);m+=len(x)
 assert d.eof and not d.unused_data and n==expected and m==expected_decoded
 out.flush();os.fsync(out.fileno())
assert part.stat().st_size==expected
assert part.open('rb').read(2)==b'\x1f\x8b'
part.rename(p)
receipt=dict(url=url,headers=headers,compressedBytes=n,decodedBytes=m,compressedSHA256=compressed.hexdigest(),decodedSHA256=decoded.hexdigest(),elapsedSeconds=time.time()-start,completedUnix=time.time(),storage='Original gzip retained; decoded in transit with CRC and SHA256 checks; no decoded duplicate')
(r/'download.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
