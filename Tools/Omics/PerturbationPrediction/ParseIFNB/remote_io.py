import os,io,json,urllib.request,urllib.parse,http.client,ssl,collections,hashlib,time
from pathlib import Path
import h5py
ROOT=Path(os.environ['NUMIVIVO_PARSE_STUDY']);source=json.loads((ROOT/'sources/asset-head.json').read_text());URL=source['url'];ETAG=source['headers']['ETag'];SIZE=int(source['headers']['Content-Length'])
class Remote(io.RawIOBase):
 def __init__(self):self.position=0;self.cache=collections.OrderedDict();self.requests=[];self.connection=None
 def fetch(self,start,stop,timeout):
  if self.connection is None:
   parsed=urllib.parse.urlsplit(URL);self.connection=http.client.HTTPSConnection(parsed.netloc,timeout=timeout,context=ssl.create_default_context())
  self.connection.request('GET',urllib.parse.urlsplit(URL).path,headers={'Range':f'bytes={start}-{stop}','If-Match':ETAG,'Accept-Encoding':'identity'})
  r=self.connection.getresponse()
  assert r.status==206 and r.getheader('Content-Range')==f'bytes {start}-{stop}/{SIZE}' and r.getheader('ETag')==ETAG
  data=r.read(stop-start+2);assert len(data)==stop-start+1
  return data
 def readable(self):return True
 def seekable(self):return True
 def seek(self,offset,whence=0):
  self.position=offset if whence==0 else self.position+offset if whence==1 else SIZE+offset
  assert 0<=self.position<=SIZE;return self.position
 def tell(self):return self.position
 def read(self,n=-1):
  assert 0<=n<=32000000,n
  end=min(SIZE,self.position+n);result=[];block=262144
  if n>=1048576:
   start=self.position;stop=end-1;req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{stop}','If-Match':ETAG,'Accept-Encoding':'identity'});before=time.time()
   data=self.fetch(start,stop,90)
   self.position=end;self.requests.append(dict(start=start,stop=stop,SHA256=hashlib.sha256(data).hexdigest(),seconds=time.time()-before));print('direct range',len(data),flush=True);return data
  while self.position<end:
   start=self.position//block*block
   if start not in self.cache:
    stop=min(SIZE,start+block)-1;req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{stop}','If-Match':ETAG,'Accept-Encoding':'identity'});before=time.time()
    data=self.fetch(start,stop,45)
    self.cache[start]=data;self.requests.append(dict(start=start,stop=stop,SHA256=hashlib.sha256(data).hexdigest(),seconds=time.time()-before));print('range',len(self.requests),start,flush=True)
    if len(self.cache)>512:self.cache.popitem(last=False)
   data=self.cache[start];self.cache.move_to_end(start);take=min(end-self.position,len(data)-(self.position-start));result.append(data[self.position-start:self.position-start+take]);self.position+=take
  return b''.join(result)
 def readinto(self,b):
  data=self.read(len(b));b[:len(data)]=data;return len(data)
