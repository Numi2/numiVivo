"""Transient errors retry the same pinned range; identity/access failures do not."""
import unittest,http.client,hashlib
from unittest.mock import patch
import run_counts as a
class Response:
 def __init__(self,status,etag=None):self.status=status;self.etag=etag or a.ETAG
 def getheader(self,name):return {'Content-Range':f'bytes 0-2/{a.SIZE}','ETag':self.etag}.get(name)
 def read(self,n):return b'abc'
class Connection:
 def __init__(self,response):self.response=response;self.requests=[];self.closed=False
 def request(self,method,path,headers):self.requests.append((method,path,headers))
 def getresponse(self):return self.response
 def close(self):self.closed=True
class Checks(unittest.TestCase):
 def tearDown(self):
  if hasattr(a.LOCAL,'connection'):del a.LOCAL.connection
 def test_transient_then_success_preserves_identity(self):
  first=Connection(Response(500));second=Connection(Response(206));records=[]
  with patch.object(a.http.client,'HTTPSConnection',side_effect=[first,second]),patch.object(a.time,'sleep') as wait:
   self.assertEqual(a.fetch(0,2,records),b'abc');wait.assert_called_once_with(1)
  self.assertTrue(first.closed);self.assertEqual(first.requests,second.requests);self.assertEqual(second.requests[0][2]['If-Match'],a.ETAG);self.assertEqual(records[0]['SHA256'],hashlib.sha256(b'abc').hexdigest());self.assertEqual(len(records[0]['transportErrors']),1)
 def test_no_retry_for_access_or_identity_failures(self):
  for status,etag in [(403,None),(412,None),(200,None),(206,'wrong')]:
   self.tearDown();connection=Connection(Response(status,etag))
   with patch.object(a.http.client,'HTTPSConnection',return_value=connection),patch.object(a.time,'sleep') as wait:
    with self.assertRaises(AssertionError):a.fetch(0,2,[])
    wait.assert_not_called();self.assertEqual(len(connection.requests),1)
 def test_transient_budget_is_finite(self):
  connections=[Connection(Response(503)) for _ in range(6)]
  with patch.object(a.http.client,'HTTPSConnection',side_effect=connections),patch.object(a.time,'sleep') as wait:
   with self.assertRaises(http.client.HTTPException):a.fetch(0,2,[])
   self.assertEqual([c.args[0] for c in wait.call_args_list],[1,2,4,8,16])
  self.assertTrue(all(c.closed for c in connections))
if __name__=='__main__':unittest.main()
