from pathlib import Path
import hashlib,json,tarfile
r=Path(__file__).parent;m=json.loads((r/'manifest.json').read_text());p=r/'evidence.tar.gz'
assert p.stat().st_size==m['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==m['sha256']
with tarfile.open(p) as t:
 def read(n):return t.extractfile(n).read()
 def obj(n):return json.loads(read(n))
 terminal=obj('terminal.json');q=obj('qualification.json');protocol=obj('protocol.json')
 assert [x['step'] for x in terminal['steps']]==['pca','mnn','verify']
 assert all(x['exitCode']==0 for x in terminal['steps'])
 assert q['status']=='PASS'
 assert [x['case'] for x in q['rejections']]==['existing-output','work-budget','tampered-anchor']
 assert all(x['exitCode']==65 for x in q['rejections'])
 assert len(q['byteEquivalence'])==4 and all(x['fresh']==x['prior'] for x in q['byteEquivalence'].values())
 assert obj('budget-plan.json')['mnn']['maximumWork']==1
 assert 'SyntaxError' in read('driver-attempt-1.log').decode()
 compile(read('driver.py'),'driver.py','exec');compile(read('rejections.py'),'rejections.py','exec')
 assert protocol['cliSHA256']=='04b2f225157bb072703b5b69c9f2a94b366358ce9b9080de7310f036726e064a'
print('PASS: retained CLI lifecycle evidence; not fresh execution')
