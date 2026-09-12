from pathlib import Path
import tempfile,json,hashlib,importlib.util
spec=importlib.util.spec_from_file_location("review",Path(__file__).with_name("verify_parse_complete.py"))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
study=Path('/Users/n/numivivo-parse-logcpm-20260912')
source=Path('/Users/n/numivivo-parse-bcell-counts-20260912')
binary=Path('/Users/n/numivivo-streamed-logcpm-20260912/stream-check')
h=hashlib.sha256((study/'run.py').read_bytes()).hexdigest()
protocol=json.loads((study/'protocol.json').read_text())
checks=[]
with tempfile.TemporaryDirectory(prefix='numivivo-logcpm-admission-') as tmp:
 p=Path(tmp);(p/'run.py').write_bytes((study/'run.py').read_bytes())
 def reject(name, expected):
  try:m.review(p,source,binary,h)
  except (ValueError,OSError) as e:
   assert expected in str(e),(name,str(e));checks.append(name)
  else:raise RuntimeError('accepted '+name)
 (p/'protocol.json').write_text(json.dumps(protocol))
 reject('missing completion','complete.json')
 (p/'complete.json').write_text(json.dumps(dict(status='PASS-all-donors',predictionFitted=False,predictionScored=False,donors=[],cells=72446)))
 reject('forged complete with no donors','completed donor coverage')
 altered=dict(protocol);altered['maximumAbsoluteErrorTolerance']=1.
 (p/'protocol.json').write_text(json.dumps(altered))
 reject('relaxed tolerance','tolerance differs')
 altered=dict(protocol);altered['driverSHA256']='0'*64
 (p/'protocol.json').write_text(json.dumps(altered))
 reject('driver mismatch','driver differs')
print(json.dumps({'status':'PASS-admission-tests','checks':checks,'positiveCompleteCohortTest':False}))
