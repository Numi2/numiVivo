"""Real-data replay under group permutation and excluded-donor rejection."""
from pathlib import Path
import gzip,hashlib,json,shlex,subprocess,sys
root=Path(sys.argv[1]);runtime=json.loads((root/'runtime-freeze.json').read_text());raw=gzip.decompress((root/'inputs/Kang-00.json.gz').read_bytes());original=json.loads(raw);out=root/'contract';out.mkdir(exist_ok=False)
def run(payload,name):
 data=json.dumps(payload,sort_keys=True,separators=(',',':')).encode();(out/(name+'-input.json.gz')).write_bytes(gzip.compress(data,mtime=0));command='ulimit -c 0; exec '+shlex.quote(runtime['nativePath']);r=subprocess.run(['ssh','-4','-C','-o','ConnectTimeout=10','macmini',command],input=data,capture_output=True);(out/(name+'-stderr.log')).write_bytes(r.stderr);(out/(name+'-output.json.gz')).write_bytes(gzip.compress(r.stdout,mtime=0));return r
permuted=dict(original);permuted['groups']=list(reversed(original['groups']));r=run(permuted,'permuted');r.check_returncode();a=json.loads(gzip.decompress((root/'native/Kang-00.json.gz').read_bytes()));b=json.loads(r.stdout);assert a['inputSHA256']!=b['inputSHA256']
for side in ['control','treated']:
 a[side].pop('trainingSource');b[side].pop('trainingSource');assert a[side]==b[side]
bad=dict(original);bad['excludedDonorID']=original['groups'][0]['donorID'];r=run(bad,'excluded-donor-present');assert r.returncode!=0 and not r.stdout and b'excluded donor or paired condition contract' in r.stderr
result={'status':'passed-real-data-training-only-contract','permutedModelFields':'exactly identical apart from input/training-source hashes','excludedDonorPresentReturncode':r.returncode,'binarySHA256':runtime['binarySHA256'],'inputSHA256':hashlib.sha256(raw).hexdigest()};(root/'contract-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
