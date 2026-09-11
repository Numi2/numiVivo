"""Verify mapped real-data groups, unchanged full-cohort replay and rejections."""
from pathlib import Path
import copy,gzip,hashlib,json,os,subprocess,sys
root=Path(sys.argv[1]);repo=Path(sys.argv[2]);full=Path(sys.argv[3]);checks={}
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_bytes())
old=full/'Kang-gap-repaired/00000-00064-receipt.json';new=root/'compat-Kang/00000-00064-receipt.json';a,b=read(old),read(new)
for key in ['expandedInputSHA256','expandedOutputSHA256','compressedInputSHA256','compressedOutputSHA256']:assert a[key]==b[key],key
checks['fullCohortReplay']='first64 Kang native input and output byte-identical'
fold=root/'Kang-00';base=read(fold/'Kang-manifest.json.gz');prepared=read(fold/'Kang-prepare-state.json');cases=root/'rejected-inputs';cases.mkdir(exist_ok=False)
for name in ['excluded-donor-present','duplicate-cache-index','wrong-cache-identity']:
 path=cases/name;path.mkdir();m=copy.deepcopy(base)
 if name=='excluded-donor-present':m['excludedDonorID']=m['groups'][0]['donorID'];expected=b"manifest['excludedDonorID']"
 elif name=='duplicate-cache-index':m['groups'][1]['index']=m['groups'][0]['index'];expected=b'len(set(indices))'
 else:m['groups'][0]['index'],m['groups'][1]['index']=m['groups'][1]['index'],m['groups'][0]['index'];expected=b"group.attrs['donorID']"
 data=json.dumps(m,sort_keys=True,separators=(',',':')).encode();(path/'Kang-manifest.json.gz').write_bytes(gzip.compress(data,mtime=0));p=dict(prepared,manifestSHA256=hashlib.sha256(data).hexdigest());(path/'Kang-prepare-state.json').write_text(json.dumps(p));(path/'runtime-freeze.json').write_bytes((fold/'runtime-freeze.json').read_bytes());os.link(fold/'Kang-cells.h5',path/'Kang-cells.h5');r=subprocess.run([sys.executable,str(repo/'Tools/Omics/CountObservation/Joint/Adaptive/Full/run.py'),str(path),'Kang','1'],capture_output=True);(path/'stderr.log').write_bytes(r.stderr);assert r.returncode!=0 and expected in r.stderr and not list((path/'Kang').glob('*-input*')) and not list((path/'Kang').glob('*-receipt.json'));checks[name]={'returncode':r.returncode,'rejectedBeforeNativeInput':True}
checks['pilots']={}
for tag,origin in [('Kang-00','Kang'),('HIRISA-00','HIRISA')]:
 v=read(root/tag/origin/'verification-state.json');assert v['genesVerified']==64 and not v['errors'];checks['pilots'][tag]=v
(root/'joint-input-checks.json').write_text(json.dumps(checks,indent=2)+'\n');print('Byte replay, three pre-native rejection cases and both real-data pilots pass')
