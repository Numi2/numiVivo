"""Real query controls and malformed identities must fail before prediction."""
from pathlib import Path
import copy,gzip,hashlib,json,shlex,subprocess,sys
root=Path(sys.argv[1]);out=root/'contract';out.mkdir(exist_ok=False);freeze=json.loads((root/'runtime-freeze.json').read_text());lines=gzip.decompress((root/'Kang-00/00000-00064-input.jsonl.gz').read_bytes()).splitlines();header=json.loads(lines[0]);gene=next(json.loads(l) for l in lines[1:] if json.loads(l)['status']=='boundedContinuousLikelihood');header.update(firstFeatureIndex=gene['featureIndex'],featureIDs=[gene['featureID']]);cases={}
h=copy.deepcopy(header);h['queryDonorID']=h['trainingDonorIDs'][0];cases['training-donor-query']=(h,gene,'prediction donor or feature header')
g=copy.deepcopy(gene);g['featureID']='incorrect-feature';cases['wrong-feature']=(header,g,'prediction feature sequence')
g=copy.deepcopy(gene);g['model']['model']['trainingDonorIDs'][0]='incorrect-training-donor';cases['wrong-model-donors']=(header,g,'prediction model training identity')
g=copy.deepcopy(gene);g.update(bins=[0,0],counts=[1,1]);cases['duplicate-count-bin']=(header,g,'prediction count bins')
report={};native=freeze['nativePath'];check='import hashlib,pathlib;assert hashlib.sha256(pathlib.Path('+repr(native)+').read_bytes()).hexdigest()=='+repr(freeze['binarySHA256'])
for name,(h,g,error) in cases.items():
 data=b'\n'.join(json.dumps(x,sort_keys=True,separators=(',',':')).encode() for x in [h,g])+b'\n';(out/(name+'-input.jsonl.gz')).write_bytes(gzip.compress(data,mtime=0));p=subprocess.run(['ssh','-4','-o','ConnectTimeout=10','macmini','python3 -c '+shlex.quote(check)+' && '+shlex.quote(native)],input=data,capture_output=True);(out/(name+'-stderr.log')).write_bytes(p.stderr);assert p.returncode!=0 and not p.stdout and error.encode() in p.stderr;report[name]={'returncode':p.returncode,'rejectedBeforePrediction':True,'expectedError':error}
(out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(report)
