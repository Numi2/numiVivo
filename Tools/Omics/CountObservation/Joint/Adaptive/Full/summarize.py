"""Audit complete source-gene coverage and retain every solver status."""
from pathlib import Path
import collections,gzip,hashlib,json,re,sys,time
root=Path(sys.argv[1]);origins=sys.argv[2:];assert origins

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1<<20):h.update(b)
 return h.hexdigest()

def audit(folder,expected):
 statuses=collections.Counter();failures=[];identity=[];peak=0;seconds=0;refs=0;values=0;leaves=0;error=0
 receipts=sorted(folder.glob('*-receipt.json'));entries=[]
 for p in receipts:
  r=json.loads(p.read_text());tag=p.name.removesuffix('-receipt.json');output=folder/(tag+'-output.jsonl.gz');inp=folder/(tag+'-input.jsonl.gz');vr=folder/(tag+'-verification.json.gz');v=json.loads(gzip.decompress(vr.read_bytes()))
  assert r['status']=='completed' and sha(inp)==r['compressedInputSHA256'] and sha(output)==r['compressedOutputSHA256'];assert v['status']=='passed-independent-numerics' and v['nativeCompressedSHA256']==r['compressedOutputSHA256']
  indices=[];counts=collections.Counter()
  with gzip.open(output,'rt') as f:
   for line in f:
    row=json.loads(line);indices.append(row['featureIndex']);identity.append((row['featureIndex'],row['featureID'],row['status']));counts[row['status']]+=1
    if row['status'] not in ['boundedContinuousLikelihood','unavailableCellDispersion']:
     a=row.get('model',{});c=a.get('certificate',{});failures.append({'featureIndex':row['featureIndex'],'featureID':row['featureID'],'status':row['status'],'certificateStatus':c.get('status'),'finiteGap':a.get('model',{}).get('meanLogLikelihoodGap'),'detail':row.get('detail')})
  assert indices==list(range(r['first'],r['last'])) and dict(counts)==r['fitStates'];statuses.update(counts);seconds+=r['seconds'];refs+=1
  text=(folder/(tag+'-stderr.log')).read_text();match=re.search(r'(\d+)\s+maximum resident set size',text);assert match;peak=max(peak,int(match[1]));entries.append({'receipt':p.name,'SHA256':sha(p),'inputSHA256':r['compressedInputSHA256'],'outputSHA256':r['compressedOutputSHA256'],'verificationSHA256':sha(vr)})
 assert [i[0] for i in identity]==list(range(expected))
 state=json.loads((folder/'verification-state.json').read_text());assert state['genesVerified']==expected and not state['errors'] and state['fitStates']==dict(statuses)
 return dict(genes=expected,fitStates=dict(statuses),nonconverged=failures,shards=refs,peakNativeRSSBytes=peak,sumShardWallSeconds=seconds,independentVerification=state,shardArtifacts=entries),identity

results={'status':'completed-listed-origins','createdUnix':time.time(),'origins':{},'scope':'All source genes and previously admitted training cells, fixed cell dispersions; no query prediction or new outcome validation','biologicalStatus':'Unchanged: bounded mean RNA response evidence; general biological outcomes and treated-interval calibration unqualified.'}
for origin in origins:
 prep=json.loads((root/(origin+'-prepare-state.json')).read_text());assert prep['status']=='completed'
 completed=root/(origin+'-gap-repaired');watch=json.loads((completed/'verification-watch.json').read_text());assert watch['status']=='completed-all-workers-and-verification'
 repaired,identity=audit(completed,prep['genes']);entry={'preparation':prep,'repaired':repaired,'runtime':json.loads((root/'gap-repair/runtime-freeze.json').read_text())}
 baseline=root/origin
 if (baseline/'state.json').exists() and json.loads((baseline/'state.json').read_text())['genesCompleted']==prep['genes']:
  b,old=audit(baseline,prep['genes']);entry['baseline']=b;assert [(i,n) for i,n,_ in old]==[(i,n) for i,n,_ in identity]
  transition=collections.Counter((a[2],b[2]) for a,b in zip(old,identity));entry['statusTransitions']=[{'before':a,'after':b,'genes':n} for (a,b),n in sorted(transition.items())]
 results['origins'][origin]=entry
(root/'results.json').write_text(json.dumps(results,indent=2,sort_keys=True)+'\n')
print(json.dumps({o:{'genes':r['repaired']['genes'],'states':r['repaired']['fitStates'],'peakNativeRSSBytes':r['repaired']['peakNativeRSSBytes'],'transitions':r.get('statusTransitions')} for o,r in results['origins'].items()}))
