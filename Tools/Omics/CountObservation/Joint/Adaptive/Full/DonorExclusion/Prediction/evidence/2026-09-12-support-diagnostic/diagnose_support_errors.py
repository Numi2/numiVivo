from pathlib import Path
import gzip,json,hashlib,math
root=Path(__file__).parent;out={'scope':'Query-control pseudobulk CPM relative to training donor NB MLE rate range; descriptive mismatch, not an NB MLE comparison or held-out outcome evaluation. Uses frozen prediction inputs and held-out scores after evaluation; exploratory association only.','folds':{}}
for tag in [f'Kang-{i:02d}' for i in range(8)]+[f'HIRISA-{i:02d}' for i in range(5)]:
 n=low=high=0;bindings={}
 scorefile=root/tag/'scores-by-gene.json.gz'; report=json.loads((root/tag/'scores.json').read_text());assert hashlib.sha256(scorefile.read_bytes()).hexdigest()==report['scoresSHA256']; scores={r['featureIndex']:r for r in json.loads(gzip.decompress(scorefile.read_bytes()))}; groups={'inside':[], 'outside':[]}
 for p in sorted((root/tag).glob('*input.jsonl.gz')):
  bindings[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
  with gzip.open(p,'rt') as f:
   h=json.loads(next(f));total=sum(a*b for a,b in zip(h['cellsPerLibrary'],h['libraryCounts']));assert total>0
   for l in f:
    x=json.loads(l)
    if x['status']!='boundedContinuousLikelihood':continue
    m=x['model']['model'];assert h['queryDonorID'] not in m['trainingDonorIDs'];q=sum(x['counts'])*1e6/total;rates=m['controlDonorMLERatesCPM'];n+=1;low+=q<min(rates);high+=q>max(rates)
    row=scores[x['featureIndex']];assert math.isclose(math.log1p(q),row['controlLog1pCPM'],abs_tol=1e-10);groups['outside' if q<min(rates) or q>max(rates) else 'inside'].append(row)
 metrics={}
 for group,rows in groups.items():
  metrics[group]={'genes':len(rows),'jointRMSE':math.sqrt(math.fsum((r['predictionLog1pCPM']-r['observedLog1pCPM'])**2 for r in rows)/len(rows)),'meanRMSE':math.sqrt(math.fsum((r['trainingMeanResponseLog1pCPM']-r['observedLog1pCPM'])**2 for r in rows)/len(rows))}
 out['folds'][tag]={'eligibleGenes':n,'belowTrainingRange':low,'aboveTrainingRange':high,'fractionOutside':(low+high)/n,'inputSHA256':bindings,'postHocErrorStrata':metrics,'scoresSHA256':report['scoresSHA256']};print(tag,round((low+high)/n,4),flush=True)
out['status']='complete';(root/'support-error-diagnosis.json').write_text(json.dumps(out,indent=2)+'\n')
