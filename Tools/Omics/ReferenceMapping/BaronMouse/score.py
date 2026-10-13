#!/usr/bin/env python3
"""Score the declared frozen annotation comparison only after all byte seals pass."""
import argparse,csv,datetime,gzip,hashlib,io,json,pathlib,tarfile,math
import numpy as np

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(p.read_text())
def metrics(truth,predicted,training):
 assert len(truth)==len(predicted) and len(truth)>0
 classes=sorted(set(training)|set(truth)|set(predicted));index={v:i for i,v in enumerate(classes)};matrix=np.zeros((len(classes),len(classes)),dtype=np.int64)
 for a,b in zip(truth,predicted):matrix[index[a],index[b]]+=1
 rows=[]
 for c,i in index.items():
  support=int(matrix[i].sum());calls=int(matrix[:,i].sum());tp=int(matrix[i,i]);precision=tp/calls if calls else 0.;recall=tp/support if support else 0.;f1=2*precision*recall/(precision+recall) if precision+recall else 0.
  rows.append(dict(label=c,support=support,predicted=calls,truePositive=tp,precision=precision,recall=recall,f1=f1))
 return dict(cells=len(truth),classes=classes,macroF1=float(np.mean([r['f1'] for r in rows])),accuracy=float(np.trace(matrix)/len(truth)),balancedAccuracy=float(np.mean([r['recall'] for r in rows if r['support']])),perClass=rows,confusionMatrix=matrix.tolist())
def score(study,archive,out):
 assert not out.exists(),'output must be new';protocol=load(study/'protocol-freeze.json');inputs=load(study/'input-freeze.json');freeze=load(study/'prediction-freeze.json')
 assert sha(study/'protocol-freeze.json')==inputs['protocolFreezeSHA256'];assert sha(study/'PROTOCOL.md')==protocol['protocolSHA256'];assert sha(study/'input-freeze.json')==freeze['inputFreezeSHA256'];assert freeze['queryLabelsRead'] is False and freeze['nativeRepeatIdentical'] is True
 assert sha(study/'numivivo-omics')==protocol['implementationSHA256'];assert sha(archive)==protocol['sourceArchiveSHA256']
 for name,h in inputs['files'].items():assert sha(study/'inputs'/name)==h,name
 for name,h in freeze['files'].items():assert sha(study/name)==h,name
 assert sha(study/'mapped-logistic/report.json')==sha(study/'mapped-logistic-repeat/report.json')
 pred=load(study/'predictions.json');model=load(study/'reference-logistic/model.json');mapped=load(study/'mapped-logistic/report.json');knn=load(study/'mapped-knn/report.json');classes=model['classes'];assert classes==sorted(classes)==pred['classes'];assert len(pred['cells'])==1064
 # These arithmetic checks occur before opening the outcome-bearing source member.
 x=np.asarray([r['scores'] for r in mapped['cells']]);fit=model['logistic'];params=np.asarray(fit['parameters']);z=(x-np.asarray(fit['means']))/np.asarray(fit['scales']);assert params.shape==(len(classes),x.shape[1]+1) and np.isfinite(z).all()
 # Native coefficients are class-major with intercept LAST; reproduce its ordered scalar additions.
 logits=[]
 for point in z:
  values=[]
  for row in params:
   value=float(row[-1])
   for j in range(len(point)):value+=float(row[j])*float(point[j])
   values.append(value)
  logits.append(values)
 logits=np.asarray(logits);exp=np.exp(logits-logits.max(axis=1,keepdims=True));prob=exp/exp.sum(axis=1,keepdims=True);stored=np.asarray([r['classProbabilities'] for r in mapped['cells']]);error=float(np.max(np.abs(prob-stored)));assert error<=1e-12
 logistic=[classes[i] if row['totalCounts'] else '__unmapped__' for i,row in zip(prob.argmax(axis=1),mapped['cells'])];assert logistic==pred['nativeLogistic']
 train=np.asarray(model['reduction']['scores']);labels=np.asarray(model['referenceLabels']);centroids=np.asarray([train[labels==c].mean(axis=0) for c in classes]);dist=((x[:,None,:]-centroids[None,:,:])**2).sum(axis=2);centroid=[classes[i] for i in dist.argmin(axis=1)];majority=sorted(classes,key=lambda c:(-int(np.count_nonzero(labels==c)),c))[0]
 assert centroid==pred['nearestCentroid'];assert [majority]*len(x)==pred['trainingMajority'];assert [r.get('candidateLabel') or '__unmapped__' for r in knn['cells']]==pred['nativeKNN15'];assert [r['cell'] for r in mapped['cells']]==[r['cell'] for r in knn['cells']]==pred['cells']
 revealUTC=datetime.datetime.now(datetime.timezone.utc).isoformat();query=next(r for r in inputs['sources'] if r['role']=='query');truth=[];identities=[]
 with tarfile.open(archive) as t:
  raw=t.extractfile(query['member']).read();assert hashlib.sha256(raw).hexdigest()==query['sha256']
  with gzip.GzipFile(fileobj=io.BytesIO(raw)) as src:
   reader=csv.reader(io.TextIOWrapper(src));header=next(reader);assert len(header[3:])==14878
   for row in reader:identities.append(dict(barcode='mouse2:'+row[0],sampleID='mouse2'));truth.append(row[2])
 assert identities==pred['cells'] and len(truth)==query['cells']==1064
 result={name:metrics(truth,pred[name],classes) for name in ['nativeLogistic','nativeKNN15','nearestCentroid','trainingMajority']};m={k:v['macroF1'] for k,v in result.items()};gates={'logistic_minus_centroid_at_least_0_02':m['nativeLogistic']-m['nearestCentroid']>=.02,'logistic_minus_knn_at_least_0_02':m['nativeLogistic']-m['nativeKNN15']>=.02,'logistic_strictly_above_majority':m['nativeLogistic']>m['trainingMajority']}
 record=dict(schema='numi.vivo.frozen-mouse-annotation-result.v1',status='PASS_declared_comparison' if all(gates.values()) else 'FAIL_declared_comparison',revealUTC=revealUTC,protocolSHA256=protocol['protocolSHA256'],predictionFreezeSHA256=sha(study/'prediction-freeze.json'),inputFreezeSHA256=sha(study/'input-freeze.json'),runtimeSHA256=protocol['implementationSHA256'],sourceArchiveSHA256=protocol['sourceArchiveSHA256'],queryMemberSHA256=query['sha256'],modelChangedAfterReveal=False,independentSoftmaxMaximumAbsoluteError=error,nativeRepeatByteIdentical=True,results=result,gates=gates,scope='One historical observational mouse sample/strain holdout; author annotation agreement only. Original CPU execution, not current full-package or Metal qualification; no independent biology or prospective claim.')
 # Recheck all seals after scoring. Nothing in the frozen study is rewritten.
 for name,h in freeze['files'].items():assert sha(study/name)==h,name
 out.mkdir(parents=True);(out/'result.json').write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n');print(json.dumps({'status':record['status'],'macroF1':m,'gates':gates,'independentSoftmaxMaximumAbsoluteError':error}))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=pathlib.Path,required=True);p.add_argument('--archive',type=pathlib.Path,required=True);p.add_argument('--out',type=pathlib.Path,required=True);a=p.parse_args();score(a.study,a.archive,a.out)
