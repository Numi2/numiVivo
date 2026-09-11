"""Complete paired Parse aggregate comparisons with frozen Bioconductor workflows."""
import argparse,csv,hashlib,json,os,platform,subprocess,time
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
 return h.hexdigest()
def write(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
def table(path,header,rows):
 with path.open('x',newline='') as f:
  writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(header);writer.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--rscript',type=Path,required=True);p.add_argument('--r-library',type=Path,required=True);a=p.parse_args();s=a.study.resolve()
 state=s/'reference-pipeline.json';assert not state.exists();freeze=json.loads((s/'reference-freeze.json').read_text())
 for path,digest in freeze['files'].items():assert sha(Path(path))==digest,path
 protocol=json.loads((s/'protocol.json').read_text());native_state=json.loads((s/'native-pipeline.json').read_text());expected_pid=native_state['pid']
 def owner():return subprocess.check_output(['ps','-p',str(expected_pid),'-o','lstart=','-o','command='],text=True).strip()
 identity=owner() if native_state['status'] in ('starting','native-comparison-running','checking-native-results') else None
 if identity is not None:assert str(s/'recipe/run_parse_native.py') in identity
 record=dict(status='waiting-for-native-comparison',pid=os.getpid(),nativeControllerPID=expected_pid,startedUnix=time.time(),predictionFitOrScoring=False);write(state,record)
 try:
  while native_state['status'] in ('starting','native-comparison-running','checking-native-results'):
   assert identity is not None and owner()==identity,'Native controller missing or changed; no automatic retry'
   record['lastNativeStatus']=native_state['status'];record['lastLiveObservationUnix']=time.time();write(state,record);time.sleep(30)
   try:native_state=json.loads((s/'native-pipeline.json').read_text())
   except json.JSONDecodeError:time.sleep(0.1);native_state=json.loads((s/'native-pipeline.json').read_text())
  assert native_state['status']=='passed-full-parse-native-aggregate-handoff',native_state
  record.update(status='preparing-complete-reference-inputs');write(state,record)
  native=json.loads((s/'native-file/report.json').read_text());counts=Path(protocol['countStudy'])/'native-file'
  assert sha(counts/'receipt.json')==protocol['sourceReceiptSHA256']['file']
  source=json.loads((counts/'report.json').read_text());features=source['featureIDs'];groups=source['groups'];raw=source['matrix']
  assert len(features)==40352 and len(groups)==24 and len(set(features))==40352
  matrix=np.zeros((24,40352),dtype=np.uint64)
  for row in range(24):
   start,stop=raw['rowOffsets'][row:row+2];matrix[row,np.asarray(raw['featureIndices'][start:stop],dtype=np.int64)]=raw['counts'][start:stop]
  design=native['design'];indices=design['sourcePseudobulkIndices'];observations=design['observations'];request=protocol['request'];n=len(observations)
  assert indices==list(range(24)) and n==24 and all(g['sourceCellCount']>=request['minimumCellsPerPseudobulk'] for g in groups)
  donors=sorted({g['donorID'] for g in observations});assert len(donors)==12
  assert all({g['condition'] for g in observations if g['donorID']==donor}=={'PBS','IFN-beta'} for donor in donors)
  assert {batch for g in observations for batch in g['batchIDs']}=={'Parse-original-release'}
  independent_design=[[1,1 if g['condition']=='IFN-beta' else 0]+[1 if g['donorID']==donor else 0 for donor in donors[1:]] for g in observations]
  assert independent_design==design['rows'] and design['columnNames']==['intercept','treatment-minus-control']+['donor:'+d for d in donors[1:]]
  assert design['contrast']==[0,1]+[0]*11 and design['residualDegreesOfFreedom']==11
  libraries=matrix.sum(axis=1);assert libraries.tolist()==design['libraryCounts'] and int(libraries.sum())==3070817047
  positive=np.all(matrix>0,axis=0);assert np.flatnonzero(positive).tolist()==design['referenceFeatureIndices']
  logged=np.log(matrix[:,positive].astype(np.float64));ratios=np.exp(logged-logged.mean(axis=0));log_factors=np.log(np.median(ratios,axis=1));independent_factors=np.exp(log_factors-log_factors.mean())
  assert np.allclose(independent_factors,design['sizeFactorValues'],rtol=1e-12,atol=0)
  assert all(len(g['sampleIDs'])==1 for g in observations);samples=[g['sampleIDs'][0] for g in observations]
  keep=(matrix.sum(axis=0)>=request['minimumFeatureCounts'])&((matrix>0).sum(axis=0)>=request['minimumExpressingPseudobulks'])
  inputs=s/'reference-inputs';inputs.mkdir(exist_ok=False)
  table(inputs/'counts.tsv',['featureID']+samples,([feature]+[str(x) for x in matrix[:,i]] for i,feature in enumerate(features)))
  table(inputs/'design.tsv',['sampleID']+design['columnNames'],([sample]+design['rows'][i] for i,sample in enumerate(samples)))
  table(inputs/'samples.tsv',['sampleID','donor','condition','sizeFactor','libraryCounts'],([sample,observations[i]['donorID'],observations[i]['condition'],format(design['sizeFactorValues'][i],'.17g'),str(libraries[i])] for i,sample in enumerate(samples)))
  metadata=dict(scope='Complete original Parse paired PBMC aggregates; descriptive inference comparison, not prediction or calibrated FDR',minimumFeatureCounts=request['minimumFeatureCounts'],minimumExpressingPseudobulks=request['minimumExpressingPseudobulks'],eligibleFeatures=int(keep.sum()),contrast=design['contrast'],files={p.name:sha(p) for p in inputs.iterdir()})
  write(inputs/'input.json',metadata);write(inputs/'independent-design-check.json',dict(status='passed',cells=725031,features=40352,donors=12,observations=24,completeCounts=3070817047,allCountsRetained=True,allPairedDesignValuesExact=True,independentSizeFactorsRelativeTolerance=1e-12,maximumSizeFactorRelativeError=float(np.max(np.abs(independent_factors/np.array(design['sizeFactorValues'])-1))),predictionFitOrScoring=False))
  command=[str(a.rscript),str(s/'recipe/reference.R'),str(inputs),str(s/'reference')]
  record.update(status='running-bioconductor',command=command,inputSHA256=sha(inputs/'input.json'),platform=platform.platform());write(state,record);started=time.monotonic()
  with (s/'reference.log').open('x') as log:
   result=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,R_LIBS_USER=str(a.r_library),VECLIB_MAXIMUM_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
  record.update(referenceExit=result.returncode,referenceSeconds=time.monotonic()-started,status='comparing-complete-reference-results');write(state,record)
  assert result.returncode==0,'Bioconductor failure retained; do not change the frozen comparison to force a pass'
  runs=json.loads((s/'reference/runs.json').read_text());assert len(runs)==6 and all(r['status']=='completed' for r in runs.values())
  native_features=native['features'];assert [f['featureID'] for f in native_features]==features
  summaries={};reference_rows={}
  def number(value):
   if value in (None,'NA','NaN',''):return float('nan')
   return float(value)
  for method in runs:
   with (s/'reference'/(method+'.tsv')).open() as f:rows={r['featureID']:r for r in csv.DictReader(f,delimiter='\t')}
   assert set(rows)=={feature for i,feature in enumerate(features) if keep[i]};reference_rows[method]=rows
   finite=[];p_pairs=[];native_sig=set();reference_sig=set()
   for i,f in enumerate(native_features):
    if not keep[i]:continue
    other=rows[f['featureID']];x=number(f.get('log2FoldChange'));y=number(other['log2FoldChange'])
    if np.isfinite(x) and np.isfinite(y):finite.append((x,y))
    x_p=number(f.get('pValue'));y_p=number(other['pValue'])
    if np.isfinite(x_p) and np.isfinite(y_p):p_pairs.append((x_p,y_p))
    if number(f.get('adjustedPValue'))<=0.05:native_sig.add(f['featureID'])
    if number(other['adjustedPValue'])<=0.05:reference_sig.add(f['featureID'])
   effects=np.array(finite);ps=np.array(p_pairs);assert len(effects)>1 and len(ps)>1
   union=native_sig|reference_sig
   summaries[method]=dict(referenceEligibleFeatures=len(rows),finiteEffectPairs=len(effects),effectRMSE=float(np.sqrt(np.mean((effects[:,0]-effects[:,1])**2))),effectMedianAbsoluteDifference=float(np.median(np.abs(effects[:,0]-effects[:,1]))),effectMaximumAbsoluteDifference=float(np.max(np.abs(effects[:,0]-effects[:,1]))),effectSpearman=float(spearmanr(effects[:,0],effects[:,1]).statistic),effectSignConcordance=float(np.mean(np.sign(effects[:,0])==np.sign(effects[:,1]))),finitePValuePairs=len(ps),pValueSpearman=float(spearmanr(ps[:,0],ps[:,1]).statistic),nativeBH005=len(native_sig),referenceBH005=len(reference_sig),sharedBH005=len(native_sig&reference_sig),BH005Jaccard=len(native_sig&reference_sig)/len(union) if union else None)
  columns=['featureID','nativeStatus','nativeLog2FoldChange','nativePValue','nativeAdjustedPValue']
  for method in runs:columns.extend(method+'-'+key for key in ['log2FoldChange','pValue','adjustedPValue'])
  def joined():
   for f in native_features:
    row=[f['featureID'],f['status'],f.get('log2FoldChange'),f.get('pValue'),f.get('adjustedPValue')]
    for method in runs:
     other=reference_rows[method].get(f['featureID'],{});row.extend(other.get(key) for key in ['log2FoldChange','pValue','adjustedPValue'])
    yield ['NA' if value is None else value for value in row]
  table(s/'all-feature-reference-comparison.tsv',columns,joined())
  statuses={}
  for f in native_features:statuses[f['status']]=statuses.get(f['status'],0)+1
  checked=dict(status='completed-six-full-parse-reference-comparisons',cells=725031,features=40352,donors=12,observations=24,commonCountFilterEligibleFeatures=int(keep.sum()),nativeStatuses=statuses,comparisons=summaries,runs=runs,policy='All original aggregate observations and features; common declared low-count filter. Shared native offsets and package-specific normalization are separate. Unlike Wald, QL F and moderated-t tests are not required to give identical p-values. Concordance is descriptive, not FDR calibration or prediction validation.',predictionFitOrScoring=False)
  write(s/'reference-comparison.json',checked)
  for path,digest in freeze['files'].items():assert sha(Path(path))==digest,path
  for name,digest in metadata['files'].items():assert sha(inputs/name)==digest,name
  record.update(status=checked['status'],finishedUnix=time.time(),comparison=str(s/'reference-comparison.json'));write(state,record);print(json.dumps(checked),flush=True)
 except BaseException as error:
  record.update(status='failed',errorType=type(error).__name__,error=str(error),finishedUnix=time.time());write(state,record);raise

if __name__=='__main__':main()
