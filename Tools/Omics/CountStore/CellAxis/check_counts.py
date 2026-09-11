"""Offline complete-cohort checker; does not replay native code or fetch source ranges."""
from pathlib import Path
import argparse,importlib.util,json,sys

def check(source,work,repository):
 # Load the exact recipe that executed the count run, not a mutable installed helper.
 freeze=json.loads((work/'count-recipe/source-freeze.json').read_text())
 sys.path.insert(0,str(work/'count-recipe'))
 from parse_support import sha,bind_source,check_axis,archive_members,CELLS,FEATURES,RECORDS,COUNTS_SHA
 for item in freeze['recipe']:assert sha(work/'count-recipe'/item['path'])==item['SHA256'],item['path']
 from compare_parse_counts import compare
 pipeline=json.loads((work/'count-pipeline.json').read_text())
 assert pipeline['status']=='passed-paired-full-count-ingestion-and-replay','Both complete native phases are required'
 assert freeze['unitTests']==14 and freeze['unitSuites']==2 and len(freeze['nativeBuildInputs'])==105
 binary=sha(work/'build/numivivo-omics');assert binary==freeze['binarySHA256']
 binding=bind_source(source,repository)
 axis=json.loads((work/'axis-execution.json').read_text());assert axis['status']=='passed-complete-cell-axis'
 assert axis['binarySHA256']==binary and axis['sourceBindings']==binding
 assert check_axis(source,work/'axis',json.loads((work/'header.json').read_text()))==axis['independentCheck']
 for name in ['header.json','receipt.json','rows.bin','strings.bin']:
  assert sha(work/'axis'/name)==sha(work/'native-file/axis'/name),name
 file_receipt=json.loads((work/'native-file/receipt.json').read_text())
 assert sha(work/'native-file/axis/receipt.json')==bytes(file_receipt['axis']['bytes']).hex()
 assert file_receipt['qualityEncoding']=='total-u64-mito-u64-nnz-u32-group-u32-le-v1'
 prior=archive_members(repository/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-counts/results.tar.gz',COUNTS_SHA)
 expected={}
 for name,item in prior.items():
  p=Path(name)
  if name.startswith('donor-execution/ingest/') and p.name.startswith('run-') and p.suffix=='.json':
   index=int(p.stem[4:]);assert index not in expected;expected[index]=item
 assert sorted(expected)==list(range(3456))
 # Read archived reference reports sequentially, without requiring the removed donor directories.
 import tarfile,hashlib
 selected={v['SHA256']:i for i,v in expected.items()};references={}
 with tarfile.open(repository/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-counts/results.tar.gz','r|gz') as container:
  for item in container:
   digest=item.name.removeprefix('objects/')
   if digest not in selected:continue
   assert item.size<=1048576;raw=container.extractfile(item).read();assert hashlib.sha256(raw).hexdigest()==digest
   references[selected[digest]]=json.loads(raw)
 assert len(references)==3456
 fields=['nonzeros','totalCounts','sourceGeneCountMismatch','sourceTscpCountMismatch','sourceMinusMatrixTotal','sourceMinusMatrixGeneCount','sortedRows']
 phases={};stream=None;results={}
 for phase in ['ingest','verify']:
  directory=work/phase;execution=json.loads((directory/'execution.json').read_text())
  assert execution['status']=='passed' and execution['phase']==phase
  assert execution['cells']==CELLS and execution['records']==RECORDS and execution['streamBytes']==RECORDS*16
  assert execution['sourceBindings']==binding and execution['predictionFitOrScoring'] is False
  assert execution['binaries']==dict(file=binary,resident='20b13e585e1526dba7be566484373a3e4775efebd66f22b3b20dff31c1ee3516')
  if stream is None:stream=execution['streamSHA256']
  assert stream==execution['streamSHA256']==bytes(file_receipt['stream']['bytes']).hex()
  assert {f.name for f in directory.glob('run-*.json')}=={f'run-{i:04d}.json' for i in range(3456)}
  totals={key:0 for key in fields};ranges=errors=0
  for i in range(3456):
   observed=json.loads((directory/f'run-{i:04d}.json').read_text());before=references[i]
   assert observed['run']==i and all(observed[key]==before[key] for key in fields),(phase,i)
   assert [(r['start'],r['stop'],r['SHA256']) for r in observed['ranges']]==[(r['start'],r['stop'],r['SHA256']) for r in before['ranges']],(phase,i)
   for key in fields:totals[key]+=observed[key]
   ranges+=len(observed['ranges']);errors+=sum(len(r['transportErrors']) for r in observed['ranges'])
  assert totals['nonzeros']==RECORDS and totals['totalCounts']==3070817047
  metrics=execution['native'];assert set(metrics)=={'file','resident'}
  for owner in metrics:
   result=metrics[owner];assert result['nativeExit']==0 and result['bytesSent']==RECORDS*16 and result['maximumNativeRSS']>0 and result['seconds']>0
   saved=json.loads((work/('native-'+owner)/'receipt.json').read_text())
   assert json.loads((directory/(owner+'.stdout')).read_text())==saved
  independent=compare(source,work,directory/'independent-counts.npz')
  assert independent==execution['independentCheck']==json.loads((directory/'independent-check.json').read_text())
  phases[phase]=dict(native=metrics,ranges=ranges,recoveredTransportErrors=errors,sourceQC=totals,allPriorRangeIdentitiesExact=True,independentCheck=independent)
  results[phase]=dict(filePeakRSSBytes=metrics['file']['maximumNativeRSS'],residentPeakRSSBytes=metrics['resident']['maximumNativeRSS'],fileToResidentPeakRSSRatio=metrics['file']['maximumNativeRSS']/metrics['resident']['maximumNativeRSS'])
 return dict(status='passed-complete-paired-count-evidence',cells=CELLS,features=FEATURES,groups=24,records=RECORDS,streamBytes=RECORDS*16,streamSHA256=stream,phases=phases,memoryComparison=results,sourceBindings=binding,binarySHA256=binary,
  measurementScope='Peak RSS of each native process, same complete canonical source stream; excludes Python adapter and OS file cache. Elapsed times include upstream streaming waits and cannot measure isolated throughput or a speedup.',
  offlineCheckScope='Checks retained complete native execution records, all original range identities, all QC/memberships and independent aggregates. This checker does not execute native code or replay the network source.',predictionFitOrScoring=False)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--source-study',type=Path,required=True);parser.add_argument('--work',type=Path,required=True);parser.add_argument('--dependencies-repo',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
 assert not args.out.exists()
 def forbid_network_and_children(event,args):
  if event.startswith('socket.') or event in ('subprocess.Popen','os.system','os.posix_spawn','os.fork'):raise RuntimeError('Offline evidence checker forbids network and subprocess operations')
 sys.addaudithook(forbid_network_and_children)
 result=check(args.source_study.resolve(),args.work.resolve(),args.dependencies_repo.resolve())
 args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
