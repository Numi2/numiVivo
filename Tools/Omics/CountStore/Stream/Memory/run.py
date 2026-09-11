from pathlib import Path
import numpy as np,json,subprocess,time,resource,sys,argparse,hashlib
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--label',required=True);a=p.parse_args();assert sys.platform=='darwin', 'RSS units qualified on macOS only';assert a.label.replace('-','').isalnum();root=a.out;root.mkdir(parents=True,exist_ok=True);binary=a.binary.resolve();binary_sha=hashlib.sha256(binary.read_bytes()).hexdigest()
features=[dict(id=f'f{i}',name=f'f{i}',mitochondrial=False) for i in range(100000)]
for rows in [100,1000]:
 out=root/f'{a.label}-{rows}';out.mkdir(exist_ok=False)
 plan=dict(schemaVersion=1,metadata=dict(id='allocator-control',evidence='synthetic',sourceDescription='Deterministic allocator stress control; not biological data',countUnit='umiCount',samples=[dict(id='s',biologicalReplicateID='d',donorID='d',condition='fixture',batchID='fixture',organism='fixture')],features=features,cells=[dict(barcode=f'c{i}',sampleID='s') for i in range(rows)]),sourceDeclaration='Every feature has one count in every synthetic row',rowNonzeros=[100000]*rows,rowTotals=[100000]*rows)
 (out/'plan.json').write_text(json.dumps(plan,separators=(',',':')))
 packed=np.empty((100000,2),dtype='<u8');packed[:,1]=1;base=np.arange(100000,dtype=np.uint64)<<np.uint64(32)
 before=time.monotonic()
 with (out/'stdout.json').open('wb') as stdout,(out/'stderr.log').open('wb') as stderr:
  p=subprocess.Popen([str(binary),'singlecell-count-stream-pseudobulk','--plan',str(out/'plan.json'),'--output',str(out/'native')],stdin=subprocess.PIPE,stdout=stdout,stderr=stderr)
  for row in range(rows):packed[:,0]=base|np.uint64(row);p.stdin.write(packed.tobytes())
  p.stdin.close();_,status,usage=__import__('os').wait4(p.pid,0);p.returncode=__import__('os').waitstatus_to_exitcode(status);assert p.returncode==0
 report=json.loads((out/'native/report.json').read_text());assert report['pseudobulk']['matrix']['counts']==[rows]*100000;assert all(x['totalCounts']==100000 for x in report['quality'])
 summary=dict(binarySHA256=binary_sha,rows=rows,features=100000,records=rows*100000,streamBytes=rows*100000*16,seconds=time.monotonic()-before,maximumRSS=usage.ru_maxrss,nativeExit=0,allCountsVerified=True)
 (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
