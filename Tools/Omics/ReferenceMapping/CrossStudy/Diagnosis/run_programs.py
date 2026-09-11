#!/usr/bin/env python3
"""Run native transcript-program diagnostics on both complete frozen sources."""
import argparse,json,os,re,subprocess,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'Logistic'))
from prepare import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 programs=[]
 for name,genes in [('hemoglobin-three',['HBB','HBA1','HBA2']),('pf4-ppbp-five',['PF4','PPBP','ITGA2B','GP9','TUBB1'])]:
  programs.append(dict(id=name,organism='NCBITaxon:9606',featureNamespace='literal-source-gene-symbol-v1',sourceURI='urn:numivivo:post-result-reference-transfer-diagnosis:20260911',sourceVersion='1',sourceDescription='Analyst-selected post-result transcript means for tracing the frozen source-label discrepancy. Not independently validated cell-identity or protein markers.',members=[dict(featureID=g,weight=1) for g in genes],minimumWeightCoverage=1))
 records={}
 for cohort,file in [('kang','kang.h5ad'),('ding','ding-query.h5ad')]:
  source=a.prior/'inputs'/file;mapping=json.loads((a.prior/'inputs'/(cohort+'-query.json')).read_text())['mapping'];plan=dict(schemaVersion=1,mapping=mapping,programs=dict(definitions=programs),normalizationTarget=10000,matchFeatureNames=False);write(a.out/(cohort+'-plan.json'),plan);records[cohort]=dict(sourcePath=str(source),sourceSHA256=sha(source),planSHA256=sha(a.out/(cohort+'-plan.json')))
 write(a.out/'input-freeze.json',dict(protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),runnerSHA256=sha(__file__),binarySHA256=sha(a.binary),hdf5SHA256=sha(os.environ['NUMIVIVO_HDF5_LIBRARY']),sources=records,modelRefitting=False,sourceRelabelling=False))
 commands=[]
 for cohort in ['kang','ding']:
  for label,args in [('score',['singlecell-h5ad-programs',records[cohort]['sourcePath'],'--plan',a.out/(cohort+'-plan.json'),'--output',a.out/cohort]),('verify',['singlecell-h5ad-programs-verify',a.out/cohort])]:
   t=time.monotonic();r=subprocess.run(['/usr/bin/time','-l',str(a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(cohort+'-'+label+'.log')).write_text(r.stdout+r.stderr);peak=re.search(r'(\d+)\s+maximum resident set size',r.stderr);commands.append(dict(cohort=cohort,operation=label,args=list(map(str,args)),exitCode=r.returncode,seconds=time.monotonic()-t,maximumResidentBytes=int(peak[1]) if peak else None));write(a.out/'commands.json',commands);assert r.returncode==0,r.stderr
 write(a.out/'execution.json',dict(status='completed',cells=68704,commands=commands,binarySHA256=sha(a.binary),inputsSHA256=sha(a.out/'input-freeze.json')));print(json.dumps(dict(status='completed',commands=len(commands),cells=68704)))
if __name__=='__main__':main()
