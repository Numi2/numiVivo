#!/usr/bin/env python3
"""Actual CLI structural controls; these do not qualify biological ATAC data."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path
import h5py
import numpy as np
import mudata

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
plan=dict(schemaVersion=1,id='structural',evidence='synthetic',sourceDescription='Synthetic RNA/ATAC structural fixture; no biology',
 sample=dict(id='s',biologicalReplicateID='s',condition='fixture',batchID='s',organism='NCBITaxon:9606'),
 assays=[dict(featureType='Gene Expression',id='rna',kind='rna',featureNamespace='fixture-genes',countUnit='umiCount',genomeAssembly='GRCh38'),
         dict(featureType='Peaks',id='atac',kind='chromatinAccessibility',featureNamespace='fixture-peaks',countUnit='fragmentCount',genomeAssembly='GRCh38')])
(a.out/'plan.json').write_text(json.dumps(plan))
records=[]
def run(name,args,success):
 r=subprocess.run([str(a.binary),*map(str,args)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (a.out/(name+'.log')).write_text(r.stdout)
 records.append(dict(name=name,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=success))
 assert (r.returncode==0)==success,(name,r.stdout)
for case in ['valid','fractional','negative','feature-index','offsets','peak-interval','unmapped-type']:
 path=a.out/(case+'.h5')
 with h5py.File(path,'w') as f:
  g=f.create_group('matrix')
  g['shape']=np.array([2,2],dtype=np.uint64)
  g['barcodes']=np.array(['b0','b1'],dtype=h5py.string_dtype())
  g['indptr']=np.array([0,2,3] if case!='offsets' else [0,3,2],dtype=np.uint64)
  g['indices']=np.array([0,1,1] if case!='feature-index' else [0,9,1],dtype=np.uint64)
  g['data']=np.array([2**53+1,2,3],dtype=np.uint64) if case not in ['fractional','negative'] else np.array([1.5,2,3] if case=='fractional' else [-1,2,3])
  v=g.create_group('features')
  for key,values in dict(id=['gene','chr1:10-20' if case!='peak-interval' else 'chr1:20-10'],name=['Gene','Peak'],
                        feature_type=['Gene Expression','Peaks' if case!='unmapped-type' else 'Unknown'],genome=['GRCh38','GRCh38']).items():
   v[key]=np.array(values,dtype=h5py.string_dtype())
 run(case,['multiassay-10x-import',path,'--plan',a.out/'plan.json','--output',a.out/(case+'-bundle')],case=='valid')
 if case!='valid':assert not (a.out/(case+'-bundle')).exists()
run('verify',['multiassay-verify',a.out/'valid-bundle'],True)
d=json.loads((a.out/'valid-bundle/dataset.json').read_text())
assert d['assays'][0]['matrix']['counts']==[2**53+1]
assert d['assays'][1]['features'][0]['interval']==dict(contig='chr1',start=10,end=20)
assert d['assays'][1]['countUnit']=='fragmentCount'
m=mudata.read_h5mu(a.out/'valid-bundle/dataset.h5mu')
assert int(m.mod['rna'].X[0,0])==2**53+1 and int(m.mod['atac'].X[1,0])==3
result=dict(status='passed',scope='synthetic structural controls only',commands=records,
 integerImportAboveFloatPrecisionExact=True,peakIntervalExact=True,assayUnitsPreserved=True,
 implementationSHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(a.out/'checks.json').write_text(json.dumps(result,indent=2)+'\n')
