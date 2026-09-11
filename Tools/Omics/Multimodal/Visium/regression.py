#!/usr/bin/env python3
"""Exercise native publication, exact replay, tamper rejection and unchanged 10x output."""
import argparse,copy,json,shutil,subprocess
from pathlib import Path
from prepare import sha,write

def main():
 p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
 shutil.copytree(a.inputs,a.out/'inputs');root=a.out/'inputs';commands=[]
 def run(label,args,success=True,binary=None):
  r=subprocess.run([str(binary or a.binary),*map(str,args)],capture_output=True,text=True);(a.out/(label+'.log')).write_text(r.stdout+r.stderr)
  commands.append(dict(label=label,arguments=list(map(str,args)),exitCode=r.returncode,expectedSuccess=success));write(a.out/'commands.json',commands);assert (r.returncode==0)==success,(label,r.stderr)
 def imp(label,plan='plan.json',success=True):
  dest=a.out/label;run(label,['multiassay-visium-import',root/'outs','--plan',root/plan,'--output',dest],success);assert dest.exists()==success;return dest
 base=imp('native');run('verify',['multiassay-verify',base]);repeat=imp('repeat')
 for p in base.iterdir():assert sha(p)==sha(repeat/p.name),p.name
 original=(root/'outs/spatial/tissue_positions_list.csv').read_bytes();(root/'outs/spatial/tissue_positions_list.csv').write_bytes(b'invalid original source now\n')
 run('snapshot-independent',['multiassay-verify',base]);(root/'outs/spatial/tissue_positions_list.csv').write_bytes(original)
 run('existing-destination',['multiassay-visium-import',root/'outs','--plan',root/'plan.json','--output',base],False)
 assert sha(base/'dataset.json')==sha(repeat/'dataset.json')
 with (repeat/'positions.csv').open('ab') as f:f.write(b'changed')
 run('tampered-positions',['multiassay-verify',repeat],False)
 for label,content in [('missing',b'b0,1,0,0,10,20\n'),('extra',original+b'extra,1,0,4,12,22\n'),('duplicate',original+b'b0,1,0,0,10,20\n'),('outside',original.replace(b'b1,1',b'b1,0')),('swapped-format',b'barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n'+original)]:
  (root/'outs/spatial/tissue_positions_list.csv').write_bytes(content);imp(label,success=False)
 (root/'outs/spatial/tissue_positions_list.csv').write_bytes(original)
 plan=json.loads((root/'plan.json').read_text());bad=copy.deepcopy(plan);bad['frame']['unit']='micrometer';write(root/'bad-plan.json',bad);imp('invalid-unit','bad-plan.json',False)
 positions=root/'outs/spatial/tissue_positions_list.csv';positions.rename(root/'positions-retained.csv');positions.symlink_to(root/'positions-retained.csv');imp('symlink-position',success=False);positions.unlink();positions.write_bytes(original)
 write(root/'tenx.json',plan['counts'])
 for label,binary in [('prior-10x',a.baseline),('current-10x',a.binary)]:run(label,['multiassay-10x-import',root/'outs/filtered_feature_bc_matrix.h5','--plan',root/'tenx.json','--output',a.out/label],binary=binary)
 for name in ['dataset.json','dataset.h5mu','plan.json','original.h5']:assert sha(a.out/'prior-10x'/name)==sha(a.out/'current-10x'/name)
 assert 'positions' not in json.loads((a.out/'current-10x/receipt.json').read_text())
 run('current-10x-verify',['multiassay-verify',a.out/'current-10x'])
 assert not list(a.out.glob('.numivivo-multiassay-*'))
 write(a.out/'checks.json',dict(status='passed',commands=len(commands),exactRepeat=True,sourceMutationDoesNotChangeSnapshot=True,positionsTamperRejected=True,rejectedMissingExtraDuplicateOutsideRows=True,ambiguousFormatRejected=True,wrongUnitsRejected=True,symlinkRejected=True,failedOutputsAbsent=True,priorNonspatialBytesExact=True,optionalReceiptFieldAbsentForNonspatial=True,binarySHA256=sha(a.binary),baselineSHA256=sha(a.baseline),checkerSHA256=sha(__file__)))
 print((a.out/'checks.json').read_text())
if __name__=='__main__':main()
