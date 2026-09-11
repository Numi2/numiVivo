#!/usr/bin/env python3
"""Reuse the retained complete 10x Visium release; no expression filtering or conversion."""
import argparse,hashlib,json,shutil,tarfile
from pathlib import Path

def sha(p):
 with Path(p).open('rb') as f:
  h=hashlib.sha256()
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--prior',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 source=a.prior/'source/original.h5';spatial=a.prior/'source/spatial.tar.gz'
 assert sha(source)=='eb883e48d9aa935d7959600153ad7ba8c7f1748302fcf56eef4c6999574bb37d'
 assert sha(spatial)=='064f7e62e43a705730c30911911940dbbe583ebc6b743dd2dc52eba5037a2ef5'
 a.out.mkdir(parents=True,exist_ok=False);outs=a.out/'outs';(outs/'spatial').mkdir(parents=True)
 shutil.copy2(source,outs/'filtered_feature_bc_matrix.h5')
 with tarfile.open(spatial) as t:
  m=t.getmember('spatial/tissue_positions_list.csv');assert m.isfile() and m.size<1000000
  positions=t.extractfile(m).read();(outs/'spatial/tissue_positions_list.csv').write_bytes(positions)
 header=b'barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n'
 # Explicit structural control for the documented newer CSV header, using the
 # same original rows. It is not a second biological dataset or newer release.
 (outs/'spatial/tissue_positions.csv').write_bytes(header+positions)
 old=json.loads((a.prior/'plan.json').read_text());mapping=old['assays'][0]
 counts={k:old[k] for k in ['schemaVersion','id','evidence','sourceDescription']};counts['sample']=old['samples'][0]
 counts['assays']=[{k:mapping[k] for k in ['id','kind','featureNamespace','countUnit','genomeAssembly']}];counts['assays'][0]['featureType']='Gene Expression'
 for label,fmt in [('legacy','legacyHeaderlessCSV'),('header','csvWithHeader')]:
  write(a.out/(label+'-plan.json'),dict(schemaVersion=1,counts=counts,positionsFormat=fmt,frame=old['spatial']['frame']))
 write(a.out/'input-freeze.json',dict(sourceArchiveSHA256=sha(spatial),originalCountsSHA256=sha(source),
  files={str(p.relative_to(a.out)):sha(p) for p in sorted(a.out.rglob('*')) if p.is_file()},preparerSHA256=sha(__file__),
  cohort='Every released filtered spot and gene; original legacy CSV plus a header-format structural control'))
 print(json.dumps(dict(status='prepared',countsSHA256=sha(source),positionsBytes=len(positions))))
if __name__=='__main__':main()
