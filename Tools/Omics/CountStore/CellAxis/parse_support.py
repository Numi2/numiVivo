"""Bounded readers for the exact already-qualified Parse preparation, not general H5AD parsing."""
from pathlib import Path
import hashlib,json,struct,tarfile,zipfile
import numpy as np

PREPARATION_SHA='f08c78ca7a74ab01eab6c3f7748615f45e217cbf95446ba35ef830eaf9c75d02'
COUNTS_SHA='ae94b3d586d2212fc09b1d451f893e7b3badc05ad57eb04550bece351af12b97'
PLAN_SHA='8958f33a675211c355c119a68503cf7cff07575e4dcbff0cd2bd0d2124c36362'
CELLS=725031;FEATURES=40352;RECORDS=1373870697
ROW=struct.Struct('<QIIIIIIQ');QC=struct.Struct('<QQII')
def sha(path):
 digest=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(1048576),b''):digest.update(chunk)
 return digest.hexdigest()
def write(path,value):Path(path).write_text(json.dumps(value,indent=2)+'\n')
def archive_members(path,expected):
 assert sha(path)==expected
 with tarfile.open(path,'r:gz') as t:return json.load(t.extractfile('members.json'))
def bind_source(source,repository):
 preparation=repository/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-preparation/preparation.tar.gz'
 counts=repository/'Tools/Omics/PerturbationPrediction/ParseIFNB/evidence/2026-09-11-counts/results.tar.gz'
 before=archive_members(preparation,PREPARATION_SHA);after=archive_members(counts,COUNTS_SHA)
 for name in ['prepared/plan.json','prepared/runs.json','sources/asset-head.json']:
  assert sha(source/name)==before[name]['SHA256'],name
 assert sha(source/'prepared/plan.json')==PLAN_SHA
 assert sha(source/'donor-complete-counts.npz')==after['donor-complete-counts.npz']['SHA256']
 for run in range(3456):
  for name in [f'axes/run-{run:04d}.npz',f'chunks/run-{run:04d}.json']:
   assert sha(source/name)==before[name]['SHA256'],name
 return dict(preparationArchiveSHA256=PREPARATION_SHA,countArchiveSHA256=COUNTS_SHA,planSHA256=PLAN_SHA,
             verifiedSourceFiles=6916,independentTotalsSHA256=after['donor-complete-counts.npz']['SHA256'])
def metadata_without_cells(plan):
 # This exact JSON ordering is pinned by PLAN_SHA. Read only the dictionary
 # prefix; do not instantiate the complete cell array for preparation.
 marker='"cells":[';prefix=''
 with plan.open() as f:
  while marker not in prefix:
   chunk=f.read(65536);assert chunk;prefix+=chunk;assert len(prefix.encode())<=67108864
 prefix=prefix[:prefix.index(marker)]+'"cells":[]}}'
 value=json.loads(prefix);assert set(value)=={'schemaVersion','metadata'} and value['schemaVersion']==1
 assert value['metadata']['cells']==[]
 return value['metadata']
def array_items(plan,name,maximum_item_chars=16384):
 marker='"'+name+'":[';decoder=json.JSONDecoder()
 with plan.open() as f:
  buffer=''
  while marker not in buffer:
   chunk=f.read(65536);assert chunk;buffer=buffer[-len(marker):]+chunk
  buffer=buffer[buffer.index(marker)+len(marker):]
  while True:
   buffer=buffer.lstrip()
   if not buffer:
    buffer=f.read(65536);assert buffer;continue
   if buffer[0]==']':return
   while True:
    try:
     value,end=decoder.raw_decode(buffer)
     # A scalar at a chunk boundary may have more digits in the next chunk.
     if end==len(buffer):
      chunk=f.read(65536);assert chunk;buffer+=chunk;continue
     break
    except json.JSONDecodeError:
     assert len(buffer)<=maximum_item_chars
     chunk=f.read(65536);assert chunk;buffer+=chunk
   yield value
   buffer=buffer[end:].lstrip()
   if not buffer:buffer=f.read(65536).lstrip();assert buffer
   if buffer[0]==']':return
   assert buffer[0]==',';buffer=buffer[1:]
def source_rows(source):
 runs=json.loads((source/'prepared/runs.json').read_text());cells=array_items(source/'prepared/plan.json','cells');cardinality=array_items(source/'prepared/plan.json','rowNonzeros')
 with zipfile.ZipFile(source/'donor-complete-counts.npz') as z,z.open('cellTotals.npy') as totals:
  version=np.lib.format.read_magic(totals);assert version in [(1,0),(2,0)]
  read_header=np.lib.format.read_array_header_1_0 if version==(1,0) else np.lib.format.read_array_header_2_0
  shape,fortran,dtype=read_header(totals);assert shape==(CELLS,) and not fortran and dtype==np.dtype('<u8')
  count=entries=total_count=0
  for run in runs:
   with np.load(source/f'axes/run-{run["run"]:04d}.npz') as a:
    assert int(a['lo'])==run['lo']==count and int(a['hi'])==run['hi']
    barcodes=a['barcodes'];nnz=np.diff(a['indptr']);raw=totals.read(8*len(barcodes));assert len(raw)==8*len(barcodes)
    row_totals=np.frombuffer(raw,dtype='<u8')
    for barcode,n,total in zip(barcodes,nnz,row_totals):
     cell=next(cells);assert set(cell)=={'barcode','sampleID'} and cell==dict(barcode=str(barcode),sampleID=run['sample'])
     assert next(cardinality)==int(n)
     yield dict(barcode=cell['barcode'],sampleID=cell['sampleID'],nonzeros=int(n),totalCounts=int(total))
     count+=1;entries+=int(n);total_count+=int(total)
  assert count==CELLS and entries==RECORDS and total_count==3070817047
  assert totals.read(1)==b'' and next(cells,None) is None and next(cardinality,None) is None

def check_axis(source,axis,header):
 saved=json.loads((axis/'header.json').read_text());assert saved==header
 receipt=json.loads((axis/'receipt.json').read_text())
 for name,key in [('header.json','header'),('rows.bin','rows'),('strings.bin','strings')]:assert sha(axis/name)==bytes(receipt[key]['bytes']).hex()
 expected_offset=0;checked=0
 with (axis/'rows.bin').open('rb') as rows,(axis/'strings.bin').open('rb') as strings:
  for expected in source_rows(source):
   raw=rows.read(ROW.size);assert len(raw)==ROW.size
   offset,barcode_length,sample_length,group_length,sample,annotation,nnz,total=ROW.unpack(raw)
   assert offset==expected_offset and group_length==0 and annotation==2**32-1
   barcode=strings.read(barcode_length);sample_id=strings.read(sample_length)
   assert barcode==expected['barcode'].encode() and sample_id==expected['sampleID'].encode()
   assert header['metadata']['samples'][sample]['id']==expected['sampleID']
   assert nnz==expected['nonzeros'] and total==expected['totalCounts']
   expected_offset+=barcode_length+sample_length;checked+=1
  assert rows.read(1)==b'' and strings.read(1)==b'' and expected_offset==receipt['stringBytes']
 assert checked==CELLS
 return dict(status='passed',cells=checked,allOriginalIdentityBytesAndRowsExact=True,allDeclaredMatrixTotalsAndCardinalitiesExact=True,
             rowBytes=ROW.size,axisRowsBytes=(axis/'rows.bin').stat().st_size,axisStringsBytes=expected_offset,predictionFitOrScoring=False)

def object_value(path,name,maximum_chars=67108864):
 marker='"'+name+'":';decoder=json.JSONDecoder()
 with path.open() as f:
  buffer=''
  while marker not in buffer:
   chunk=f.read(65536);assert chunk;buffer=buffer[-len(marker):]+chunk
  buffer=buffer[buffer.index(marker)+len(marker):].lstrip()
  while True:
   try:return decoder.raw_decode(buffer)[0]
   except json.JSONDecodeError:
    assert len(buffer)<=maximum_chars
    chunk=f.read(65536);assert chunk;buffer+=chunk
