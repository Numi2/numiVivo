from pathlib import Path
import sys,json,collections,gzip
r=Path(__file__).parent;sys.path.insert(0,str(r/'deps'))
import indexed_gzip,h5py,numpy as np
with indexed_gzip.IndexedGzipFile(str(r/'GSE181897_concat.4.raw.h5ad.gz'),index_file=str(r/'source.gzidx')) as gz,h5py.File(gz,'r') as f:
 def column(path):
  d=f[path]
  if 'categories' in d.attrs:
   cats=f[d.attrs['categories']].asstr()[:].tolist();codes=d[:];assert np.all((codes>=-1)&(codes<len(cats)))
   return [None if i<0 else cats[i] for i in codes]
  if d.dtype.kind in ['O','S','U']:return d.asstr()[:].tolist()
  return d[:].tolist()
 obs={k:column('obs/'+k) for k in ['_index','DROPLET.TYPE','batch','cond','exp_id','free_id','pool_code','ct1','ct2','ct3']}
 var={k:column('var/'+k) for k in ['_index','gene_ids','feature_types','genome','batch']}
 for key,d in [('obs',obs),('var',var)]:
  with gzip.open(r/(key+'-metadata.json.gz'),'wt') as o:json.dump(d,o)
  for k,v in d.items():
   if len(set(v))<100:print(key,k,collections.Counter(v))
   else:print(key,k,'unique',len(set(v)),'head',v[:5],'tail',v[-100:-85])
 print('feature_examples',json.dumps([{k:var[k][i] for k in var} for i in list(range(5))+list(range(len(var['_index'])-95,len(var['_index'])))]))
