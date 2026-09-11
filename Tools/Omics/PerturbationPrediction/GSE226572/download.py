#!/usr/bin/env python3
"""Download every metadata-listed raw source with hashes and a capacity reserve."""
import argparse,concurrent.futures,hashlib,json,re,shutil,urllib.request
from pathlib import Path

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
def main():
 p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);a=p.parse_args();root=a.study;out=root/'sources';out.mkdir(exist_ok=False)
 records=[]
 for block in (root/'series.soft').read_text().split('^SAMPLE = ')[1:]:
  lines=block.splitlines();accession=lines[0];title=next(x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_title = '));m=re.fullmatch(r'(D\d+)_(\d+)h(?:_2)?',title);assert m
  characteristics=[x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_characteristics_ch1 = ')];time=next(x.removeprefix('time of_stimulation: ') for x in characteristics if x.startswith('time of_stimulation: '));assert time==m[2]+'h'
  treated=int(m[2])>0;assert ('treatment: human recombinant IFN-β' in characteristics)==treated
  url=next(x.split(' = ',1)[1] for x in lines if x.startswith('!Sample_supplementary_file_1 = ')).replace('ftp://','https://');assert url.endswith('.h5')
  records.append(dict(accession=accession,title=title,donor=m[1],hours=int(m[2]),condition='IFNB' if treated else 'control',url=url,path=url.rsplit('/',1)[1]))
 assert len(records)==24 and {r['donor'] for r in records}=={'D34','D38','D39'} and sum(r['hours']==0 for r in records)==6
 for r in records:
  with urllib.request.urlopen(urllib.request.Request(r['url'],method='HEAD'),timeout=30) as f:r['expectedBytes']=int(f.headers['Content-Length'])
 total=sum(r['expectedBytes'] for r in records);free=shutil.disk_usage(root).free;assert free>total+1500000000,(free,total)
 write(root/'source-roster-freeze.json',dict(status='frozen-before-count-download',protocolSHA256=sha(Path(__file__).with_name('PROTOCOL.md')),metadataSHA256=sha(root/'series.soft'),runnerSHA256=sha(Path(__file__)),sourceBytes=total,freeBytes=free,records=records))
 def get(r):
  dest=out/r['path'];part=dest.with_suffix('.partial');h=hashlib.sha256();size=0
  with urllib.request.urlopen(r['url'],timeout=60) as f,part.open('xb') as g:
   while b:=f.read(1048576):
    assert shutil.disk_usage(root).free>1000000000
    g.write(b);h.update(b);size+=len(b)
  assert size==r['expectedBytes'];part.rename(dest);result=dict(r,SHA256=h.hexdigest(),bytes=size);write(out/(r['accession']+'.download.json'),result);print(json.dumps(dict(accession=r['accession'],bytes=size)),flush=True);return result
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:completed=list(pool.map(get,records))
 write(root/'downloads.json',dict(status='completed',files=completed,bytes=sum(r['bytes'] for r in completed),sourceRosterSHA256=sha(root/'source-roster-freeze.json')));print('All24 raw sources downloaded',flush=True)
if __name__=='__main__':main()
