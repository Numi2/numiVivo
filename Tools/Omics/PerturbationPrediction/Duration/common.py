"""Bounded experiment bookkeeping shared by preparation, execution and scoring."""
import hashlib,json
from pathlib import Path

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def write(path,value):
 Path(path).write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')

def verify_files(root,files):
 for name,h in files.items():
  p=Path(name);assert not p.is_absolute() and '..' not in p.parts,name
  assert sha(root/p)==h,name

def hours(condition):
 assert condition.startswith('IFNB:') and condition.endswith('h'),condition
 return float(condition[5:-1])
