from pathlib import Path
import csv,io,gzip,json,hashlib,collections
base=Path('/Users/home/numivivo-parse-bcell-admission-20260912')
hg=Path('/Users/home/numivivo-parse-ifnb-20260911/feature-identity')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((base/'audit.json').read_text());source=json.loads((base/'source-roster.json').read_text())['sourceFeatures']
raw=gzip.decompress((hg/'hgnc_complete_set.txt.gz').read_bytes());identity=json.loads((hg/'source.json').read_text())
assert hashlib.sha256(raw).hexdigest()==identity['SHA256']
rows=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'))
lookup=collections.defaultdict(list);coords=collections.defaultdict(list)
for i,s in enumerate(source):coords[s].append(i)
for r in rows:
 for s in set([r['symbol']]+r['prev_symbol'].split('|'))- {''}:lookup[s].append(r)
dest=collections.defaultdict(list)
for m in a['mapping']:dest[m['sourceIndex']].append(m['modelID'])
records=[]
for m in a['missing']:
 sym=m['modelID'].removeprefix('symbol|');matches=lookup[sym];c=[]
 for r in matches:
  for s in sorted(set([r['symbol']]+r['prev_symbol'].split('|')) & set(coords)):
   for i in coords[s]:c.append(dict(hgncID=r['hgnc_id'],approvedSymbol=r['symbol'],sourceSymbol=s,sourceIndex=i,ensemblGeneID=r['ensembl_gene_id']))
 rec=dict(**m,matchedHGNCIDs=sorted({r['hgnc_id'] for r in matches}),candidates=c)
 if not matches:status='no-HGNC-approved-or-previous-symbol'
 elif len(matches)>1:status='HGNC-identity-ambiguous'
 elif not c:status='no-corresponding-source-symbol'
 elif len(c)>1:status='source-correspondence-ambiguous'
 else:status='candidate';dest[c[0]['sourceIndex']].append(m['modelID'])
 rec['status']=status;records.append(rec)
for r in records:
 if r['status']=='candidate':r['status']='unique-candidate' if len(dest[r['candidates'][0]['sourceIndex']])==1 else 'source-coordinate-conflict'
# Independent direct row scan checks every identity relation.
for r in records:
 sym=r['modelID'].removeprefix('symbol|')
 direct=[v for v in rows if v['symbol']==sym or sym in v['prev_symbol'].split('|')]
 assert {v['hgnc_id'] for v in direct}==set(r['matchedHGNCIDs'])
 pairs={(v['hgnc_id'],i) for v in direct for i,s in enumerate(source) if s==v['symbol'] or s in v['prev_symbol'].split('|')}
 assert pairs=={(c['hgncID'],c['sourceIndex']) for c in r['candidates']}
out=dict(status='complete-context-feature-correspondence-audit',modelFeatures=11800,exactMatches=11600,missing=200,classification=dict(collections.Counter(r['status'] for r in records)),records=records,sources={str(p):sha(p) for p in [base/'audit.json',base/'source-roster.json',hg/'hgnc_complete_set.txt.gz',hg/'source.json']},independentIdentityAndCoordinateCheck=True,modelPanelModified=False,sourceCountsModified=False,predictionScored=False,policy='Frozen retained HGNC approved/previous symbols only; no alias_symbol or expression outcomes. Candidates do not establish sequence/annotation equivalence or authorize remapping.')
output=Path('/Users/home/numivivo-expression-memory-20260911/Tools/Omics/PerturbationPrediction/ParseIFNB/BCellAdmission/context-feature-correspondence.json')
output.write_text(json.dumps(out,indent=2)+'\n')
print(out['classification'])

