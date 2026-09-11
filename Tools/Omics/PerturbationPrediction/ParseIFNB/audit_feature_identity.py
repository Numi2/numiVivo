"""Audit absent fixed-panel symbols using HGNC-approved/previous symbols only.

This emits correspondence candidates, never modifies counts or the model panel.
"""
import argparse,csv,gzip,hashlib,io,json
from pathlib import Path
from collections import defaultdict,Counter


def audit(panel,features,rows):
 assert len(panel)==len(set(panel)) and len(features)==len(set(features))
 assert all(p.startswith('symbol|') and p.count('|')==1 for p in panel)
 feature_set=set(features);symbols={p:p.removeprefix('symbol|') for p in panel};exact={p:s for p,s in symbols.items() if s in feature_set};index=defaultdict(list)
 assert len({r['hgnc_id'] for r in rows})==len(rows)
 for row in rows:
  assert row['status']=='Approved'
  for symbol in {row['symbol']}|set(filter(None,row['prev_symbol'].split('|'))):index[symbol].append(row)
 records=[]
 for name,symbol in symbols.items():
  if name in exact:continue
  matches=index.get(symbol,[]);candidates=[]
  for row in matches:
   accepted={row['symbol']}|set(filter(None,row['prev_symbol'].split('|')))
   for candidate in sorted(accepted & feature_set):
    candidates.append(dict(sourceSymbol=candidate,hgncID=row['hgnc_id'],approvedSymbol=row['symbol'],panelRole='approved' if symbol==row['symbol'] else 'previous',sourceRole='approved' if candidate==row['symbol'] else 'previous',ensemblGeneID=row['ensembl_gene_id'],authority='HGNC symbol/prev_symbol, not alias_symbol'))
  status='no-HGNC-approved-or-previous-symbol' if not matches else ('HGNC-identity-ambiguous' if len(matches)!=1 else ('no-corresponding-source-symbol' if not candidates else ('source-correspondence-ambiguous' if len(candidates)!=1 else 'unique-candidate')))
  records.append(dict(panelFeature=name,status=status,matchedHGNCIDs=[r['hgnc_id'] for r in matches],candidates=candidates))
 claimed=defaultdict(list)
 for name,source in exact.items():claimed[source].append(name)
 for r in records:
  if r['status']=='unique-candidate':claimed[r['candidates'][0]['sourceSymbol']].append(r['panelFeature'])
 for r in records:
  if r['status']=='unique-candidate':
   owners=claimed[r['candidates'][0]['sourceSymbol']]
   if len(owners)!=1:r.update(status='source-coordinate-conflict',conflicts=owners)
 statuses=dict(sorted(Counter(r['status'] for r in records).items()));return dict(status='annotation-only-correspondence-audit',panelCount=len(panel),sourceFeatureCount=len(features),exactMatches=len(exact),exactMissing=len(records),classification= statuses,records=records,modelPanelModified=False,sourceCountsModified=False,predictionFitOrScoring=False,policy='Exact approved or explicitly curated previous symbols, no alias_symbol, case folding, suffix stripping, fuzzy match, count-dependent selection or zero fill. Candidates require one HGNC identity, one source coordinate and no collision with any exact or candidate panel coordinate. Source H5AD supplies no stable gene IDs; an HGNC symbol correspondence alone does not establish sequence/annotation equivalence. The original exact-name gate remains unchanged.')


def main():
 p=argparse.ArgumentParser();p.add_argument('--panel',type=Path,required=True);p.add_argument('--roster',type=Path,required=True);p.add_argument('--hgnc',type=Path,required=True);p.add_argument('--hgnc-sha256',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 raw=gzip.decompress(a.hgnc.read_bytes());assert hashlib.sha256(raw).hexdigest()==a.hgnc_sha256
 rows=list(csv.DictReader(io.StringIO(raw.decode()),delimiter='\t'));panel=json.loads(a.panel.read_text());features=json.loads(a.roster.read_text())['sourceFeatures'];result=audit(panel,features,rows)
 result['sources']={key:dict(path=str(path),SHA256=hashlib.sha256(path.read_bytes()).hexdigest()) for key,path in [('panel',a.panel),('sourceRoster',a.roster),('hgncGzip',a.hgnc)]};result['HGNCDecodedSHA256']=a.hgnc_sha256;result['HGNCApprovedRows']=len(rows)
 assert not a.out.exists();a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
if __name__=='__main__':main()
