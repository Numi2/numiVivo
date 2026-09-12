from pathlib import Path
import json,hashlib
s=Path('/Users/n/numivivo-parse-bcell-counts-20260912')
out=Path('/Users/n/numivivo-parse-endpoint-audit-20260912');out.mkdir(exist_ok=False)
rows=[]
for i in range(1,13):
 p=s/'execution/ingest'/('Donor'+str(i))/'published.json'; receipt=json.loads(p.read_text()); b=s/receipt['bundle']
 report=b/'report.json'; a=json.loads(report.read_text())
 assert set(a)=={'canonicalNonzeros','method','pseudobulk','quality','schemaVersion'}
 assert set(a['pseudobulk'])=={'countUnit','featureIDs','groups','matrix','method'}
 assert len(a['pseudobulk']['featureIDs'])==40352
 assert hashlib.sha256((b/'receipt.json').read_bytes()).hexdigest()==receipt['receiptSHA256']
 rows.append({'donor':receipt['donor'],'cells':receipt['cells'],'reportSHA256':hashlib.sha256(report.read_bytes()).hexdigest(),'reportPath':str(report),'reportKeys':sorted(a),'pseudobulkKeys':sorted(a['pseudobulk']),'hasPerCellLogNormalizedMeans':False})
result={'status':'complete-endpoint-inventory','donors':rows,'cells':sum(r['cells'] for r in rows),'requiredEndpoint':'Per-donor mean of per-cell log1p(1e6 * gene count / full-source cell total)','availableEndpoint':'Integer donor/condition feature count sums plus cell QC','interchangeable':False,'predictionFitted':False,'predictionScored':False,'nextImplementation':'Stream cell counts using full-source cell totals and accumulate per-donor/condition log-normalized feature means; independently compare with sparse reference before model use. Do not substitute log-normalized pseudobulk.'}
assert result['cells']==72446
(out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
print({'status':result['status'],'donors':len(rows),'cells':result['cells'],'interchangeable':False})
