"""Admit exactly the qualified source version; never silently follow a replacement."""
from pathlib import Path
import os, json, urllib.request
root=Path(os.environ['NUMIVIVO_PARSE_STUDY']);(root/'sources').mkdir(parents=True,exist_ok=True)
url='https://parse-wget.s3.us-west-2.amazonaws.com/10m/Parse_10M_PBMC_cytokines.h5ad'
with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=60) as response:
 headers=dict(response.headers);assert headers['ETag']=='"f0644aed6d5db8e2fdf394523f3e08f2-6780"' and int(headers['Content-Length'])==227497986816
 record={'url':url,'headers':headers}
path=root/'sources/asset-head.json'
if path.exists():
 old=json.loads(path.read_text());assert old['url']==url and old['headers']['ETag']==headers['ETag'] and old['headers']['Content-Length']==headers['Content-Length']
else:path.write_text(json.dumps(record,indent=2)+'\n')
print('Pinned source admitted; CC BY-NC 4.0, Parse Biosciences')
