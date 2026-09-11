"""Inspect literal feature correspondence without fitting or looking at outcomes."""
from pathlib import Path
import argparse,os,json,hashlib
p=argparse.ArgumentParser();p.add_argument('--panel',type=Path,required=True);a=p.parse_args();root=Path(os.environ['NUMIVIVO_PARSE_STUDY']);panel=json.loads(a.panel.read_text());features=json.loads((root/'sources/roster.json').read_text())['sourceFeatures']
assert len(panel)==len(set(panel)) and all(x.startswith('symbol|') for x in panel)
namespaced={'symbol|'+x for x in features};missing=[x for x in panel if x not in namespaced]
record=dict(status='identifier-compatibility-only',sourceFeatureCount=len(features),fixedDurationPanelCount=len(panel),present=len(panel)-len(missing),missing=missing,panelSHA256=hashlib.sha256(a.panel.read_bytes()).hexdigest(),symbolCorrespondence='Add the established symbol| namespace wrapper to the unchanged Parse source symbol; exact character equality only, no gene synonym remapping or outcome-driven panel change',predictionFitOrScoring=False)
(root/'feature-panel-compatibility.json').write_text(json.dumps(record,indent=2)+'\n');print({k:v for k,v in record.items() if k!='missing'})
