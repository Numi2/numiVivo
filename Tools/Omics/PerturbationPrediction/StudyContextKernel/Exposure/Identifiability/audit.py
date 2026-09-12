from pathlib import Path
import argparse,json,hashlib
import numpy as np
parser=argparse.ArgumentParser();parser.add_argument('metadata');parser.add_argument('output');args=parser.parse_args()
p=Path(args.metadata);m=json.loads(p.read_text());durations={'Kang':6,'GSE181897':9,'HIRISA':21};assert set(m['studies'])==set(durations)
studies=sorted(m['studies']);rows=[{'donor':d,'study':s,'protocolDurationHours':durations[s]} for s in studies for d in m['studies'][s]]
assert len(rows)==75 and len({r['donor'] for r in rows})==75
s=np.array([[float(r['study']==v) for v in studies] for r in rows]);t=np.array([r['protocolDurationHours'] for r in rows],dtype=float)
# Study-indicator coefficients reproduce every duration exactly. This algebraic
# check supplements numerical rank and does not depend on a rank tolerance.
coeff=np.array([durations[v] for v in studies],dtype=float);residual=s@coeff-t;assert np.array_equal(residual,np.zeros_like(t))
aug=np.column_stack([s,t]);rank=int(np.linalg.matrix_rank(s));aug_rank=int(np.linalg.matrix_rank(aug));assert rank==aug_rank==3
out={'scope':'Design identifiability audit using actual retained donor membership and provisional study-level protocol durations; not donor-level exposure admission','rows':rows,'columns':studies+['durationHours'],'studyIndicatorRank':rank,'augmentedRank':aug_rank,'columnCount':4,'nullDirection':[-float(durations[v]) for v in studies]+[1.0],'exactDurationReconstruction':True,'independentDurationEffectIdentifiableWithStudyEffects':False,'metadataSHA256':hashlib.sha256(p.read_bytes()).hexdigest(),'limitations':['Duration mapping is provisional at study level; unresolved library linkage is unchanged.','No dose units were converted and no treated values were read.','Removing study effects would impose a structural assumption, not resolve confounding through data.']}
with Path(args.output).open('x') as f:json.dump(out,f,indent=2)
print('PASS: exact duration/study dependence; rank remains 3 after adding duration to three study indicators')
