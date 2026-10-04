"""Native inference in batch-aligned chunks; bounds disk without changing the model."""
import hashlib,tempfile
from pathlib import Path
import numpy as np
from safetensors.numpy import load_file,save_file
from train_intervention_design import native
from wetlab import read

def predict(binary,plan,data,weights):
    a=load_file(str(data));p=read(plan);size=p['batchSize']*8;pieces=[]
    for start in range(0,len(a['context']),size):
        with tempfile.TemporaryDirectory(prefix='spatial-chunk-') as tmp:
            root=Path(tmp);q={k:(v if k in ('prior','responseScale') else v[start:start+size]) for k,v in a.items() if k in ('context','contextMask','descriptor','target','known','prior','responseScale')};save_file(q,str(root/'query.safetensors'))
            try:native(binary,'predict',plan,root/'query.safetensors',root/'prediction',weights)
            except Exception as e:raise RuntimeError(str(e)+' '+(root/'prediction.log').read_text()[-4000:]) from e
            pieces.append(load_file(str(root/'prediction/prediction.safetensors')))
    result={k:np.concatenate([x[k] for x in pieces]) for k in pieces[0]};digest=hashlib.sha256()
    for k in sorted(result):digest.update(k.encode());digest.update(result[k].tobytes())
    return result,digest.hexdigest()
