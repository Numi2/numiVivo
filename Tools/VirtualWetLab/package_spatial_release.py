#!/usr/bin/env python3
"""Package qualified weights, native inference and sealed query artifacts."""
import argparse,shutil,subprocess,tarfile
from pathlib import Path
from wetlab import read,write,sha
from evaluate_spatial_response import VARIANTS

def package(evidence,output):
 evidence=Path(evidence);output=Path(output);output.mkdir(exist_ok=False);inputs=evidence/'learning-v2';runtime=evidence/'native-build/Build/Products/Release';(output/'runtime').mkdir()
 shutil.copy2(runtime/'numivivo',output/'runtime/numivivo');shutil.copytree(runtime/'mlx-swift_Cmlx.bundle',output/'runtime/mlx-swift_Cmlx.bundle')
 for variant in VARIANTS:
  source=inputs/variant;folder=output/'models'/variant;folder.mkdir(parents=True);selection=read(source/'selection.json')
  for src,dst in [(source/'training'/('weights-'+str(selection['selected']['step'])+'.safetensors'),'weights.safetensors'),(source/'plan.json','plan.json'),(source/'chip1.safetensors','query.safetensors'),(source/'test-prediction/prediction.safetensors','prediction.safetensors'),(source/'selection.json','selection.json'),(source/'prediction-seal.json','original-seal.json')]:shutil.copy2(src,folder/dst)
 for name in ('preregistration.json','features.json','reference-pca.npz','prepared.json'):shutil.copy2(inputs/name,output/name)
 for name in ('representation-admission.json','representation-receipt.json','gene-mapping.json','frozen-pca.npz'):shutil.copy2(inputs/'frozen-Nicheformer'/name,output/name)
 records=output/'qualification';records.mkdir()
 for name in ('cohort-qualification-v2.json','lifecycle-final.json','training-reproduction.json','browser-performance.json','browser-workflow-final.json','v02-replay.json'):shutil.copy2(evidence/name,records/name)
 shutil.copytree(evidence/'evaluation-v1',records/'spatial-evaluation');shutil.copytree(evidence/'arc-local-v1/official-evaluation-repaired',records/'arc-local-official');shutil.copy2(evidence/'arc-local-v1/successful-execution.json',records/'arc-local-execution.json')
 record=Path(read(evidence/'lifecycle-final.json')['run']);shutil.copytree(record,output/'experiment-record')
 (output/'README.txt').write_text('Virtual Wet Lab v0.3 experimental spatial-response candidate\n\nBiological promotion: FAILED. No useful neighborhood advantage established; no-change beats this model. Four isolated test controls are below the required five. Cfap410 is unseen in training but exposed during validation selection; not clean unseen-perturbation qualification.\n\nNative inference on Apple silicon macOS:\n  runtime/numivivo spatial-response predict models/neighborhood/plan.json models/neighborhood/query.safetensors replay models/neighborhood/weights.safetensors\nCompare replay/prediction.safetensors against models/neighborhood/prediction.safetensors as float tensors. Same-host inference and complete retraining were bit-exact; cross-hardware bit identity is not claimed.\n\nEvery model is the existing NumiVivo MLX cell-response network with reference-cell and neighborhood conditioning. Query tensors contain no measured responses. No-change, training-mean, cell-type-matched and upstream State PerturbMean results are in qualification/spatial-evaluation.\n\nFull biological source and cohort regeneration: see Tools/VirtualWetLab/V03.md at the manifest source commit. The original GEO GEF files and full spatial cohort are retained separately; no source download is replaced by generated biology. Full Python experiment-record replay requires those original source paths and the pinned owner implementation. This package supports portable native query replay independently.\n\nArc results here are local official-evaluator scores on Replogle-2020 K562 CRISPRi, NOT VCC2026 submission conformance, a challenge score, or spatial validation.\n')
 revision=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();write(output/'manifest.json',{'format':'numivivo-spatial-candidate-release/v1','sourceCommit':revision,'biologicalPromotion':False,'files':{str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}})
 archive=output.with_suffix('.tar.gz')
 with tarfile.open(archive,'w:gz',compresslevel=6) as t:t.add(output,arcname=output.name)
 print(archive,sha(archive))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--evidence',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();package(a.evidence,a.output)
