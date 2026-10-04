#!/usr/bin/env python3
"""Retain the executed design experiment and portable native query replay.

Packaging never changes source-bound records or their historical verdicts.
Full observation re-extraction still requires the original source paths.
"""
import argparse
import shutil
import subprocess
import tarfile
from pathlib import Path
from wetlab import read, write, sha, inventory, require


def package(evidence, output):
    require(not output.exists(), 'Package exists; do not overwrite evidence')
    output.mkdir()
    runtime = output / 'runtime'
    runtime.mkdir()
    shutil.copy2(evidence / 'runtime/numivivo', runtime / 'numivivo')
    shutil.copytree(evidence / 'runtime/mlx-swift_Cmlx.bundle', runtime / 'mlx-swift_Cmlx.bundle')
    for name in ('learning-v3', 'campaign-v1', 'receivers-v3', 'target-sources', 'experiments',
                 'followup/learning-v2', 'followup/campaign', 'followup/experiments', 'browser', 'research',
                 'training-group-diagnostic', 'followup/training-group-diagnostic'):
        shutil.copytree(evidence / name, output / name)
    for name in ('preregistration.json', 'selection-evaluation.json', 'diagnostics-final.json',
                 'training-fit-diagnostics.json', 'lifecycle.json', 'v03-regression.json',
                 'v02-replay.json', 'safari-session.json', 'receiver-qualification.json',
                 'followup/preregistration.json', 'followup/sources-manifest.json',
                 'followup/objective-admission.json', 'followup/evaluation.json',
                 'followup/diagnostics-final.json'):
        shutil.copy2(evidence / name, output / name)
    source = output / 'owner-source'
    source.mkdir()
    for p in Path(__file__).parent.glob('*.py'):
        shutil.copy2(p, source / p.name)
    historical = source / 'initial-campaign-27ec9f27'
    historical.mkdir()
    for name in ('prepare_intervention_design.py', 'train_intervention_design.py'):
        (historical / name).write_bytes(subprocess.check_output(
            ['git', 'show', '27ec9f27:Tools/VirtualWetLab/' + name]))
    require(sha(historical / 'train_intervention_design.py') ==
            read(evidence / 'campaign-v1/registration.json')['ownerSHA256'], 'Initial trainer source changed')
    shutil.copy2(Path(__file__).with_name('V04.md'), output / 'V04.md')
    models = []
    for prefix in ('campaign-v1', 'followup/campaign', 'receivers-v3'):
        variants = ('receiver-pretrained', 'anchor-reference-pretrained', 'receiver-spatial-only') if prefix == 'receivers-v3' else ('target-descriptor', 'target-ID', 'shuffled-target-descriptor')
        for variant in variants:
            folder = output / prefix / variant
            step = read(folder / 'selection.json')['selected']['step']
            receiver = prefix == 'receivers-v3'
            models.append({'name': prefix + '/' + variant,
                           'plan': str((folder / 'plan.json').relative_to(output)),
                           'weights': str((folder / 'training' / f'weights-{step}.safetensors').relative_to(output)),
                           'query': str((folder / ('chip1.safetensors' if receiver else 'reserved.safetensors')).relative_to(output)),
                           'prediction': str((folder / ('development-prediction' if receiver else 'reserved-prediction') / 'prediction.safetensors').relative_to(output))})
    write(output / 'native-models.json', models)
    (output / 'README.txt').write_text(
        'Virtual Wet Lab v0.4 DEVELOPMENT INCREMENT\n'
        'Full requested milestone: NOT COMPLETE. Biological promotion: FAILED / UNAVAILABLE.\n'
        'No useful target-aware selection advantage established. Spatial edge admission and preservation campaign unavailable. Safari qualification blocked.\n\n'
        'Portable same-Mac native replay (Python needs numpy and safetensors):\n'
        '  python owner-source/verify_design_package.py PACKAGE NEW_REPLAY_DIRECTORY\n'
        'All nine selected native candidates must reproduce mean and variance tensors bit-exactly. Cross-hardware bit identity is not claimed.\n\n'
        'Training inputs, all fixed-budget checkpoints, plans, target sources, experiment records and ablations are retained. See V04.md for executed invocations, source admission and failures.\n'
        'Original absolute paths remain in sealed records. Portable native replay uses native-models.json; it does not rewrite those records. Full source re-extraction and cohort rebuild require original GEO/scPerturb files and the preserved v0.3 cohort. Those large raw files are not duplicated in this archive.\n')
    write(output / 'manifest.json', {'format': 'numivivo-design-development-package/v1',
          'sourceCommit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
          'completeV04Milestone': False, 'biologicalPromotion': False, 'files': inventory(output)})
    archive = output.with_suffix('.tar.gz')
    with tarfile.open(archive, 'w:gz', compresslevel=6) as t:
        t.add(output, arcname=output.name)
    print(archive, sha(archive))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    package(a.evidence, a.output)
