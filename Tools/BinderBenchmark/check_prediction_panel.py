#!/usr/bin/env python3
"""Qualify native outcome-blind surface analysis against Biopython on a pinned panel.

This is numerical qualification, not a trained ranker or improved binder selection.
Biopython is a test reference and sequence mapper, not a production dependency.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import platform
import subprocess
import time

import Bio
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
import numpy as np


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def write_json(path: Path, value: object) -> None:
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def fasta(text: str) -> dict[str, str]:
    result, key = {}, None
    for line in text.splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            if key in result:
                raise ValueError('duplicate construct reference')
            result[key] = ''
        elif line.strip():
            if key is None:
                raise ValueError('sequence before FASTA identity')
            result[key] += line.strip()
    return result


def surface_plan(points: int) -> dict:
    return dict(schemaVersion=1, radiusProfile='explicit-elemental-CNOS-v1',
                radiiNM={'C': .17, 'N': .155, 'O': .152, 'S': .18}, probeRadiusNM=.14,
                pointsPerAtom=points, overlapThresholdNM=.04, maximumNeighborPairs=2_000_000,
                maximumPointTests=100_000_000, maximumNeighborSearchTests=10_000_000)


def prepare(model: dict, declared: dict, references: dict, contents: str, points: int) -> tuple:
    record = model['published']
    structure = MMCIFParser(QUIET=True).get_structure('reference', io.StringIO(contents))
    models = list(structure)
    if len(models) != 1:
        raise ValueError('qualification panel expects one coordinate model')
    chains = {c.id: seq1(''.join(r.resname for r in c)) for c in models[0]}
    construct = record['target_construct_id']
    target_sequence = references[construct + '|protein_chain_1']
    binder = [c for c, s in chains.items() if s == declared['sequence']]
    target = [c for c, s in chains.items() if s == target_sequence]
    if len(chains) != 2 or len(binder) != 1 or len(target) != 1 or binder == target:
        raise ValueError('ambiguous or unsupported chain assignment; do not guess')
    if any(a.element not in ('C', 'N', 'O', 'S') for a in structure.get_atoms()):
        raise ValueError('qualification panel must contain canonical CNOS heavy atoms')
    request = dict(schemaVersion=1,
        identity=dict(candidateID=model['candidateID'], predictor=record['predictor'],
                      sampleID='published-seed-' + (record['seed'] or 'unspecified'), targetConstructID=construct),
        binderSequence=declared['sequence'],
        source=dict(candidateID=model['candidateID'], target=model['target'],
                    sourceLabel='Anthropic pinned published model: ' + record['resolvedPath'], format='mmcif',
                    contents=contents, sha256=digest(contents.encode()), targetChainSequences={target[0]: target_sequence},
                    interfacePlan=dict(schemaVersion=1, binderChains=binder, targetChains=target,
                        conformerID='model-1', contactDistanceNM=.45, shortDistanceNM=.2,
                        maximumPairEvaluations=25_000_000)), surfacePlan=surface_plan(points))
    return request, structure


def compare(report: dict, structure: object, request: dict) -> dict:
    plan = request['surfacePlan']
    sr = ShrakeRupley(probe_radius=plan['probeRadiusNM'] * 10, n_points=plan['pointsPerAtom'],
                     radii_dict={k: v * 10 for k, v in plan['radiiNM'].items()})
    atoms = list(structure.get_atoms())
    sr.compute(structure)
    bound = [a.sasa / 100 for a in atoms]
    for chain in structure[0]:
        sr.compute(chain)
    free = [a.sasa / 100 for a in atoms]
    native = report['surface']['atoms']
    if [a['atomIndex'] for a in native] != list(range(len(atoms))):
        raise ValueError('reference/native atom population differs')
    maximum_error = max(max(abs(row['complexAreaNM2'] - b), abs(row['isolatedAreaNM2'] - f))
                        for row, b, f in zip(native, bound, free))
    # Reference sphere directions are float32; allow at most one quadrature
    # point of area per atom, plus a tight total-area tolerance.
    point_area = max(4 * np.pi * (r + plan['probeRadiusNM'])**2 / plan['pointsPerAtom']
                     for r in plan['radiiNM'].values())
    if maximum_error > point_area + 1e-10:
        raise ValueError(f'native/reference atomic SASA mismatch: {maximum_error}')
    native_buried = report['surface']['binderBuriedAreaNM2'] + report['surface']['targetBuriedAreaNM2']
    ref_buried = sum(free) - sum(bound)
    if abs(native_buried - ref_buried) > max(.01, 1e-3 * max(native_buried, ref_buried)):
        raise ValueError('native/reference total buried-area mismatch')
    left_ids = report['interface']['binderHeavyAtomIndices']; right_ids = report['interface']['targetHeavyAtomIndices']
    binder_chain = request['source']['interfacePlan']['binderChains'][0]
    target_chain = request['source']['interfacePlan']['targetChains'][0]
    if left_ids != [i for i, a in enumerate(atoms) if a.parent.parent.id == binder_chain] or right_ids != [i for i, a in enumerate(atoms) if a.parent.parent.id == target_chain]:
        raise ValueError('native/reference partner atom mapping differs')
    left = np.array([atoms[i].coord for i in left_ids], dtype=float) / 10
    right = np.array([atoms[i].coord for i in right_ids], dtype=float) / 10
    lr = np.array([plan['radiiNM'][atoms[i].element] for i in left_ids])
    rr = np.array([plan['radiiNM'][atoms[i].element] for i in right_ids])
    count, maximum = 0, 0.0
    for start in range(0, len(left), 128):
        distances = np.linalg.norm(left[start:start+128, None, :] - right[None, :, :], axis=2)
        overlaps = lr[start:start+128, None] + rr[None, :] - distances
        count += int((overlaps > plan['overlapThresholdNM']).sum())
        maximum = max(maximum, float(overlaps.max()))
    if count != report['surface']['crossPartnerOverlapPairs'] or abs(maximum - report['surface']['maximumCrossPartnerOverlapNM']) > 1e-5:
        raise ValueError('native/reference overlap mismatch')
    return dict(atomCount=len(atoms), maximumAtomAreaErrorNM2=maximum_error,
                nativeBuriedAreaSumNM2=native_buried, referenceBuriedAreaSumNM2=ref_buried,
                maximumAllowedAtomAreaErrorNM2=point_area + 1e-10, overlapPairs=count)


def run(panel: Path, binary: Path, output: Path, points: int) -> dict:
    if output.exists():
        raise FileExistsError(output)
    source = json.loads((panel / 'SOURCE.json').read_text())
    for name, checksum in source['files'].items():
        p = Path(name)
        if p.is_absolute() or '..' in p.parts or any((panel / Path(*p.parts[:i+1])).is_symlink() for i in range(len(p.parts))):
            raise ValueError('unsafe source manifest path')
        if digest((panel / p).read_bytes()) != checksum:
            raise ValueError('source manifest checksum mismatch: ' + name)
    selected = json.loads((panel / 'panel.json').read_text())
    declared = {r['uuid']: r for r in selected['requested']}
    if len(declared) != len(selected['requested']) or any(not isinstance(k, str) or not k for k in declared):
        raise ValueError('duplicate or empty candidate identity')
    if len({r['candidateID'] for r in selected['models']}) != len(selected['models']):
        raise ValueError('duplicate candidate model in source panel')
    refs = fasta((panel / 'target_constructs.fasta').read_text())
    output.mkdir(); (output / 'inputs').mkdir(); (output / 'bundles').mkdir(); (output / 'logs').mkdir()
    protocol = dict(schemaVersion=1, points=points, surfacePlan=surface_plan(points),
        selection=selected['selection'], sourceManifestSHA256=digest((panel/'SOURCE.json').read_bytes()),
        requested=[r['uuid'] for r in selected['requested']], binarySHA256=digest(binary.read_bytes()),
        scope='source-bound numerical qualification only; no labels, training or selection improvement',
        checkerSHA256=digest(Path(__file__).read_bytes()))
    write_json(output / 'protocol.json', protocol)
    completed, failures = [], []
    for model in selected['models']:
        cid = model['candidateID']
        safe_name = digest(cid.encode())
        try:
            if model['sourcePath'] not in source['files'] or model['target'] != declared[cid]['target']:
                raise ValueError('unverified model path or target mismatch')
            raw = (panel / model['sourcePath']).read_bytes()
            if digest(raw) != model['published']['sha256']:
                raise ValueError('model/source mismatch')
            request, structure = prepare(model, declared[cid], refs, raw.decode('utf-8'), points)
            input_file = output / 'inputs' / (safe_name + '.json'); bundle = output / 'bundles' / safe_name
            write_json(input_file, request)
            start = time.perf_counter()
            result = subprocess.run([str(binary), 'binder-analyze-prediction', str(input_file), str(bundle)],
                                    capture_output=True, text=True, timeout=180)
            (output/'logs'/(safe_name+'.stdout')).write_text(result.stdout)
            (output/'logs'/(safe_name+'.stderr')).write_text(result.stderr)
            result.check_returncode()
            native_seconds = time.perf_counter() - start
            report = json.loads((bundle / 'analysis.json').read_text())
            checked = compare(report, structure, request)
            subprocess.run([str(binary), 'binder-verify', str(bundle)], check=True, capture_output=True, timeout=180)
            completed.append(dict(candidateID=cid, target=model['target'], inputSHA256=report['inputSHA256'],
                sourceSHA256=report['sourceSHA256'], bundleName=safe_name, nativePublishSeconds=native_seconds, **checked))
            print(cid, 'PASS', flush=True)
        except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
            failures.append(dict(candidateID=cid, reason=str(error)))
            print(cid, 'FAILED', str(error), flush=True)
    summary = dict(schemaVersion=1, scope=protocol['scope'], requestedCount=len(declared),
        availableCount=len(selected['models']), successCount=len(completed), failures=failures,
        unavailable=selected['unavailable'], completed=completed, python=platform.python_version(), biopython=Bio.__version__,
        originalRevision=source['revision'], binarySHA256=protocol['binarySHA256'])
    write_json(output / 'summary.json', summary)
    if failures or len(completed) != len(selected['models']):
        raise ValueError('not all acquired predictions passed; see retained failure report')
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('panel', type=Path); p.add_argument('binary', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--points', type=int, default=960)
    args = p.parse_args()
    if not 32 <= args.points <= 4096:
        p.error('points must be in 32...4096')
    run(args.panel.resolve(), args.binary.resolve(), args.output.resolve(), args.points)


if __name__ == '__main__':
    main()
