#!/usr/bin/env python3
"""Run fixed, outcome-blind reference rankings on the pinned three-target panel.

No new inference, structural features, fitting, threshold search or label-policy changes.
This reused panel is development evidence, not independent validation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import tempfile

from check_rankings import FEATURES, ROOT, command, verify_ranking, verify_assessment
from fetch_source import BLOB, REVISION, TARGETS, URL

SOURCE_SHA256 = 'a96dc3a4fcc293dd393a8c56fa367dadcac4db326c30f1c3ea3d8c2b9f72367a'
METHODS = ('boltz2', 'confidence3')


def sha(path: Path) -> str:
    with path.open('rb') as h:
        digest = hashlib.sha256()
        for block in iter(lambda: h.read(1 << 20), b''): digest.update(block)
    return digest.hexdigest()


def run(source: Path, binary: Path, output: Path) -> dict:
    if output.exists() or not output.parent.is_dir():
        raise ValueError('new output directory with existing parent required')
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SOURCE_SHA256:
        raise ValueError('source SHA-256 differs from pinned table')
    if hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() != BLOB:
        raise ValueError('source Git blob identity differs from pin')
    with tempfile.TemporaryDirectory(prefix='.binder-ranking-', dir=output.parent) as temporary:
        root = Path(temporary)/'campaign'; root.mkdir()
        (root/'source.csv').write_bytes(raw)
        (root/'plans').mkdir()
        plans = {}
        for method in METHODS:
            destination = root/f'plans/{method}.json'
            shutil.copyfile(ROOT/f'Tools/BinderBenchmark/plans/ranking-{method}.json', destination)
            plans[method] = sha(destination)
        protocol = {'schemaVersion':1,'sourceSHA256':SOURCE_SHA256,'sourceGitBlobSHA1':BLOB,
            'sourceRevision':REVISION,'sourceURL':URL,'targets':TARGETS,'planSHA256':plans,
            'scope':'fixed reference scoring on already-inspected development panel; no fitting',
            'normalization':'within-target population SD over all complete query candidates, including untested outcomes',
            'selection':'top 10 before outcome join; no replacement of unavailable outcomes; fractional cutoff ties',
            'confidence3':'equal standardized seed-best ipSAE from ef2fast, ef2full, ptxv2; no optional scDockQ terms',
            'notClaimed':['campaign-time inference reproduction','equal GPU inference cost','Numi physical uplift','independent biological validation']}
        # Freeze plans/protocol before imports/evaluation; timestamps alone cannot attest pre-registration.
        (root/'protocol.json').write_text(json.dumps(protocol,sort_keys=True,indent=2)+'\n')
        shutil.copy2(binary,root/'native-cli')
        implementation = sha(binary)
        results = []; reference_query = None; reference_ranks = {}
        for assay in ('adaptyv','twist'):
            config = root/f'import-{assay}.json'
            config.write_text(json.dumps({'schemaVersion':1,'assay':assay,'targets':TARGETS,'features':FEATURES},sort_keys=True)+'\n')
            imported = root/f'input-{assay}'; query = root/f'query-{assay}'
            command(binary,'binder-import',root/'source.csv',config,imported)
            command(binary,'binder-ranking-query',imported,query)
            query_bytes = (query/'query.json').read_bytes()
            if reference_query is not None and query_bytes != reference_query:
                raise AssertionError('laboratory outcomes changed query projection')
            reference_query = query_bytes
            q = json.loads(query_bytes)
            records = json.loads((imported/'imported.json').read_text())['dataset']['records']
            if len(q['candidates']) != 270: raise AssertionError('unexpected candidate population')
            for method in METHODS:
                plan = root/f'plans/{method}.json'; ranked = root/f'ranking-{assay}-{method}'
                assessed = root/f'assessment-{assay}-{method}'
                command(binary,'binder-rank',query,plan,ranked)
                command(binary,'binder-assess-ranking',imported,ranked,assessed)
                for bundle in (query,ranked,assessed):
                    command(binary,'binder-verify',bundle)
                    receipt=json.loads((bundle/'receipt.json').read_text())
                    if receipt['implementationSHA256'] != implementation: raise AssertionError('implementation mismatch')
                    for name, digest in receipt['files'].items():
                        if sha(bundle/name) != digest: raise AssertionError('artifact checksum mismatch')
                ranking_bytes = (ranked/'ranking.json').read_bytes()
                if method in reference_ranks and ranking_bytes != reference_ranks[method]:
                    raise AssertionError('laboratory outcomes changed fixed ranking')
                reference_ranks[method] = ranking_bytes
                r = json.loads(ranking_bytes); a = json.loads((assessed/'assessment.json').read_text())
                verify_ranking(q,json.loads(plan.read_text()),r)
                verify_assessment(records,r,a)
                for t in a['targets']: results.append({'assay':assay,'method':method,**t})
        code = {}
        for folder in ('Sources/NumiVivoKit/Binder','Tools/BinderBenchmark'):
            for path in sorted((ROOT/folder).glob('*')):
                if path.is_file() and path.suffix in ('.swift','.py'): code[str(path.relative_to(ROOT))]=sha(path)
        summary = {'schemaVersion':1,'status':'completed','evidenceStatus':'reused-retrospective-development',
            'sourceSHA256':SOURCE_SHA256,'implementationSHA256':implementation,'protocolSHA256':sha(root/'protocol.json'),
            'independentNumericalChecks':'passed','labIndependentRankingBytes':'identical',
            'sourceHashes':code,'results':results,'biologicalImprovement':'not-established'}
        (root/'summary.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
        manifest = {str(p.relative_to(root)):sha(p) for p in sorted(root.rglob('*')) if p.is_file()}
        (root/'manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
        if output.exists(): raise FileExistsError(output)
        os.rename(root,output)
        return summary


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path); p.add_argument('binary',type=Path); p.add_argument('output',type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.source.resolve(),args.binary.resolve(),args.output.absolute()),sort_keys=True))
