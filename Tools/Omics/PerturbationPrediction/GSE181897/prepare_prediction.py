"""Admit primary author condition identities and freeze complete prediction inputs.

Run in the existing qualified GSE181897 study using the scientific Python env.
No model is fitted and no prediction is scored here. H5AD query rows represent
source-qualified donor aggregates, not synthetic individual cells.
"""
from pathlib import Path
import gzip, hashlib, json, shutil, time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

root = Path(__file__).resolve().parent
history = root / 'primary-history'
out = root / 'prediction-inputs'
out.mkdir(exist_ok=False)
prior = Path('/Users/home/numivivo-cross-study-ifnb-20260911/inputs')

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        while b := f.read(1048576): h.update(b)
    return h.hexdigest()

def write(p, value):
    p.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')

commit = '5e1f3e2a3ec9ddbf4d594e7b0c0a2f5f2a5fe6f0'
tree = json.loads((history / 'tree.json').read_bytes())
assert not tree['truncated']
assert tree['sha'] == 'e23ff74d73ee97c68b8cdca444e7500e30ecfcff'
sources = {}
notebooks = {}
for upstream, local in [
    ('general.ipynb', 'general.ipynb'),
    ('prod/ifn_compare/ifn.ipynb', 'ifn.ipynb'),
    ('prod/deseq/fix_read_counts.ipynb', 'fix_read_counts.ipynb'),
    ('old/production.run/aggr.combat/diff.ex/fix.read.counts.ipynb', 'old-fix.read.counts.ipynb'),
]:
    p = history / local
    raw = p.read_bytes()
    blob = next(x['sha'] for x in tree['tree'] if x['path'] == upstream)
    assert hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() == blob
    sources[upstream] = dict(gitBlob=blob, SHA256=sha(p), bytes=len(raw),
        URL=f'https://raw.githubusercontent.com/yelabucsf/clue/{commit}/{upstream}')
    notebooks[local] = [''.join(c.get('source', [])) for c in json.loads(raw)['cells']]

# These are source-code identities, not expression-based inference or execution
# of the author's notebook. Cell numbers are zero-based notebook JSON positions.
assert "comp = ('B', 'cM'), ('G', 'cM')" in notebooks['ifn.ipynb'][12]
assert "set_xlabel('Log2FC, IFN-β')" in notebooks['ifn.ipynb'][16]
assert "set_ylabel('Log2FC, IFN-ɣ')" in notebooks['ifn.ipynb'][16]
assert "self.cols = (l0s, l1s, p0s, p1s)" in notebooks['ifn.ipynb'][11]
assert "['A', 'B', 'G' ,'R', 'P']" in notebooks['general.ipynb'][18]
for local, cell in [('fix_read_counts.ipynb', 13), ('old-fix.read.counts.ipynb', 10)]:
    assert 'conditions and controls' in notebooks[local][cell]
    assert "(counts['COND'] == 'C')" in notebooks[local][cell]
    assert "['A','B','G','P','R']" in notebooks[local][cell]
roles = dict(status='primary-author-analysis-code-resolved', repository='https://github.com/yelabucsf/clue',
    commit=commit, tree=tree['sha'], conditionMeaningQualified=True,
    mapping={'B': 'IFNB', 'C': 'control'}, sources=sources,
    evidenceCells={'ifn.ipynb': [11, 12, 16], 'general.ipynb': [18],
                   'fix_read_counts.ipynb': [12, 13], 'old-fix.read.counts.ipynb': [9, 10]},
    scope='Historical primary author analysis code identifies the condition roles. No new exposure measurement, participant-overlap check, or cell-label validation.',
    outcomeExpressionUsedForRoleAssignment=False, originalAdmissionPreserved=True)
write(out / 'primary-condition-roles.json', roles)

src = root / 'handoff'
original = json.loads((src / 'input-freeze.json').read_bytes())
for name, digest in original['files'].items(): assert sha(src / name) == digest, name
native = root / 'native-handoff-qualified'
execution = json.loads((native / 'execution.json').read_bytes())
assert sha(native / 'aggregate/report.json') == execution['reportSHA256']
groups = json.loads((src / 'reference-groups.json').read_bytes())
counts = sparse.load_npz(src / 'reference-pseudobulk.npz')
v = json.loads(gzip.decompress((root / 'var-metadata.json.gz').read_bytes()))
symbols = [s for s, genome in zip(v['_index'], v['genome']) if genome == 'GRCh38']
assert counts.shape == (379, 20303) and len(symbols) == len(set(symbols)) == 20303
assert np.array_equal(np.asarray(counts.sum(axis=1)).ravel(), [g['libraryCounts'] for g in groups])
panel = json.loads((src / 'candidate-panel.json').read_bytes())
assert len(panel) == len(set(panel)) == 11800
assert all(x.removeprefix('symbol|') in symbols for x in panel)
lookup = {(g['expID'], g['conditionCode']): i for i, g in enumerate(groups)}
donors = sorted({g['expID'] for g in groups}, key=int)
eligible, excluded = [], []
for donor in donors:
    missing = [c for c in ['C', 'B'] if (donor, c) not in lookup]
    if missing:
        excluded.append(dict(donor=donor, missingSourceCodes=missing)); continue
    for c in ['C', 'B']:
        g = groups[lookup[donor, c]]
        assert g['cells'] > 0 and g['libraryCounts'] > 0
    eligible.append(donor)
assert len(eligible) == 62 and {g['donor'] for g in excluded} == {'5', '23'}
write(out / 'cohort.json', dict(eligibleDonors=eligible, excludedDonors=excluded,
    sourceGroups=[{k: val for k, val in groups[lookup[d, c]].items() if k != 'sourceCellIndices'}
                  for d in eligible for c in ['C', 'B']],
    responseFeatures=11800, originalRNAFeatures=20303,
    sourceCells=136142, sourceBLineageCells=15272, treatedLabelsConditionEndpoint=True))
write(out / 'panel.json', panel)
shutil.copyfile(src / 'unsupported-prior-panel.json', out / 'unsupported-prior-panel.json')
cross_freeze = json.loads((prior / 'input-freeze.json').read_bytes())
origins = {'Kang': 'cross-Kang-to-HIRISA-2616BW', 'HIRISA': 'cross-HIRISA-to-Kang-patient_101'}
for origin, name in origins.items():
    dest = out / origin; dest.mkdir()
    for filename in ['training.h5ad', 'training.json']:
        assert sha(prior / name / filename) == cross_freeze['files'][name + '/' + filename]
    shutil.copyfile(prior / name / 'training.h5ad', dest / 'training.h5ad')
    plan = json.loads((prior / name / 'training.json').read_bytes())
    plan['responseFeatureIDs'] = panel
    plan['donorResponseIntervalCoverage'] = 0.95
    plan['provenance'] = 'Frozen GSE181897 protocol: complete original training donors, full source RNA denominators, exact 11800-symbol panel; independent query collection, participant overlap unverified.'
    write(dest / 'training.json', plan)
    a = ad.read_h5ad(dest / 'training.h5ad')
    assert all(x in a.var_names for x in panel)
    assert a.n_obs == (16 if origin == 'Kang' else 10)

query_dir = out / 'queries'; query_dir.mkdir()
chunks = []
# Four donors per native publication bounds temporary disk, not biological scope.
for start in range(0, len(eligible), 4):
    selected = eligible[start:start+4]
    key = f'{start // 4:02d}'
    sample_ids = ['GSE181897-control-' + x for x in selected]
    matrix = counts[[lookup[d, 'C'] for d in selected]].copy()
    obs = pd.DataFrame({'sample': sample_ids, 'comparison': ['Bcell-context-transfer']*len(selected)}, index=sample_ids)
    var = pd.DataFrame(index=['symbol|' + x for x in symbols])
    ad.AnnData(matrix, obs=obs, var=var).write_h5ad(query_dir / (key + '.h5ad'), compression='gzip')
    mapping = dict(schemaVersion=1, id='GSE181897-query-' + key, evidence='measured', countUnit='umiCount',
        sourceDescription='Complete source-qualified B-lineage donor CONTROL aggregates. Rows are donor aggregates, not individual cells. Original RNA library denominators retained.',
        matrixPath='X', sampleColumn='sample', groupColumn='comparison',
        samples=[dict(id=s, donorID='GSE181897:exp_id:'+d, biologicalReplicateID='GSE181897:exp_id:'+d,
            condition='control', batchID='multiple-source-technical-pools', organism='NCBITaxon:9606')
            for d, s in zip(selected, sample_ids)])
    write(query_dir / (key + '.json'), dict(schemaVersion=1, mapping=mapping,
        featureNamespace='explicit-source-symbol-correspondence', perturbationID='IFNB-combined-study-context-transfer'))
    chunks.append(dict(id=key, donors=selected))
write(out / 'chunks.json', chunks)
write(out / 'input-freeze.json', dict(schemaVersion=1, createdUnix=time.time(), fitStarted=False, scoringStarted=False,
    primaryConditionMeaningQualified=True, protocolSHA256=sha(root / 'PROTOCOL.md'),
    scriptSHA256=sha(__file__), originalHandoffFreezeSHA256=sha(src / 'input-freeze.json'),
    nativeAggregationExecutionSHA256=sha(native / 'execution.json'), sourceCountValidationSHA256=sha(root / 'source-count-validation.json'),
    originalTrainingFreezeSHA256=sha(prior / 'input-freeze.json'),
    files={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file()}))
print(json.dumps(dict(status='frozen-before-fit-and-score', donors=len(eligible), trainingOrigins=list(origins),
    panel=len(panel), chunks=len(chunks), inputFreezeSHA256=sha(out / 'input-freeze.json'))))
