#!/usr/bin/env python3
"""Bind the existing frozen HIRISA folds to verified native source group indices."""
import argparse
import hashlib
import json
import time
from pathlib import Path
import ijson

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--repo',type=Path,required=True);a=p.parse_args();root=a.root;repo=a.repo
    protocol=repo/'Tools/Omics/Benchmarks/HIRISA'; frozen=root/'prediction-folds.json'
    assert sha(frozen)=='2d87d76ef9489798707eed33938bf4cd062e69f8b37745a7568432bbc15562a7'
    folds=json.loads(frozen.read_text()); audit=json.loads((root/'native-independent-verification.json').read_text())
    assert audit['status']=='passed' and audit['nativeRawIngestionLinked'] and audit['sourceLibrariesVerified']==131
    assert sha(protocol/'PREDICTION_EXECUTION.md')==folds['predictionProtocolSHA256']
    assert sha(protocol/'PROTOCOL.md')==folds['protocolSHA256']
    assert audit['sourceSHA256']==folds['sourceH5ADSHA256']
    report=root/'native-full/report.json';assert sha(report)==audit['reportSHA256']
    groups={}
    with report.open('rb') as f:
        for i,g in enumerate(ijson.items(f,'pseudobulk.groups.item')):
            assert len(g['sampleIDs'])==1
            accession=g['sampleIDs'][0];assert accession not in groups
            groups[accession]=(i,{k:g[k] for k in ('donorID','condition','batchIDs','biologicalReplicateID')})
    assert len(groups)==131
    native=[]; scoring=[]
    for fold in folds['folds']:
        training=[]
        for pair in fold['trainingPairs']:
            assert pair['donor']!=fold['heldOutDonor']
            for key,condition in [('control','none'),('treated',fold['treatment'])]:
                row,g=groups[pair[key]]
                assert g['donorID']==pair['donor'] and g['condition']==condition and g['batchIDs']==[pair['batch']+'-'+pair['pool']], (pair,g)
                training.append(row)
        query,q=groups[fold['queryControlAccession']];target,t=groups[fold['scoringTargetAccession']]
        assert q['donorID']==t['donorID']==fold['heldOutDonor'] and q['condition']=='none' and t['condition']==fold['treatment']
        assert q['batchIDs']==t['batchIDs'] and query not in training and target not in training
        native.append(dict(id=fold['id'],perturbationID='hirisa-'+fold['treatment'],controlCondition='none',
            treatmentCondition=fold['treatment'],cellGroup=fold['population'],trainingGroupIndices=training,queryGroupIndices=[query]))
        scoring.append(dict(id=fold['id'],queryControlAccession=fold['queryControlAccession'],scoringTargetAccession=fold['scoringTargetAccession']))
    assert len(native)==folds['foldCount']==79
    destination=root/'prediction-batch-plan.json';assert not destination.exists()
    plan=dict(schemaVersion=1,sourceReport={'bytes':list(bytes.fromhex(audit['reportSHA256']))},featureNamespace='ensembl-gene-id',
        provenance='HIRISA GSE306664 frozen 79 donor-held-out folds; original enriched populations; PREDICTION_EXECUTION.md and PREDICTION_TRANSPORT.md; no treated query input.',folds=native)
    destination.write_text(json.dumps(plan,sort_keys=True,indent=2)+'\n')
    (root/'prediction-scoring-targets.json').write_text(json.dumps(scoring,sort_keys=True,indent=2)+'\n')
    paths=['Sources/NumiVivoKit/Omics/VivoPerturbation.swift','Sources/NumiVivoKit/Omics/VivoPerturbationIO.swift',
           'Sources/NumiVivoKit/Omics/VivoPerturbationAggregateBatch.swift','Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift']
    receipt=dict(schemaVersion=1,createdAtUnix=time.time(),foldCount=len(native),frozenFoldsSHA256=sha(frozen),
        planSHA256=sha(destination),scoringTargetsSHA256=sha(root/'prediction-scoring-targets.json'),
        supplementSHA256=sha(protocol/'PREDICTION_TRANSPORT.md'),preparationSHA256=sha(Path(__file__)),
        sourceReportSHA256=audit['reportSHA256'],nativeVerificationSHA256=sha(root/'native-independent-verification.json'),
        originalOwnerSourceSHA256=folds['ownerSourceSHA256'],currentOwnerSourceSHA256={name:sha(repo/name) for name in paths},
        syntheticProductChecksSHA256=sha(root/'prediction-product-check/checks.json'),responseModelFittingStarted=False)
    assert receipt['currentOwnerSourceSHA256'][paths[1]]==receipt['originalOwnerSourceSHA256'][paths[1]]
    (root/'prediction-transport-receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'foldCount':len(native),'planSHA256':receipt['planSHA256']}))

if __name__=='__main__':main()
