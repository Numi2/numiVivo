"""Hard task admission for Arc 2026; no local surrogate is an official score."""
from wetlab import require
import re
ARC_2026={
 'id':'arc-vcc-2026','task':'zero-shot across cellular contexts',
 'officialSource':'https://arcinstitute.org/news/virtual-cell-challenge-2026',
 'input':'unperturbed expression plus CRISPRi target IDs',
 'split':'three validation cell lines and three final-test cell lines; no challenge training set',
 'scoringOwner':'https://github.com/ArcInstitute/cell-eval2',
 'scoringStatus':'pinned six-metric official evaluator executed locally on public human CRISPRi; official challenge artifacts require login',
 'scorerRevision':'0ce25cf4c495e76e60dfca4ddc8f822507aa0e18',
 'localTask':'Replogle 2020 K562; eight targets; measured-control no-change baseline; not unseen cellular-context qualification',
 'officialSubmissionConformance':False,
 'returnedChallengeScore':None,
 'dimensions':['context','target','gene','pathway','biologicalReplicate','dataset','baseline','externalMethod'],
 'result':'UNAVAILABLE',
 'limits':'The Kang regression and known Clu transfer are not evidence of Arc zero-shot performance.'}

def admit_zero_shot(plan):
    require(plan['benchmarkID']=='arc-vcc-2026','Benchmark identity')
    require(plan['inputState']=='unperturbed' and plan['interventionKind']=='CRISPRi','Arc 2026 input/intervention contract')
    require(plan['trainingProvenanceComplete'] is True,'Incomplete training provenance cannot establish zero-shot status')
    train=set(plan['trainingPerturbedContexts']);query=set(plan['queryContexts'])
    require(query and train.isdisjoint(query),'Query context has a known perturbed training response')
    require(set(plan['queryUnits']).isdisjoint(plan['trainingUnits']),'Biological unit leakage')
    require(plan['heldOutObservationsUsed'] is False,'Held-out outcomes were exposed to model selection/training')
    require(all(isinstance(plan[k],str) and re.fullmatch('[0-9a-f]{64}',plan[k]) for k in ('modelSHA256','inputSHA256')) and isinstance(plan['scorerRevision'],str) and re.fullmatch('[0-9a-f]{40}',plan['scorerRevision']),'Pinned model/input SHA256 and scorer commit required')
    require(plan['queryUnits'] and plan['queryTargets'] and all(isinstance(x,str) and x for x in plan['queryUnits']+plan['queryTargets']),'Query biological units and targets required')
    return {'status':'contract-admitted-only','performance':'UNAVAILABLE','zeroShotDemonstrated':False,
        'unseenContext':True,'unseenPerturbation':set(plan['queryTargets']).isdisjoint(plan['trainingTargets']),
        'requiredComparisons':['no-change','training-mean-where-defined','registered-external-methods'],
        'officialSubmission':False}
