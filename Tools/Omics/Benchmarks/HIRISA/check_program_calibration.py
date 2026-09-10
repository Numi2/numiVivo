#!/usr/bin/env python3
"""Post-result development: fit within-library program variation explicitly.

Reuses previously verified full-source sufficient statistics; original mixed
objective results remain immutable. All outcomes and insensitive controls stay.
"""
import argparse,copy,json,time
from pathlib import Path
import numpy as np
from check_integration_response import sha,validate_design
from check_integration_programs import fit as original_fit,library_metrics,compare


def fit_within(stats, selected, ridge):
    assert ridge>0 and np.all(stats['number'][selected]>0)
    x=stats['x'][selected];y=stats['y'][selected]
    covxx=(stats['xx'][selected]-np.einsum('si,sj->sij',x,x)).mean(axis=0)
    covxy=(stats['xy'][selected]-np.einsum('si,sj->sij',x,y)).mean(axis=0)
    variance=np.diag(covxx)
    assert variance.min()>=-1e-10*max(1.,float(np.abs(stats['xx'][selected]).max()))
    scale=np.sqrt(np.maximum(variance,0));scale[scale<1e-12]=1
    gram=covxx/np.outer(scale,scale);rhs=covxy/scale[:,None]
    w=np.linalg.solve(gram+ridge*np.eye(len(scale)),rhs)/scale[:,None]
    # This intercept only reports an uncentered diagnostic. Within-library score
    # centering removes it; no held-out target mean enters model fitting.
    mean=y.mean(axis=0);intercept=mean-x.mean(axis=0)@w
    assert np.isfinite(w).all() and np.isfinite(intercept).all()
    return w,intercept,mean


def evaluate(previous, ids, protocol, erase=False):
    stats={k:np.asarray(v) for k,v in previous['moments'].items()}
    if erase:
        stats['xx']=np.einsum('si,sj->sij',stats['x'],stats['x'])
        stats['xy']=np.einsum('si,sj->sij',stats['x'],stats['y'])
    index={v:i for i,v in enumerate(ids)};folds=[]
    for old in previous['folds']:
        fold={k:old[k] for k in ['population','treatment','donor','trainingLibraries','evaluationLibraries']}
        training=[index[v] for v in fold['trainingLibraries']];query=[index[v] for v in fold['evaluationLibraries']]
        if np.any(stats['number'][training+query]==0):
            fold.update(status='unavailable-empty-library',measurements=[dict(programID=p,withinLibraryR2=None) for p in protocol['programs']])
        else:
            w,b,mean=fit_within(stats,training,protocol['ridge'])
            metrics=[library_metrics(stats,k,w,b,mean,protocol['margins']['minimumProgramVariance']) for k in query]
            measurements=[]
            for j,p in enumerate(protocol['programs']):
                values=[m[j]['withinLibraryR2'] for m in metrics]
                measurements.append(dict(programID=p,libraries=[m[j] for m in metrics],withinLibraryR2=float(np.mean(values)) if all(v is not None for v in values) else None))
            fold.update(status='measured',weight=w.tolist(),intercept=b.tolist(),measurements=measurements)
        folds.append(fold)
    result={k:copy.deepcopy(previous[k]) for k in ['cells','matchedOriginalCells','matchedProgramTargetCells','cellsOutsideMatchedFolds','originalLibraryCells','programTargetLibraryCells','unavailableEmptyCellTargets']}
    result.update(folds=folds,moments={k:v.tolist() for k,v in stats.items()},libraryMeanOnly=erase)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--previous',type=Path,required=True);p.add_argument('--protocol',type=Path,required=True);p.add_argument('--design',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();assert not a.out.exists()
    start=time.monotonic();protocol=json.loads(a.protocol.read_text());_,ids=validate_design(json.loads(a.design.read_text()))
    assert sha(a.design)==protocol['designSHA256']
    inputs={name:a.previous/(name+'.json') for name in protocol['representations']}
    assert {n:sha(path) for n,path in inputs.items()}==protocol['previousResultsSHA256']
    previous={n:json.loads(path.read_text()) for n,path in inputs.items()}
    for n,old in previous.items():
        assert old['status']=='measured' and not old['libraryMeanOnly'] and old['cells']==1612594
        assert old['bindings']['protocol']['SHA256']==protocol['previousProtocolSHA256']
        for k in ['source','metadata','design','programReference','programScores']:assert old['bindings'][k]==previous['baseline']['bindings'][k]
        assert old['protocol']['ridge']==protocol['ridge'] and old['protocol']['margins']==protocol['margins'] and old['protocol']['programs']==protocol['programs']
        # Confirm exact moment reuse reconstructs the already-qualified old fit.
        stats={k:np.asarray(v) for k,v in old['moments'].items()};index={v:i for i,v in enumerate(ids)}
        for fold in old['folds']:
            if fold['status']=='measured':
                w,b,_=original_fit(stats,[index[v] for v in fold['trainingLibraries']],protocol['ridge'])
                np.testing.assert_array_equal(w,fold['weight']);np.testing.assert_array_equal(b,fold['intercept'])
    a.out.mkdir(parents=True,exist_ok=False)
    results={n:evaluate(old,ids,protocol) for n,old in previous.items()}
    results['identity']=evaluate(previous['baseline'],ids,protocol)
    results['library-mean-only']=evaluate(previous['baseline'],ids,protocol,True)
    assert results['identity']==results['baseline']
    bindings={k:v for k,v in previous['baseline']['bindings'].items() if k not in ['baseline','scores']}
    for name,result in results.items():
        source_name=name if name in previous else 'baseline';old=previous[source_name]
        result.update(status='measured',protocol=protocol,scope=protocol['scope'],bindings={**bindings,'scores':old['bindings']['scores']},
                      previousResultSHA256=sha(inputs[source_name]),calibrationProtocolSHA256=sha(a.protocol),evaluatorSHA256=sha(Path(__file__)),
                      dependencies={n:sha(Path(__file__).with_name(n)) for n in ['check_integration_programs.py','check_integration_response.py']})
        if name!='baseline':
            result['comparisons']=compare(results['baseline'],result,protocol)
            result['allProgramGradientGatesPassed']=all(all(v['gates'].values()) for v in result['comparisons'])
        if name=='library-mean-only':
            assert all(abs(m['withinLibraryR2'])<1e-10 for f in result['folds'] for m in f['measurements'] if m['withinLibraryR2'] is not None)
        (a.out/(name+'.json')).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    assert {n:sha(path) for n,path in inputs.items()}==protocol['previousResultsSHA256']
    summary={}
    for name,result in results.items():
        if 'comparisons' not in result:continue
        comparisons=result['comparisons'];failures=[c for c in comparisons if not c['gates']['meanPreserved'] or not c['gates']['everyFoldPreserved']]
        summary[name]=dict(sensitive=sum(c['controlSensitive'] for c in comparisons),comparisons=len(comparisons),allGates=result['allProgramGradientGatesPassed'],marginFailures=[{k:c[k] for k in ['population','treatment','programID','controlSensitive','meanWithinLibraryR2Loss','maximumFoldWithinLibraryR2Loss']} for c in failures])
    (a.out/'complete.json').write_text(json.dumps(dict(status='completed-development-evaluation',seconds=time.monotonic()-start,summary=summary,protocolSHA256=sha(a.protocol),results={n:sha(a.out/(n+'.json')) for n in results}),indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
