"""Additive registry: legacy v1 adapter files stay byte-identical for replay."""
import adapters
from adapters import WetLabExperimentAdapter
import wetlab as rna
from pathlib import Path

def adapter_for_family(family):
    if family=='measured-cellular-cohort':
        from measured_cellular_cohort import MeasuredCellularCohort
        return MeasuredCellularCohort()
    if family=='learned-cellular-response':
        from cellular_response_adapter import CellularResponseAdapter
        return CellularResponseAdapter()
    if family=='learned-spatial-response':
        from learned_spatial import LearnedSpatialResponseAdapter
        return LearnedSpatialResponseAdapter()
    if family=='molecular-perturbation':
        from molecular_adapter import MolecularPerturbationAdapter
        return MolecularPerturbationAdapter()
    return adapters.adapter_for_family(family)

def adapter_for_config(config):return adapter_for_family(rna.read(config)['family'])
def adapter_for_run(run):
    reg=rna.read(Path(run)/'registration.json');adapter=adapter_for_family(reg.get('family',reg['assay']['family']))
    rna.require(reg.get('adapterID',adapter.adapter_id)==adapter.adapter_id,'Adapter version mismatch')
    return adapter

# New shared operations wrap historical owners without changing their sealed code.
def capabilities(config):
    adapter=adapter_for_config(config);a=rna.read(config)
    if hasattr(adapter,'capabilities'):return adapter.capabilities(config)
    cat=adapter.catalog(config)
    return {'presentation':'tissue' if a['family'] in ('learned-spatial-response','spatial-tissue') else 'legacy',
      'geometry':'optional' if 'spatial' in a['family'] else None,
      'specimens':cat.get('specimens',[]),'populations':cat.get('populations',[]),
      'conditions':a.get('conditions',[{'id':'measured-endpoint','title':'Source measured endpoint'}]),
      'targets':cat.get('targets',[]),'features':cat.get('features',[]),
      'prediction':a['family']=='learned-spatial-response','observationComparison':a['family']=='learned-spatial-response',
      'readouts':['gene'] if a['family']=='learned-spatial-response' else [],'preservationSupported':False}

def catalog(config):
    c=adapter_for_config(config).catalog(config);cap=capabilities(config)
    return {**c,'presentation':cap['presentation'],'geometry':cap['geometry'],'conditions':cap['conditions'],'capabilities':cap}

def validate_selection(config,selection):
    adapter=adapter_for_config(config)
    if hasattr(adapter,'validate_selection'):return adapter.validate_selection(config,selection)
    cap=capabilities(config);rna.require(cap['prediction'],'Assay does not support shared prediction')
    rna.require(selection.get('condition') in {c['id'] for c in cap['conditions']},'Unsupported condition: no registered matched-control/model support')
    rna.require(set(selection)<= {'specimen','population','condition','targets','objective','preregistration'},'Unknown selection field')
    rna.require(not selection['objective']['preserveGenes'] and selection['objective']['penalty']==0,'Spatial preservation unsupported')
    # Historical compiler accepts only its historical fields; validate the complete request first.
    adapter.compile(config,{k:selection[k] for k in ('specimen','population','targets')})
    coverage=inspect_objective(config,selection,selection['objective']);rna.require(coverage['canExecute'],'Objective has unsupported features')
    support=population_support(config,selection,selection['objective']['genes'])
    rna.require(all(any(c['target']==t and c['role']=='direct' and c['canExecute'] for c in support['candidates']) for t in selection['targets']),'Insufficient selected-population reference/outcome support')
    return selection

def inspect_objective(config,selection,objective):
    adapter=adapter_for_config(config)
    if hasattr(adapter,'inspect_objective'):return adapter.inspect_objective(config,selection,objective)
    from objective_inspection import inspect
    return inspect(config,selection,objective['genes'])

def population_support(config,selection,genes):
    adapter=adapter_for_config(config)
    if hasattr(adapter,'support'):return adapter.support(config,selection,genes)
    from investigation import population_support as legacy
    return legacy(config,selection,genes)

def predict(config,runtime,workspace,selection):
    selection=validate_selection(config,selection);adapter=adapter_for_config(config)
    return adapter.predict(config,runtime,workspace,selection if hasattr(adapter,'validate_selection') else {k:selection[k] for k in ('specimen','population','targets')})

def readout(config,record,query,selection):
    adapter=adapter_for_config(config)
    if hasattr(adapter,'readout'):return adapter.readout(config,record,query,selection)
    rna.require(query.get('kind','gene')=='gene','Unsupported readout')
    return adapter.feature(config,selection['specimen'],query['gene'],record)

def evaluate_objective(config,record,objective):
    adapter=adapter_for_config(config)
    if hasattr(adapter,'evaluate_objective'):return adapter.evaluate_objective(config,record,objective)
    from objective_inspection import evaluate
    return evaluate(config,record,objective['genes'])

def scoring_owner(config):
    if rna.read(config)['family']=='learned-cellular-response':
        import cellular_response_adapter
        return Path(cellular_response_adapter.__file__)
    import objective_inspection
    return Path(objective_inspection.__file__)
