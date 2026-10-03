"""Additive registry: legacy v1 adapter files stay byte-identical for replay."""
import adapters
from adapters import WetLabExperimentAdapter
import wetlab as rna
from pathlib import Path

def adapter_for_family(family):
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
