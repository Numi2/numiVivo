"""Versioned experiment lifecycle; assay owners define quantities and evidence.

RNA v1 records retain their original verifier unchanged. New adapters must not
reinterpret another assay's observable, units or qualification verdict.
"""
from abc import ABC, abstractmethod
from pathlib import Path
import argparse
import json
import wetlab as rna


class WetLabExperimentAdapter(ABC):
    adapter_id: str
    family: str

    @abstractmethod
    def catalog(self, config): pass
    @abstractmethod
    def predict(self, config, runtime, workspace, selection): pass
    @abstractmethod
    def reveal(self, run, runtime): pass
    @abstractmethod
    def verify(self, run, runtime): pass
    @abstractmethod
    def summary(self, run): pass


class RNAResponseAdapter(WetLabExperimentAdapter):
    adapter_id = 'rna-response/v1'
    family = 'cell-response'

    def catalog(self, config):
        return {**rna.catalog(config), 'adapterID': self.adapter_id}
    def predict(self, config, runtime, workspace, selection):
        return rna.predict(config, runtime['binary'], workspace, selection['donor'], selection['hours'], selection['intervention'])
    def reveal(self, run, runtime):
        return rna.reveal(run, runtime['binary'])
    def verify(self, run, runtime):
        return rna.verify(run, runtime['binary'])
    def summary(self, run):
        return {**rna.summary(run), 'adapterID': self.adapter_id, 'family': self.family}


def adapter_for_family(family):
    if family == 'cell-response': return RNAResponseAdapter()
    if family == 'spatial-tissue':
        from spatial_adapter import SpatialTissueAdapter
        return SpatialTissueAdapter()
    raise ValueError('Unsupported assay family: ' + str(family))


def adapter_for_config(config):
    return adapter_for_family(rna.read(config)['family'])


def adapter_for_run(run):
    reg = rna.read(Path(run) / 'registration.json')
    family = reg.get('family', reg.get('assay', {}).get('family'))
    adapter = adapter_for_family(family)
    rna.require(reg.get('adapterID', adapter.adapter_id) == adapter.adapter_id, 'Record adapter version differs')
    return adapter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['catalog','predict','reveal','verify','inspect'])
    parser.add_argument('path', type=Path)
    parser.add_argument('--binary', type=Path)
    parser.add_argument('--workspace', type=Path)
    parser.add_argument('--selection', type=Path, help='Adapter-specific selection JSON')
    args=parser.parse_args();runtime={'binary':args.binary}
    if args.command in ('catalog','predict'): adapter=adapter_for_config(args.path)
    else: adapter=adapter_for_run(args.path)
    if args.command=='catalog': result=adapter.catalog(args.path)
    elif args.command=='predict':
        rna.require(args.workspace is not None and args.selection is not None,'predict requires workspace and selection')
        result={'run':str(adapter.predict(args.path,runtime,args.workspace,rna.read(args.selection)))}
    elif args.command=='reveal': result=adapter.reveal(args.path,runtime)
    elif args.command=='verify': result=adapter.verify(args.path,runtime)
    else: result=adapter.summary(args.path)
    print(json.dumps(result,allow_nan=False))


if __name__=='__main__': main()
