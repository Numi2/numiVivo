#!/usr/bin/env python3
"""Qualify seals, exact native replay and observation-checkpoint recovery."""
import argparse
from pathlib import Path
import intervention_design as design
from wetlab import read,write,require

def run(campaign,workspace,output):
    c=design.catalog(campaign);spec=c['objectives'][0];selection={'context':c['contexts'][0]['id'],'objective':{k:spec[k] for k in ('genes','preserveGenes','penalty')}}
    p=design.preview(campaign,selection);require('no-intervention' in p['candidates'],'Missing reference option');require(not p['reliableWinner'],'False confident winner')
    try:design.preview(campaign,{'context':selection['context'],'objective':{**selection['objective'],'preserveGenes':['FOS']}});raise AssertionError('Invented receiver support')
    except ValueError:pass
    r=design.seal(campaign,workspace,selection);root=workspace/r['id'];require(not (root/'observations.safetensors').exists(),'Observation opened before seal')
    f=root/'prediction.json';original=f.read_bytes();f.write_bytes(original+b' ')
    try:design.check(root);raise AssertionError('Changed objective accepted')
    except ValueError:pass
    finally:f.write_bytes(original)
    original_write=design.write
    def interrupted(path,data):
        if Path(path).name=='comparison.json':raise RuntimeError('qualification interruption after observation snapshot')
        return original_write(path,data)
    design.write=interrupted
    try:design.reveal(root);raise AssertionError('Interruption missing')
    except RuntimeError:pass
    finally:design.write=original_write
    require((root/'observations.safetensors').exists() and not (root/'comparison.json').exists(),'Recovery checkpoint absent')
    original_observations=design.observations
    def forbidden(*args):raise AssertionError('Recovery re-extracted source outcomes')
    design.observations=forbidden
    try:revealed=design.reveal(root)
    finally:design.observations=original_observations
    verification=design.verify(root);require(not revealed['comparison']['biologicalPromotion'],'Replay promoted biology')
    write(output,{'status':'passed','record':r['id'],'checks':['no-intervention candidate','unqualified ranking remains explicit','unsupported preservation rejected','no observation before seal','changed objective rejected','interrupted reveal retains source snapshot','recovery does not re-extract or refit','all three native models replay bit-exactly','selection metrics reconstruct'],'replay':verification,'scope':'software qualification on already exposed evaluation cohort'})
    print('PASS',r['id'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--campaign',type=Path,required=True);p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();run(a.campaign,a.workspace,a.output)
