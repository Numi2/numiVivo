#!/usr/bin/env python3
"""Bounded build of the shared native tissue data contract (no model replacement)."""
import argparse,subprocess
from pathlib import Path
from wetlab import write,sha
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
root=Path(__file__).resolve().parents[2]
names=['Sources/NumiVivoKit/Tissue/VivoTissueSpecimen.swift','Sources/NumiVivoKit/Tissue/VivoMechanobiologyExchange.swift','Sources/NumiVivoTissueCLI/TissueCLI.swift']
subprocess.run(['swiftc','-O','-parse-as-library',*[str(root/n) for n in names],'-o',str(a.output/'numivivo-tissue')],check=True)
write(a.output/'receipt.json',{'sources':{n:sha(root/n) for n in names},'binarySHA256':sha(a.output/'numivivo-tissue'),'scope':'native tissue contract; not biological model validation'})
