#!/usr/bin/env python3
"""Build a bounded CLI from the owning native partition sources, with receipts.

The temporary Swift package copies exact owner sources; no solver is reimplemented.
The full suite also exposes this executable through its root Package.swift.
"""
import argparse, hashlib, json, shutil, subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=Path(__file__).resolve().parents[2];a.output.mkdir(parents=True,exist_ok=False)
files=['Sources/NumiVivoKit/Physiology/VivoPhysiologicalPartitionModel.swift',
       'Sources/NumiVivoKit/Physiology/VivoPhysiologicalPartitionRuntime.swift','Sources/NumiVivoSpatialCLI/SpatialCLI.swift',
       'Sources/NumiVivoShaders/ShaderResources.swift','Sources/NumiVivoShaders/Resources/NumiVivoPartitionKernels.metal']
hashes={}
for name in files:
    source=root/name;hashes[name]=hashlib.sha256(source.read_bytes()).hexdigest()
    target=a.output/('Sources/NumiVivoShaders/'+('/'.join(name.split('/')[2:])) if '/NumiVivoShaders/' in name else 'Sources/Spatial/'+source.name)
    target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
(a.output/'Package.swift').write_text('''// swift-tools-version: 6.3
import PackageDescription
let package=Package(name:"NumiVivoSpatial",platforms:[.macOS(.v15)],products:[.executable(name:"numivivo-spatial",targets:["Spatial"])],targets:[.target(name:"NumiVivoShaders",resources:[.copy("Resources/NumiVivoPartitionKernels.metal")]),.executableTarget(name:"Spatial",dependencies:["NumiVivoShaders"],linkerSettings:[.linkedFramework("Metal")])],swiftLanguageModes:[.v6])
''')
subprocess.run(['swift','build','-c','release','--package-path',str(a.output),'-j','3'],check=True)
binary=a.output/'.build/release/numivivo-spatial'
(a.output/'build-receipt.json').write_text(json.dumps({'sources':hashes,'binarySHA256':hashlib.sha256(binary.read_bytes()).hexdigest(),'scope':'exact partition owner source; scoped native CPU/Metal build'},indent=2)+'\n')
print(binary)
