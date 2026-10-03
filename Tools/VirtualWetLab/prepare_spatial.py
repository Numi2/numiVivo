#!/usr/bin/env python3
"""Register a numerical spatial specimen using the native owners, no biological claims."""
import argparse
from pathlib import Path
import shutil
from wetlab import read,write,sha,native,require
from spatial_adapter import FORMAT
p=argparse.ArgumentParser(description=__doc__)
for name in ('mesh','material','plan','matter-binary','matter-shader','spatial-binary','partition-shader','output'):
    p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--observations',type=Path)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
for key,name in [('mesh','mesh.json'),('material','material.nmatter'),('plan','plan.json')]:shutil.copy2(getattr(a,key),a.output/name)
native(a.matter_binary,[a.output/'mesh.json',a.output/'material.nmatter',a.output/'geometry.json'],a.output/'geometry-log.json')
s={'id':read(a.plan)['specimenID']}
for key in ('geometry','mesh','material','plan'):
    name='material.nmatter' if key=='material' else key+'.json';s[key]=name;s[key+'SHA256']=sha(a.output/name)
if a.observations:
    o=read(a.observations);require(o['evidenceClass'] in ('measured','numerical-reference'),'Observation evidence class')
    require(o['geometrySHA256']==s['geometrySHA256'],'Observation geometry does not match native accepted state')
    shutil.copy2(a.observations,a.output/'observations.json')
    s.update(observations='observations.json',observationsSHA256=sha(a.output/'observations.json'),observationEvidenceClass=o['evidenceClass'])
config={'format':FORMAT,'id':'spatial-passive-tracer','title':'Spatial tissue transport','family':'spatial-tissue',
        'specimens':[s],'provenance':'Authored deformable scaffold, hypothetical diffusion and cell permeability; source-qualified anatomy and biological calibration unavailable.',
        'runtime':{'matterBinary':str(a.matter_binary.resolve()),'matterShader':str(a.matter_shader.resolve()),
                   'spatialBinary':str(a.spatial_binary.resolve()),'partitionShader':str(a.partition_shader.resolve())}}
write(a.output/'assay.json',config);print(a.output/'assay.json')
