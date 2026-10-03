#!/usr/bin/env python3
"""Native RNA/protein table admission and missingness regression, not biological validation."""
import argparse,json
from pathlib import Path
from import_spatial_table import ingest
from wetlab import read,write
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
manifest={'id':'protein-fixture','title':'Importer qualification fixture','organism':'fixture','modality':'protein','unit':'instrument-intensity','timepoint':'endpoint','state':'MEASURED','sourceURI':'urn:qualification:fixture','featureIDs':['A','B'],'coordinateSystem':{'id':'microscope','unit':'micrometre','axes':['x','y'],'description':'Declared fixture coordinates','micrometresPerUnit':1},'description':'Authored fixture; no biological evidence','assumptions':['Fixture only'],'uncertainty':'not quantified','heldOut':True,'limitations':['No real proteomics qualification']}
table=a.output/'table.csv';table.write_text('entity_id,x,y,biological_unit,A,B\nc1,1,2,sample,0,\nc2,3,4,sample,5,2\n');write(a.output/'manifest.json',manifest)
path=ingest(table,a.output/'manifest.json',a.output/'protein',a.binary);specimen=read(path);m=specimen['measurements'][0]
assert m['absentValue']=='missing' and m['values']==[0,5,2] and m['rowOffsets']==[0,1,3] and m['evidence']['state']=='MEASURED'
assert specimen['entities'][1]['position']==[3,4]
manifest.update(modality='rna',unit='umiCount');write(a.output/'rna-manifest.json',manifest)
ingest(table,a.output/'rna-manifest.json',a.output/'rna',a.binary)
table.write_text('entity_id,x,y,biological_unit,A,B\nc1,1,2,sample,0.5,\n')
try:ingest(table,a.output/'rna-manifest.json',a.output/'invalid-rna',a.binary)
except (ValueError,RuntimeError):pass
else:raise AssertionError('Fractional UMI admitted')
write(a.output/'checks.json',{'rnaAndProteinNativeAdmission':True,'missingAndZeroDistinct':True,'coordinatesExact':True,'fractionalUMIRejected':True,'biologicalQualification':False})
print(json.dumps(read(a.output/'checks.json')))
