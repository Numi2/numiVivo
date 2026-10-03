#!/usr/bin/env python3
"""Import measured spatial RNA/protein tables into the shared native contract.

CSV columns: entity_id,x,y,biological_unit,<feature...>. JSON manifest declares
modality, units, count/value semantics and coordinate system. Blank values remain
missing; no feature is imputed or reclassified as measured. No prediction occurs.
"""
import argparse,csv,math
from pathlib import Path
from wetlab import read,write,sha,require,native

def ingest(table,manifest,output,binary):
    table=Path(table);m=read(manifest);require(m['modality'] in ('rna','protein'),'RNA or protein modality required')
    require(m['state']=='MEASURED' and m['sourceURI'] and m['featureIDs'],'Measured source provenance and named panel required')
    with table.open(newline='') as stream:
        reader=csv.DictReader(stream);require(set(reader.fieldnames or [])=={'entity_id','x','y','biological_unit',*m['featureIDs']},'CSV axis differs from declared panel')
        rows=list(reader)
    require(0<len(rows)<=100_000,'Ingestion cell budget')
    ids=[r['entity_id'] for r in rows];require(len(set(ids))==len(ids),'Duplicate entity identity')
    e={'state':'MEASURED','sourceID':m['id'],'modelID':None,'assumptions':m['assumptions'],'validationDomain':None,'uncertainty':m['uncertainty'],'heldOut':m['heldOut']}
    entities=[];offsets=[0];columns=[];values=[]
    for r in rows:
        xy=[float(r['x']),float(r['y'])];require(all(math.isfinite(v) for v in xy),'Nonfinite coordinate')
        entities.append({'id':r['entity_id'],'level':'cell','parentID':None,'memberIDs':[],'position':xy,'biologicalUnitID':r['biological_unit'],'annotations':{},'evidence':e})
        for i,name in enumerate(m['featureIDs']):
            if r[name]=='':continue
            v=float(r[name]);require(math.isfinite(v),'Nonfinite measurement');columns.append(i);values.append(v)
        offsets.append(len(values))
    specimen={'format':'numivivo-tissue-specimen/v1','id':m['id'],'title':m['title'],'organism':m['organism'],
        'coordinateSystem':m['coordinateSystem'],'sources':[{'id':m['id'],'uri':m['sourceURI'],'sha256':sha(table),'description':m['description']}],
        'entities':entities,'measurements':[{'id':m['id']+'-measurements','modality':m['modality'],'unit':m['unit'],'timepoint':m['timepoint'],'entityIDs':ids,'featureIDs':m['featureIDs'],'rowOffsets':offsets,'featureIndices':columns,'values':values,'absentValue':'missing','evidence':e}], 'limitations':m['limitations']}
    output=Path(output);output.mkdir(parents=True,exist_ok=False);write(output/'specimen.json',specimen)
    native(binary,[output/'specimen.json'],output/'validation.json')
    write(output/'provenance.json',{'sourceSHA256':sha(table),'manifestSHA256':sha(manifest),'importerSHA256':sha(__file__),'nativeValidatorSHA256':sha(binary),'imputedValues':0})
    return output/'specimen.json'
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('table',type=Path);p.add_argument('manifest',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--binary',type=Path,required=True);a=p.parse_args();print(ingest(a.table,a.manifest,a.output,a.binary))
