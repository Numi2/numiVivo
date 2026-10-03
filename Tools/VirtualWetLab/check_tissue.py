import json,subprocess,copy
from pathlib import Path
import argparse,tempfile
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);r=a.output;sha='0'*64
e={'state':'MEASURED','sourceID':'s','modelID':None,'assumptions':[],'validationDomain':None,'uncertainty':'unknown','heldOut':False}
s={'format':'numivivo-tissue-specimen/v1','id':'test','title':'test','organism':'NCBITaxon:10090','coordinateSystem':{'id':'xy','unit':'micrometre','axes':['x','y'],'description':'test','micrometresPerUnit':1},'sources':[{'id':'s','uri':'fixture://numeric','sha256':sha,'description':'synthetic fixture'}],'entities':[{'id':'cell','level':'cell','parentID':None,'memberIDs':[],'position':[0,0],'biologicalUnitID':'unit','annotations':{},'evidence':e}], 'measurements':[{'id':'rna','modality':'rna','unit':'umiCount','timepoint':'endpoint','entityIDs':['cell'],'featureIDs':['A'],'rowOffsets':[0,1],'featureIndices':[0],'values':[1],'absentValue':'zero','evidence':e}], 'limitations':['synthetic contract check only']}
cases=[]
def check(name,mutate,expect):
 d=copy.deepcopy(s);mutate(d);p=r/(name+'.json');p.write_text(json.dumps(d));v=subprocess.run([str(a.binary),str(p)],capture_output=True,text=True);assert (v.returncode==0)==expect,(name,v.stderr);cases.append({'name':name,'passed':True,'status':v.returncode})
check('valid',lambda d:None,True)
check('invalid-count',lambda d:d['measurements'][0].update(values=[-1]),False)
check('invalid-coordinate',lambda d:d['entities'][0].update(position=[0]),False)
check('invalid-cycle',lambda d:d['entities'][0].update(parentID='cell'),False)
check('false-validation',lambda d:d['measurements'][0]['evidence'].update(state='SIMULATED — VALIDATED DOMAIN'),False)
check('inferred-unidentified',lambda d:d['measurements'][0]['evidence'].update(state='MODEL INFERENCE'),False)
check('unknown-source',lambda d:d['sources'][0].update(id='missing'),False)
check('protein-values',lambda d:d['measurements'][0].update(modality='protein',unit='fluorescence-AU',values=[.25],absentValue='missing'),True)
(r/'native-contract-checks.json').write_text(json.dumps(cases,indent=2));print('8 native contract checks passed')
