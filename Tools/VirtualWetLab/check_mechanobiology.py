#!/usr/bin/env python3
"""Compile the exchange contract as an external Swift client and check admission."""
import argparse,subprocess,json
from pathlib import Path
from wetlab import sha,write
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);out=a.output.resolve();root=Path(__file__).resolve().parents[2]
sources=[root/'Sources/NumiVivoKit/Tissue'/name for name in ['VivoTissueSpecimen.swift','VivoMechanobiologyExchange.swift']]
subprocess.run(['swiftc','-emit-module','-emit-library','-module-name','NumiTissueExchange',*[str(s) for s in sources],'-emit-module-path',str(out/'NumiTissueExchange.swiftmodule'),'-o',str(out/'libNumiTissueExchange.dylib')],check=True)
client=out/'client.swift';client.write_text('''import Foundation
import NumiTissueExchange
let id=UUID()
func evidence(_ state:VivoEvidenceState = .hypothesis)->VivoTissueEvidence {
    VivoTissueEvidence(state:state,sourceID:"fixture",modelID:"test-law",assumptions:["Not biological evidence"],validationDomain:nil,uncertainty:"unquantified",heldOut:false)
}
func proposal(direction:String="biologyToMechanics", values:[Double]=[1], state:VivoEvidenceState = .hypothesis)->VivoMechanobiologyProposal {
    let channel=VivoMechanobiologyChannel(identifier:"permeability",direction:direction,quantity:"permeability",unit:"m/s",entityIDs:["cell1"],components:["scalar"],values:values,modelFingerprint:"fixture-law",evidence:evidence(state))
    return VivoMechanobiologyProposal(transactionID:id,specimenFingerprint:"specimen",mechanicsCheckpoint:"m0",biologyCheckpoint:"b0",timeBeforeSeconds:0,proposedTimeAfterSeconds:1,channels:[channel],convergenceTolerance:1e-6)
}
func prepared(_ owner:String, checkpoint:String, time:Double=1, residual:Double=0)->VivoMechanobiologyPreparedState {
    VivoMechanobiologyPreparedState(transactionID:id,owner:owner,previousCheckpoint:checkpoint,candidateCheckpoint:owner+"1",acceptedTimeSeconds:time,converged:true,residual:residual)
}
var rejected=0
func reject(_ body:()throws->Void) { do {try body();fatalError("Invalid proposal accepted")} catch {rejected += 1} }
let p=proposal(),m=prepared("mechanics",checkpoint:"m0"),b=prepared("biology",checkpoint:"b0")
try p.validatePrepared(mechanics:m,biology:b)
reject {try p.validatePrepared(mechanics:prepared("mechanics",checkpoint:"stale"),biology:b)}
reject {try p.validatePrepared(mechanics:m,biology:prepared("biology",checkpoint:"b0",time:2))}
reject {try p.validatePrepared(mechanics:m,biology:prepared("biology",checkpoint:"b0",residual:1))}
reject {try proposal(values:[.nan]).validate()}
reject {try proposal(values:[]).validate()}
reject {try proposal(state:.unavailable).validate()}
reject {try proposal(direction:"mechanicsToBiology").validate()}
print("{\\"externalSwiftClient\\":true,\\"validProposalAdmitted\\":true,\\"invalidProposalsRejected\\":\\(rejected),\\"bidirectionalRuntimeExecuted\\":false}")
''')
subprocess.run(['swiftc',str(client),'-I',str(out),'-L',str(out),'-lNumiTissueExchange','-Xlinker','-rpath','-Xlinker',str(out),'-o',str(out/'check')],check=True)
result=json.loads(subprocess.check_output([str(out/'check')],text=True));assert result['invalidProposalsRejected']==7
write(out/'checks.json',{**result,'sourceSHA256':{s.name:sha(s) for s in sources},'clientSHA256':sha(client)})
print(json.dumps(result))
