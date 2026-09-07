import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct LatentProtonTopologyTests {
    private func source() -> VivoMolecularStructure {
        let atoms:[VivoMolecularAtom] = [
            .init(index:0,name:"C",element:.init(atomicNumber:6,symbol:"C"),residueIndex:0),
            .init(index:1,name:"O",element:.init(atomicNumber:8,symbol:"O"),residueIndex:0),
            .init(index:2,name:"H",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0)]
        return .init(identifier:"synthetic-acid", atoms:atoms,
            bonds:[.init(atomA:0,atomB:1),.init(atomA:1,atomB:2)],
            residues:[.init(index:0,name:"AH",atomIndices:[0,1,2])],
            conformers:[.init(identifier:"model",positionsNM:[.zero,.init(0.13,0,0),.init(0.225,0,0)])])
    }

    private func library() -> VivoForceFieldLibrary {
        let protonated = VivoForceFieldResidueTemplate(residueNames:["AH"],
            atoms:[.init(name:"C",typeIdentifier:"C",chargeE:0),
                   .init(name:"O",typeIdentifier:"O",chargeE:0),
                   .init(name:"H",typeIdentifier:"H",chargeE:0)],
            bonds:[.init(atomA:"C",atomB:"O"),.init(atomA:"O",atomB:"H")], impropers:[])
        let deprotonated = VivoForceFieldResidueTemplate(residueNames:["A"],
            atoms:[.init(name:"C",typeIdentifier:"C",chargeE:0),
                   .init(name:"O",typeIdentifier:"OM",chargeE:-1)],
            bonds:[.init(atomA:"C",atomB:"O")], impropers:[])
        return .init(identifier:"synthetic-latent-proton",version:"1",
            atomTypes:[.init(identifier:"C",elementAtomicNumber:6,massDa:12,sigmaNM:0.3,epsilonKJPerMol:0.1),
                       .init(identifier:"O",elementAtomicNumber:8,massDa:16,sigmaNM:0.28,epsilonKJPerMol:0.2),
                       .init(identifier:"OM",elementAtomicNumber:8,massDa:16,sigmaNM:0.29,epsilonKJPerMol:0.25),
                       .init(identifier:"H",elementAtomicNumber:1,massDa:1.008,sigmaNM:0.1,epsilonKJPerMol:0.01)],
            residueTemplates:[protonated,deprotonated],
            bondParameters:[.init(typeA:"C",typeB:"O",lengthNM:0.13,forceConstant:1000),
                            .init(typeA:"C",typeB:"OM",lengthNM:0.13,forceConstant:1000),
                            .init(typeA:"O",typeB:"H",lengthNM:0.095,forceConstant:1500)],
            angleParameters:[.init(typeA:"C",typeB:"O",typeC:"H",angleRadians:2.0,forceConstant:100)],
            provenance:["purpose":"synthetic latent proton topology test only"])
    }

    private func catalog() throws -> (VivoChemicalStateCatalogRequest,VivoChemicalStateCatalogResult) {
        let structure=source(), ff=library()
        let deprotonated=VivoMolecularPreparationRequest(structure:structure,microstateIdentifier:"A-",
            protonationSourceIdentifier:"synthetic explicit state",pH:7,expectedFormalCharge:-1,
            removeHydrogens:[2],chargeEdits:[.init(sourceAtom:1,formalCharge:-1)],
            residueEdits:[.init(sourceResidue:0,name:"A")],forceField:ff)
        let protonated=VivoMolecularPreparationRequest(structure:structure,microstateIdentifier:"AH",
            protonationSourceIdentifier:"synthetic explicit state",pH:7,expectedFormalCharge:0,forceField:ff)
        let evidence=VivoKineticEvidence(source:"synthetic latent topology fixture",locator:"explicit states")
        let request=VivoChemicalStateCatalogRequest(identifier:"acid-catalog",entries:[
            .init(identifier:"A-",boundProtonOffset:0,preparation:deprotonated,neighbors:["AH"],stateOrigin:.assumed,stateEvidence:evidence),
            .init(identifier:"AH",boundProtonOffset:1,preparation:protonated,neighbors:["A-"],stateOrigin:.assumed,stateEvidence:evidence)])
        return (request,try VivoChemicalStateCatalog.calculate(request))
    }

    @Test func compilerRetainsPositiveMassGhostAndCommonValenceAnchor() throws {
        let (request,catalog)=try catalog()
        let result=try VivoLatentProtonTopologyCompiler.compile(catalogRequest:request,catalog:catalog)
        #expect(result.states.count==2)
        #expect(result.latentProtons.count==1)
        let slot=try #require(result.latentProtons.first)
        #expect(slot.identity == .source(2))
        #expect(slot.anchorBondCount==1)
        #expect(slot.anchorAngleCount==1)
        let deprotonated=try #require(result.states.first(where:{$0.stateIdentifier=="A-"}))
        let protonated=try #require(result.states.first(where:{$0.stateIdentifier=="AH"}))
        #expect(deprotonated.ghostParticleIndices == [slot.unionParticleIndex])
        #expect(protonated.ghostParticleIndices.isEmpty)
        let ghost=deprotonated.system.particles[Int(slot.unionParticleIndex)]
        #expect(ghost.massDa > 0)
        #expect(ghost.chargeE == 0)
        #expect(ghost.sigmaNM == 0)
        #expect(ghost.epsilonKJPerMol == 0)
        #expect(deprotonated.system.bonds == protonated.system.bonds)
        #expect(deprotonated.system.angles == protonated.system.angles)
        #expect(deprotonated.system.structureFingerprint == protonated.system.structureFingerprint)
    }

    @Test func latentSystemsSatisfyConstantPHCommonManifoldFactory() throws {
        let (request,catalog)=try catalog()
        let result=try VivoLatentProtonTopologyCompiler.compile(catalogRequest:request,catalog:catalog)
        let config=VivoMDConfiguration(timeStepPS:0.001,cutoffNM:1,neighborSkinNM:0,
            electrostatics:.cutoff,ensemble:.nve,thermostat:.none,targetTemperatureK:nil,
            frictionPerPS:nil,neighborListEnabled:false,randomSeed:77)
        let specifications=result.states.map {
            VivoConstantPHMetalStateSpecification(identifier:$0.stateIdentifier,system:$0.system,configuration:config)
        }
        let executable=try VivoConstantPHMetalStateFactory.make(specifications:specifications,samplingTemperatureK:300)
        #expect(executable.count==2)
        #expect(Set(executable.map(\.physicalManifoldFingerprint)).count==1)
        #expect(Set(executable.map(\.hamiltonianFingerprint)).count==2)
    }
}
