import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ChemicalStateCatalogTests {
    private func water() -> VivoMolecularStructure {
        let atoms:[VivoMolecularAtom]=[
            .init(index:0,name:"O",element:.init(atomicNumber:8,symbol:"O"),residueIndex:0),
            .init(index:1,name:"H1",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0),
            .init(index:2,name:"H2",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0)]
        return .init(identifier:"catalog-water",atoms:atoms,
            bonds:[.init(atomA:0,atomB:1),.init(atomA:0,atomB:2)],
            residues:[.init(index:0,name:"HOH",atomIndices:[0,1,2])],
            conformers:[.init(positionsNM:[.zero,.init(0.096,0,0),.init(-0.025,0.093,0)])])
    }
    private func evidence(_ name:String)->VivoKineticEvidence {
        .init(source:"synthetic explicit state catalog",locator:name)
    }
    private func request() -> VivoChemicalStateCatalogRequest {
        let source=water()
        let water=VivoMolecularPreparationRequest(structure:source,microstateIdentifier:"water",
            protonationSourceIdentifier:"synthetic explicit catalog",pH:7,expectedFormalCharge:0)
        let hydronium=VivoMolecularPreparationRequest(structure:source,microstateIdentifier:"hydronium",
            protonationSourceIdentifier:"synthetic explicit catalog",pH:7,expectedFormalCharge:1,
            addHydrogens:[.init(identifier:"H3",name:"H3",parent:0,angleReference:1,planeReference:2,
                lengthNM:0.1,angleRadians:1.9,dihedralRadians:1.2)],
            chargeEdits:[.init(sourceAtom:0,formalCharge:1)],residueEdits:[.init(sourceResidue:0,name:"H3O")])
        return .init(identifier:"water-protonation-catalog",entries:[
            .init(identifier:"water",boundProtonOffset:0,preparation:water,neighbors:["hydronium"],
                  stateOrigin:.assumed,stateEvidence:evidence("water")),
            .init(identifier:"hydronium",boundProtonOffset:1,preparation:hydronium,neighbors:["water"],
                  stateOrigin:.assumed,stateEvidence:evidence("hydronium"))])
    }

    @Test func materializesPersistentHeavyAtomsAndProtonStoichiometry() throws {
        let source=request(),result=try VivoChemicalStateCatalog.calculate(source)
        try VivoChemicalStateCatalog.validate(result,request:source)
        #expect(result.states.map(\.identifier) == ["water","hydronium"])
        #expect(result.states.map(\.hydrogenCount) == [2,3])
        #expect(result.states.map(\.expectedFormalCharge) == [0,1])
        #expect(result.states[0].preparationResult.sourceToPrepared == [0,1,2])
        #expect(result.states[1].preparationResult.sourceToPrepared == [0,1,2])
        #expect(result.states[1].preparationResult.preparedToSource == [0,1,2,nil])
        #expect(result.states[0].preparedStructureFingerprint != result.states[1].preparedStructureFingerprint)
    }

    @Test func rejectsForgedProtonOffsetChargeAndSourceIdentity() throws {
        var badOffset=request();badOffset.entries[1].boundProtonOffset=0
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateCatalog.calculate(badOffset) }
        var badCharge=request();badCharge.entries[1].preparation.expectedFormalCharge=0
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateCatalog.calculate(badCharge) }
        var badSource=request();badSource.entries[1].preparation.structure.identifier="another-source"
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateCatalog.calculate(badSource) }
    }

    @Test func rejectsNonreciprocalCatalogGraph() throws {
        var source=request();source.entries[1].neighbors=["hydronium"]
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateCatalog.calculate(source) }
    }
}