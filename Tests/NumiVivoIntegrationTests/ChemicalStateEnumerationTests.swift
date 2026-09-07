import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ChemicalStateEnumerationTests {
    private func structure() -> VivoMolecularStructure {
        let atoms:[VivoMolecularAtom]=[
            .init(index:0,name:"O1",element:.init(atomicNumber:8,symbol:"O"),residueIndex:0),
            .init(index:1,name:"H11",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0),
            .init(index:2,name:"H12",element:.init(atomicNumber:1,symbol:"H"),residueIndex:0),
            .init(index:3,name:"O2",element:.init(atomicNumber:8,symbol:"O"),residueIndex:1),
            .init(index:4,name:"H21",element:.init(atomicNumber:1,symbol:"H"),residueIndex:1),
            .init(index:5,name:"H22",element:.init(atomicNumber:1,symbol:"H"),residueIndex:1)]
        return .init(identifier:"two-independent-water-sites",atoms:atoms,
            bonds:[.init(atomA:0,atomB:1),.init(atomA:0,atomB:2),
                   .init(atomA:3,atomB:4),.init(atomA:3,atomB:5)],
            residues:[.init(index:0,name:"W1",atomIndices:[0,1,2]),
                      .init(index:1,name:"W2",atomIndices:[3,4,5])],
            conformers:[.init(positionsNM:[.zero,.init(0.096,0,0),.init(-0.025,0.093,0),
                                          .init(0.5,0,0),.init(0.596,0,0),.init(0.475,0.093,0)])])
    }
    private func evidence(_ value:String)->VivoKineticEvidence {
        .init(source:"explicit synthetic enumeration rule",locator:value)
    }
    private func request(maximumStates:Int=16) -> VivoChemicalStateEnumerationRequest {
        let a=VivoChemicalStateSite(identifier:"siteA",variants:[
            .init(identifier:"neutral",boundProtonOffset:0,origin:.assumed,evidence:evidence("A-neutral")),
            .init(identifier:"protonated",boundProtonOffset:1,
                  addHydrogens:[.init(identifier:"HA3",name:"HA3",parent:0,angleReference:1,planeReference:2,
                    lengthNM:0.1,angleRadians:1.9,dihedralRadians:1.2)],
                  chargeEdits:[.init(sourceAtom:0,formalCharge:1)],origin:.assumed,evidence:evidence("A-protonated"))])
        let b=VivoChemicalStateSite(identifier:"siteB",variants:[
            .init(identifier:"neutral",boundProtonOffset:0,origin:.assumed,evidence:evidence("B-neutral")),
            .init(identifier:"protonated",boundProtonOffset:1,
                  addHydrogens:[.init(identifier:"HB3",name:"HB3",parent:3,angleReference:4,planeReference:5,
                    lengthNM:0.1,angleRadians:1.9,dihedralRadians:1.2)],
                  chargeEdits:[.init(sourceAtom:3,formalCharge:1)],origin:.assumed,evidence:evidence("B-protonated"))])
        return .init(identifier:"two-site-cartesian",structure:structure(),
            protonationSourceIdentifier:"synthetic supplied rules",referencePH:7,
            sites:[a,b],maximumStates:maximumStates)
    }

    @Test func enumeratesCartesianProductAndSingleSiteAdjacency() throws {
        let source=request(),result=try VivoChemicalStateEnumeration.calculate(source)
        try VivoChemicalStateEnumeration.validate(result,request:source)
        #expect(result.combinationCount == 4)
        #expect(result.catalog.states.count == 4)
        #expect(Set(result.catalog.states.map(\.boundProtonOffset)) == Set([0,1,2]))
        for state in result.catalogRequest.entries {
            #expect(state.neighbors.count == 2)
            for neighbor in state.neighbors {
                let left=Dictionary(uniqueKeysWithValues:state.identifier.split(separator:";").map{
                    let pair=$0.split(separator:"=",maxSplits:1);return(String(pair[0]),String(pair[1]))})
                let right=Dictionary(uniqueKeysWithValues:neighbor.split(separator:";").map{
                    let pair=$0.split(separator:"=",maxSplits:1);return(String(pair[0]),String(pair[1]))})
                #expect(left.keys.filter{left[$0] != right[$0]}.count == 1)
            }
        }
        let doubly=try #require(result.catalog.states.first{$0.boundProtonOffset==2})
        #expect(doubly.hydrogenCount == 6) // two waters (four H) plus two declared protons
        #expect(doubly.expectedFormalCharge == 2)
    }

    @Test func declaredCapacityRejectsCombinatorialExplosion() throws {
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateEnumeration.calculate(request(maximumStates:3)) }
    }

    @Test func conflictingLocalEditsRejectInsteadOfLastWriteWins() throws {
        var source=request()
        source.sites[1].variants[1].chargeEdits=[.init(sourceAtom:0,formalCharge:1)]
        #expect(throws:(any Error).self) { _ = try VivoChemicalStateEnumeration.calculate(source) }
    }
}