import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct SolventPromotionTests {
    private func fixture() throws -> (VivoMolecularStructureDocument,VivoClassicalSystem,[VivoVector3D]) {
        let C=try #require(VivoElement.from(symbol:"C")),O=try #require(VivoElement.from(symbol:"O")),H=try #require(VivoElement.from(symbol:"H"))
        let atoms:[VivoMolecularAtom]=[
            .init(index:0,name:"C1",element:C,residueIndex:0),
            .init(index:1,name:"O",element:O,residueIndex:1,isHetero:true),
            .init(index:2,name:"H1",element:H,residueIndex:1,isHetero:true),
            .init(index:3,name:"H2",element:H,residueIndex:1,isHetero:true),
            .init(index:4,name:"O",element:O,residueIndex:2,isHetero:true),
            .init(index:5,name:"H1",element:H,residueIndex:2,isHetero:true),
            .init(index:6,name:"H2",element:H,residueIndex:2,isHetero:true)]
        let residues=[VivoMolecularResidue(index:0,name:"LIG",atomIndices:[0]),
            .init(index:1,name:"HOH",atomIndices:[1,2,3]),.init(index:2,name:"HOH",atomIndices:[4,5,6])]
        let positions:[VivoVector3D]=[.init(0.05,0,0),
            .init(1.95,0,0),.init(1.96,0.01,0),.init(1.94,-0.01,0),
            .init(1.0,1.0,1.0),.init(1.01,1.0,1.0),.init(0.99,1.0,1.0)]
        let structure=VivoMolecularStructure(identifier:"periodic-solvent-shell",atoms:atoms,
            bonds:[.init(atomA:1,atomB:2),.init(atomA:1,atomB:3),.init(atomA:4,atomB:5),.init(atomA:4,atomB:6)],
            residues:residues,conformers:[.init(positionsNM:positions)],
            periodicCell:.init(a:.init(2,0,0),b:.init(0,2,0),c:.init(0,0,2)))
        let document=try VivoMolecularStructureDocument(structure:structure)
        let particles=atoms.map { atom in VivoClassicalParticle(index:atom.index,atomIndex:atom.index,
            typeIdentifier:atom.element.symbol,massDa:atom.element.atomicNumber==1 ? 1.008:(atom.element.atomicNumber==8 ? 15.999:12.011),
            chargeE:0,sigmaNM:0,epsilonKJPerMol:0) }
        let system=VivoClassicalSystem(identifier:"periodic-solvent-shell",structureFingerprint:document.structureFingerprint,particles:particles)
        return (document,system,positions)
    }

    @Test func periodicWholeWaterPromotionUpdatesElectronSector() throws {
        let (document,system,positions)=try fixture()
        let base=VivoQMMMRegionRequest(qmAtomIndices:[0],alphaElectrons:3,betaElectrons:3)
        let periodic=try VivoQMMMCompiler.promoteSolventRegion(document:document,system:system,
            particlePositionsNM:positions,request:base,policy:.init(radiusNM:0.15))
        #expect(periodic.promotedResidues.map(\.residueIndex)==[1])
        #expect(periodic.promotedResidues[0].atomIndices==[1,2,3])
        #expect(abs(periodic.promotedResidues[0].minimumDistanceNM-0.1)<1e-12)
        #expect(periodic.qmAtomIndices==[0,1,2,3])
        #expect(periodic.addedAlphaElectrons==5 && periodic.addedBetaElectrons==5)
        #expect(periodic.request.alphaElectrons==8 && periodic.request.betaElectrons==8)
        #expect(periodic.periodicCellUsed && periodic.requiresUnwrappedFiniteClusterForQMMM)
        #expect(periodic.interpretation.contains("does not itself unwrap"))

        let nonperiodic=try VivoQMMMCompiler.promoteSolventRegion(document:document,system:system,
            particlePositionsNM:positions,request:base,
            policy:.init(radiusNM:0.15,periodicMinimumImage:false))
        #expect(nonperiodic.promotedResidues.isEmpty && nonperiodic.qmAtomIndices==[0])
        #expect(nonperiodic.addedAlphaElectrons==0 && nonperiodic.addedBetaElectrons==0)
    }

    @Test func capsPartialMoleculesAndOddElectronSolventFailExplicitly() throws {
        let (document,system,positions)=try fixture()
        let base=VivoQMMMRegionRequest(qmAtomIndices:[0],alphaElectrons:3,betaElectrons:3)
        #expect(throws: (any Error).self) {
            _=try VivoQMMMCompiler.promoteSolventRegion(document:document,system:system,particlePositionsNM:positions,
                request:base,policy:.init(radiusNM:0.15,maximumPromotedAtoms:2))
        }
        let partial=VivoQMMMRegionRequest(qmAtomIndices:[0,1],alphaElectrons:7,betaElectrons:7)
        #expect(throws: (any Error).self) {
            _=try VivoQMMMCompiler.promoteSolventRegion(document:document,system:system,particlePositionsNM:positions,
                request:partial,policy:.init(radiusNM:0.15))
        }
        var charged=document.structure
        charged.atoms[1].formalCharge=1
        let chargedDocument=try VivoMolecularStructureDocument(structure:charged)
        let chargedSystem=VivoClassicalSystem(identifier:"charged-water",structureFingerprint:chargedDocument.structureFingerprint,
            particles:system.particles)
        #expect(throws: (any Error).self) {
            _=try VivoQMMMCompiler.promoteSolventRegion(document:chargedDocument,system:chargedSystem,
                particlePositionsNM:positions,request:base,policy:.init(radiusNM:0.15))
        }
    }
}
