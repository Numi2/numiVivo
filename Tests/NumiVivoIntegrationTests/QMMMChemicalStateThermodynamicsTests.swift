import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMChemicalStateThermodynamicsTests {
    private func fingerprint(_ text:String)throws->VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }

    private func state(_ id:String,_ protons:Int,_ freeEnergy:Double,_ evidenceID:String) throws -> VivoQMMMChemicalThermodynamicState {
        .init(identifier:id,boundProtonOffset:protons,relativeSemigrandFreeEnergyKJPerMol:freeEnergy,
              origin:.calculated,evidence:.init(source:"synthetic state free energy",locator:id,
                  sourceFingerprint:(try fingerprint(evidenceID)).hex))
    }

    @Test func onePHUnitProducesExactTenfoldProtonReservoirReweighting() throws {
        let states=[try state("protonated",1,0,"P"),try state("deprotonated",0,0,"D")]
        let reference=VivoQMMMChemicalStateThermodynamicsRequest(identifier:"reference",
            temperatureK:300,referencePH:7,targetPH:7,states:states)
        let atReference=try VivoQMMMChemicalStateThermodynamics.calculate(reference)
        #expect(abs((atReference.population(identifier:"protonated") ?? 0)-0.5)<1e-12)
        #expect(abs((atReference.population(identifier:"deprotonated") ?? 0)-0.5)<1e-12)
        try VivoQMMMChemicalStateThermodynamics.validate(atReference,request:reference)
        #expect(atReference.populationEvidence.sourceFingerprint==atReference.evidenceFingerprint.hex)

        var basic=reference;basic.targetPH=8
        let shifted=try VivoQMMMChemicalStateThermodynamics.calculate(basic)
        let protonated=try #require(shifted.population(identifier:"protonated"))
        let deprotonated=try #require(shifted.population(identifier:"deprotonated"))
        #expect(abs(deprotonated/protonated-10)<1e-11)
        #expect(abs(protonated-1.0/11.0)<1e-12)
        #expect(abs(deprotonated-10.0/11.0)<1e-12)
        #expect(shifted.mostPopulatedStateIdentifier=="deprotonated")
    }

    @Test func commonFreeEnergyAndProtonReferencesCancelExactly() throws {
        let base=[try state("A",1,-2,"A-base"),try state("B",0,3,"B-base")]
        let request=VivoQMMMChemicalStateThermodynamicsRequest(identifier:"base",temperatureK:310,
            referencePH:6,targetPH:8.5,states:base)
        let result=try VivoQMMMChemicalStateThermodynamics.calculate(request)

        let shifted=[try state("A",8,48,"A-shift"),try state("B",7,53,"B-shift")]
        let equivalent=VivoQMMMChemicalStateThermodynamicsRequest(identifier:"shifted",temperatureK:310,
            referencePH:6,targetPH:8.5,states:shifted)
        let shiftedResult=try VivoQMMMChemicalStateThermodynamics.calculate(equivalent)
        #expect(abs((result.population(identifier:"A") ?? -1)-(shiftedResult.population(identifier:"A") ?? -2))<1e-12)
        #expect(abs((result.population(identifier:"B") ?? -1)-(shiftedResult.population(identifier:"B") ?? -2))<1e-12)
    }

    @Test func equalProtonStoichiometryLeavesConformerRatioIndependentOfPH() throws {
        let states=[try state("conf-A",3,0,"conf-A"),try state("conf-B",3,2.5,"conf-B")]
        let acidic=try VivoQMMMChemicalStateThermodynamics.calculate(.init(identifier:"acidic",temperatureK:298.15,
            referencePH:7,targetPH:3,states:states))
        let basic=try VivoQMMMChemicalStateThermodynamics.calculate(.init(identifier:"basic",temperatureK:298.15,
            referencePH:7,targetPH:11,states:states))
        #expect(abs((acidic.population(identifier:"conf-A") ?? -1)-(basic.population(identifier:"conf-A") ?? -2))<1e-12)
        #expect(abs((acidic.population(identifier:"conf-B") ?? -1)-(basic.population(identifier:"conf-B") ?? -2))<1e-12)
    }

    @Test func calculatedStateRequiresImmutableEvidence() throws {
        let invalid=VivoQMMMChemicalThermodynamicState(identifier:"bad",boundProtonOffset:0,
            relativeSemigrandFreeEnergyKJPerMol:0,origin:.calculated,
            evidence:.init(source:"calculation",locator:"missing fingerprint"))
        let request=VivoQMMMChemicalStateThermodynamicsRequest(identifier:"invalid",temperatureK:300,
            referencePH:7,targetPH:7,states:[invalid])
        #expect(throws:(any Error).self) { _ = try VivoQMMMChemicalStateThermodynamics.calculate(request) }
    }
}
