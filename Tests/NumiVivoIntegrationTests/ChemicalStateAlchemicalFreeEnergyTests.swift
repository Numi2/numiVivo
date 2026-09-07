import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ChemicalStateAlchemicalFreeEnergyTests {
    private func evidence(_ locator:String) -> VivoKineticEvidence {
        .init(source:"synthetic alchemical chemical-state fixture", locator:locator)
    }

    private func samples(offsets:[Double], perState:Int = 32) -> [VivoMBARReducedPotentialSample] {
        var result:[VivoMBARReducedPotentialSample] = []
        for origin in offsets.indices {
            for index in 0..<perState {
                result.append(.init(identifier:"s-\(origin)-\(index)", originStateIndex:origin,
                                    reducedPotentials:offsets))
            }
        }
        return result
    }

    private func mbar() -> VivoMBARConfiguration {
        .init(residualTolerance:1e-12, maximumIterations:10_000,
              minimumBidirectionalAcceptanceOverlap:0.5,
              minimumTargetEffectiveSamples:10,
              bootstrapReplicates:8, bootstrapSeed:991,
              maximumWorkElements:1_000_000)
    }

    @Test func sameProtonStatesNeedNoReservoirCalibration() throws {
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"tautomer-pair", temperatureK:300, referencePH:7,
            states:[.init(identifier:"A", boundProtonOffset:0),
                    .init(identifier:"B", boundProtonOffset:0)],
            samples:samples(offsets:[0,1.25]), mbar:mbar(), protonReservoir:nil,
            samplingOrigin:.calculated, samplingEvidence:evidence("same-proton"))
        let result = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        #expect(result.mbar.converged)
        #expect(result.estimates.count == 2)
        let rt = 0.00831446261815324 * 300.0
        #expect(abs(result.estimates[1].relativePhysicalFreeEnergyKJPerMol - 1.25 * rt) < 1e-9)
        #expect(abs(result.estimates[1].relativeSemigrandFreeEnergyKJPerMolAtReferencePH - 1.25 * rt) < 1e-9)
        try VivoChemicalStateAlchemicalFreeEnergy.validate(result, request:request)
    }

    @Test func protonChangingStatesRejectMissingReservoirCalibration() throws {
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"acid-pair", temperatureK:300, referencePH:7,
            states:[.init(identifier:"deprotonated", boundProtonOffset:0),
                    .init(identifier:"protonated", boundProtonOffset:1)],
            samples:samples(offsets:[0,0]), mbar:mbar(), protonReservoir:nil,
            samplingOrigin:.calculated, samplingEvidence:evidence("missing-reservoir"))
        #expect(throws:(any Error).self) {
            _ = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        }
    }

    @Test func explicitReservoirPromotesPhysicalMBARToReferencePHSemigrandFreeEnergy() throws {
        let calibration = VivoProtonReservoirCalibration(
            referencePH:7, contributionPerBoundProtonKJPerMol:5.0,
            standardDeviationKJPerMol:0.4, origin:.fitted,
            evidence:evidence("reservoir-calibration"))
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"acid-pair-calibrated", temperatureK:300, referencePH:7,
            states:[.init(identifier:"deprotonated", boundProtonOffset:0),
                    .init(identifier:"protonated", boundProtonOffset:1)],
            samples:samples(offsets:[0,0]), mbar:mbar(), protonReservoir:calibration,
            samplingOrigin:.calculated, samplingEvidence:evidence("physical-mbar"))
        let result = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        #expect(abs(result.estimates[0].relativeSemigrandFreeEnergyKJPerMolAtReferencePH) < 1e-12)
        #expect(abs(result.estimates[1].relativeSemigrandFreeEnergyKJPerMolAtReferencePH - 5.0) < 1e-9)
        #expect(result.estimates[1].conditionalStandardDeviationKJPerMol != nil)
        #expect(abs(result.estimates[1].conditionalStandardDeviationKJPerMol! - 0.4) < 1e-9)
    }

    @Test func thermodynamicsAdapterProducesExpectedTenfoldPHShiftWhenReferenceStatesAreDegenerate() throws {
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"acid-degenerate", temperatureK:300, referencePH:7,
            states:[.init(identifier:"deprotonated", boundProtonOffset:0),
                    .init(identifier:"protonated", boundProtonOffset:1)],
            samples:samples(offsets:[0,0]), mbar:mbar(),
            protonReservoir:.init(referencePH:7, contributionPerBoundProtonKJPerMol:0,
                                  standardDeviationKJPerMol:nil, origin:.assumed,
                                  evidence:evidence("zero-reference-calibration")),
            samplingOrigin:.calculated, samplingEvidence:evidence("degenerate-physical"))
        let alchemical = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        let thermoRequest = try alchemical.thermodynamicsRequest(identifier:"pH8", temperatureK:300,
                                                                  referencePH:7, targetPH:8)
        let thermo = try VivoQMMMChemicalStateThermodynamics.calculate(thermoRequest)
        let deprotonated = try #require(thermo.population(identifier:"deprotonated"))
        let protonated = try #require(thermo.population(identifier:"protonated"))
        #expect(abs(deprotonated / protonated - 10.0) < 1e-8)
    }

    @Test func calibrationReferencePHMustMatchRequest() throws {
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"mismatch", temperatureK:300, referencePH:7,
            states:[.init(identifier:"A", boundProtonOffset:0),
                    .init(identifier:"B", boundProtonOffset:1)],
            samples:samples(offsets:[0,0]), mbar:mbar(),
            protonReservoir:.init(referencePH:6, contributionPerBoundProtonKJPerMol:0,
                                  origin:.assumed, evidence:evidence("wrong-pH")),
            samplingOrigin:.calculated, samplingEvidence:evidence("mismatch"))
        #expect(throws:(any Error).self) {
            _ = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        }
    }
}
