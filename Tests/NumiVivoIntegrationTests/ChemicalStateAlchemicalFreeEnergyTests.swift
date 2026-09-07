import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ChemicalStateAlchemicalFreeEnergyTests {
    private func evidence(_ locator:String) throws -> VivoKineticEvidence {
        let snapshot = try VivoCanonicalJSON.fingerprint(Data(("synthetic mathematical fixture; " + locator).utf8))
        return .init(source:"synthetic alchemical chemical-state fixture", locator:locator, sourceFingerprint:snapshot.hex)
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
            samplingOrigin:.calculated, samplingEvidence:try evidence("same-proton"))
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
            samplingOrigin:.calculated, samplingEvidence:try evidence("missing-reservoir"))
        #expect(throws:(any Error).self) {
            _ = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        }
    }

    @Test func explicitReservoirPromotesPhysicalMBARToReferencePHSemigrandFreeEnergy() throws {
        let calibration = VivoProtonReservoirCalibration(
            referencePH:7, temperatureK:300, contributionPerBoundProtonKJPerMol:5.0,
            standardDeviationKJPerMol:0.4, origin:.fitted,
            evidence:try evidence("reservoir-calibration"))
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(
            identifier:"acid-pair-calibrated", temperatureK:300, referencePH:7,
            states:[.init(identifier:"deprotonated", boundProtonOffset:0),
                    .init(identifier:"protonated", boundProtonOffset:1)],
            samples:samples(offsets:[0,0]), mbar:mbar(), protonReservoir:calibration,
            samplingOrigin:.calculated, samplingEvidence:try evidence("physical-mbar"))
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
            protonReservoir:.init(referencePH:7, temperatureK:300, contributionPerBoundProtonKJPerMol:0,
                                  standardDeviationKJPerMol:nil, origin:.assumed,
                                  evidence:try evidence("zero-reference-calibration")),
            samplingOrigin:.calculated, samplingEvidence:try evidence("degenerate-physical"))
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
            protonReservoir:.init(referencePH:6, temperatureK:300, contributionPerBoundProtonKJPerMol:0,
                                  origin:.assumed, evidence:try evidence("wrong-pH")),
            samplingOrigin:.calculated, samplingEvidence:try evidence("mismatch"))
        #expect(throws:(any Error).self) {
            _ = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        }
    }

    @Test func adapterRejectsThermodynamicContextTransplantation() throws {
        let request = VivoChemicalStateAlchemicalFreeEnergyRequest(identifier: "bound-context", temperatureK: 300, referencePH: 7,
            states: [.init(identifier: "A", boundProtonOffset: 0), .init(identifier: "B", boundProtonOffset: 0)],
            samples: samples(offsets: [0,1]), mbar: mbar(), samplingOrigin: .calculated,
            samplingEvidence: try evidence("bound-context"))
        let result = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        #expect(result.temperatureK == 300 && result.referencePH == 7 && result.origin == .calculated)
        #expect(throws: (any Error).self) {
            _ = try result.thermodynamicsRequest(identifier: "wrong-temperature", temperatureK: 310, referencePH: 7, targetPH: 8)
        }
        #expect(throws: (any Error).self) {
            _ = try result.thermodynamicsRequest(identifier: "wrong-reference", temperatureK: 300, referencePH: 6, targetPH: 8)
        }
    }

    @Test func unknownComponentUncertaintyAndAssumptionsAreNotUpgraded() throws {
        var request = VivoChemicalStateAlchemicalFreeEnergyRequest(identifier: "unknown-reservoir-SD", temperatureK: 300, referencePH: 7,
            states: [.init(identifier: "A", boundProtonOffset: 0), .init(identifier: "B", boundProtonOffset: 1)],
            samples: samples(offsets: [0,0]), mbar: mbar(),
            protonReservoir: .init(referencePH: 7, temperatureK: 300, contributionPerBoundProtonKJPerMol: 0,
                standardDeviationKJPerMol: nil, origin: .assumed, evidence: try evidence("assumed-reservoir")),
            samplingOrigin: .calculated, samplingEvidence: try evidence("sampled-MBAR"))
        let result = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        #expect(result.estimates[1].conditionalStandardDeviationKJPerMol == nil)
        #expect(result.origin == .assumed)
        let adapter = try result.thermodynamicsRequest(identifier: "retained-assumption", temperatureK: 300, referencePH: 7, targetPH: 8)
        #expect(adapter.states.allSatisfy { $0.origin == .assumed })
        request.protonReservoir?.standardDeviationKJPerMol = 0.4
        request.protonReservoir?.origin = .fitted
        request.mbar.bootstrapReplicates = 0
        let noBootstrap = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request)
        #expect(noBootstrap.estimates[1].conditionalStandardDeviationKJPerMol == nil)
        request.protonReservoir?.temperatureK = 310
        #expect(throws: (any Error).self) { _ = try VivoChemicalStateAlchemicalFreeEnergy.calculate(request) }
        request.protonReservoir?.temperatureK = 300
        request.samplingOrigin = .assumed
        #expect(try VivoChemicalStateAlchemicalFreeEnergy.calculate(request).origin == .assumed)
    }
}
