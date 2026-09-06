import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ReactionConnectivityTests {
    private func rejects(_ work: () throws -> Void) {
        do { try work(); Issue.record("Expected an explicit connectivity contract rejection") }
        catch is VivoChemistryError { }
        catch { Issue.record("Unexpected error type: \(error)") }
    }
    private func retain<T: Encodable>(_ value: T, _ filename: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let directory = URL(fileURLWithPath: path, isDirectory: true)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: directory.appendingPathComponent(filename), options: .atomic)
    }
    @Test func properRotationPreservesMappedChirality() throws {
        let coordinates: [SIMD3<Double>] = [.init(0,0,0), .init(1,0,0), .init(0,2,0), .init(0,0,3)]
        let masses: [Double] = [1,2,3,4]
        let rotated = coordinates.map { SIMD3<Double>(-$0.y + 5, $0.x - 2, $0.z + 7) }
        let mirrored = coordinates.map { SIMD3<Double>(-$0.x, $0.y, $0.z) }
        #expect(try VivoMappedGeometry.properRotationRMSD(coordinates, rotated, masses: masses) < 1e-7)
        #expect(try VivoMappedGeometry.properRotationRMSD(coordinates, mirrored, masses: masses) > 0.1)
        var exchanged = coordinates; exchanged.swapAt(0, 1)
        #expect(try VivoMappedGeometry.properRotationRMSD(coordinates, exchanged, masses: masses) > 0.1)
    }
    @Test func linearAndAtomicAlignmentRemainWellDefined() throws {
        let line: [SIMD3<Double>] = [.init(-1,0,0), .init(0,0,0), .init(2,0,0)]
        let moved = line.map { SIMD3<Double>(3, $0.x + 2, -4) }
        #expect(try VivoMappedGeometry.properRotationRMSD(line, moved, masses: [1,2,1]) < 1e-7)
        #expect(try VivoMappedGeometry.properRotationRMSD([.zero], [.init(100,200,300)], masses: [1]) == 0)
        rejects { _ = try VivoMappedGeometry.properRotationRMSD([], [], masses: []) }
        rejects { _ = try VivoMappedGeometry.properRotationRMSD(line, moved, masses: [1,-1,1]) }
        rejects { _ = try VivoMappedGeometry.properRotationRMSD(line, [.zero], masses: [1,1,1]) }
    }
    @Test func refinementContractRejectsUnboundedOrUnresolvedProfiles() throws {
        try VivoReactionConnectivityConfiguration().validate()
        rejects { var cfg = VivoReactionConnectivityConfiguration(); cfg.stepMassWeighted = 0; try cfg.validate() }
        rejects { var cfg = VivoReactionConnectivityConfiguration(); cfg.comparisonSamples = 1; try cfg.validate() }
        rejects { var cfg = VivoReactionConnectivityConfiguration(); cfg.endpointEnergyToleranceHartree = .nan; try cfg.validate() }
        rejects { var cfg = VivoReactionConnectivityConfiguration(); cfg.maximumArcMassWeighted = 0.02; try cfg.validate() }
        rejects { var cfg = VivoReactionConnectivityConfiguration(); cfg.stepMassWeighted = 1e-8; try cfg.validate() }
    }
    @Test func mappedH3BranchesHaveIndependentRefinementEvidence() throws {
        func qualify(_ name: String) throws -> VivoNuclearQualifiedPoint {
            let calculation = try VivoReactionQualificationWorkflow.template(name).calculation
            guard case .qualify(let request) = calculation else {
                throw VivoChemistryError.invalid("reaction conformance template is not a nuclear request")
            }
            return try VivoNuclearQualification.run(request)
        }
        let saddle = try qualify("h3-saddle"), h2 = try qualify("h2-minimum"), atom = try qualify("h-atom")
        let left = VivoMappedReactionEndpoint(identifier: "H0-H1_plus_H2", components: [
            .init(atomIndices: [0,1], point: h2), .init(atomIndices: [2], point: atom)])
        let right = VivoMappedReactionEndpoint(identifier: "H0_plus_H1-H2", components: [
            .init(atomIndices: [0], point: atom), .init(atomIndices: [1,2], point: h2)])
        let cfg = VivoReactionConnectivityConfiguration(initialDisplacementMassWeighted: 0.015,
            stepMassWeighted: 0.04, maximumArcMassWeighted: 20, endpointMaximumGradient: 1e-7,
            endpointRMSDBohr: 0.02, endpointEnergyToleranceHartree: 1e-5,
            minimumIntercomponentSeparationBohr: 8, comparisonSamples: 128,
            minimumComparisonArcMassWeighted: 1, profileDistanceToleranceBohr: 0.005,
            profileEnergyToleranceHartree: 2e-5, maximumDescentElectronicEvaluations: 400000)
        let request = VivoReactionConnectivityRequest(atomIdentifiers: ["H0", "H1", "H2"],
            saddle: saddle, endpoints: [left, right], configuration: cfg)
        try retain(request, "mapped-h3-request.json")
        let result = try VivoReactionConnectivity.run(request)
        try retain(result, "mapped-h3-connectivity.json")
        try #require(result.converged)
        #expect(result.trials.count == 4)
        #expect(result.comparisons.count == 12)
        #expect(result.comparisons.allSatisfy { $0.passed })
        #expect(result.trials.allSatisfy {
            $0.reverseAssignment.endpointIdentifier != $0.forwardAssignment.endpointIdentifier
        })
        #expect(Set(result.trials.map(\.displacementScale)) == Set([1.0, 0.5]))
        #expect(Set(result.trials.map(\.stepScale)) == Set([1.0, 0.5]))
        #expect(result.descentElectronicEvaluations <= cfg.maximumDescentElectronicEvaluations)

        // One extra replay exercises the actual scientific gate without making
        // the test perform several redundant full connectivity reconstructions.
        let tstRequest = VivoTransitionStateTheoryRequest(connectivity: result,
            reactantEndpointIdentifier: left.identifier,
            transmission: .init(coefficient: 1, kind: .assumed,
                                sourceIdentifier: "classical no-recrossing integration test"))
        let tst = try VivoTransitionStateTheory.estimate(tstRequest)
        #expect(tst.rateConstant.isFinite && tst.rateConstant > 0)
        #expect(tst.molecularity == 2)
        #expect(tst.rateUnits == "Pa^-1 s^-1")
        #expect(!tst.hasIndependentTransmissionEvidence)
        #expect(tst.barrier.reactantMolecularity == 2)
        #expect(abs(tst.activationExponent + tst.barrier.activationGibbsHartree /
                    (VivoNuclearUnits.kHartree * tst.temperatureK)) < 1e-12)
        try retain(tst, "mapped-h3-assumed-tst.json")
        let dynamicEvidence = VivoTransmissionCoefficientEvidence(coefficient: 0.8, kind: .computedDynamics,
            sourceIdentifier: "synthetic recrossing-control identifier")
        try dynamicEvidence.validate()
        #expect(dynamicEvidence.kind == .computedDynamics && dynamicEvidence.coefficient == 0.8)

        let invalidMap = VivoMappedReactionEndpoint(identifier: "duplicate_atom", components: [
            .init(atomIndices: [0,1], point: h2), .init(atomIndices: [1], point: atom)])
        rejects {
            try VivoReactionConnectivity.validateRequest(.init(atomIdentifiers: ["H0","H1","H2"],
                saddle: saddle, endpoints: [left, invalidMap], configuration: cfg))
        }
        rejects {
            try VivoReactionConnectivity.validateRequest(.init(atomIdentifiers: ["H0","H0","H2"],
                saddle: saddle, endpoints: [left, right], configuration: cfg))
        }
        var tooSmall = cfg; tooSmall.maximumDescentElectronicEvaluations = 4
        rejects {
            try VivoReactionConnectivity.validateRequest(.init(atomIdentifiers: ["H0","H1","H2"],
                saddle: saddle, endpoints: [left, right], configuration: tooSmall))
        }
        let different = VivoReactionConnectivityRequest(atomIdentifiers: ["A","B","C"],
            saddle: saddle, endpoints: [left, right], configuration: cfg)
        rejects { try VivoReactionConnectivity.validate(result, request: different) }
        let forgedConnectivity = VivoReactionConnectivityResult(schema: result.schema, request: result.request,
            trials: result.trials, comparisons: result.comparisons, descentElectronicEvaluations: result.descentElectronicEvaluations,
            converged: false, interpretation: result.interpretation)
        rejects { _ = try VivoTransitionStateTheory.estimate(.init(connectivity: forgedConnectivity,
            reactantEndpointIdentifier: left.identifier)) }
        rejects { _ = try VivoTransitionStateTheory.estimate(.init(connectivity: result,
            reactantEndpointIdentifier: "unknown-endpoint")) }
    }
}
