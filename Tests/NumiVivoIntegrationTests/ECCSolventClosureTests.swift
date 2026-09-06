import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct ECCSolventClosureTests {
    private func rejects(_ body: () throws -> Void) {
        do { try body(); Issue.record("expected explicit ECC/solvent rejection") } catch {}
    }
    private func mutate<T: Codable>(_ value: T, _ edit: (inout [String:Any]) -> Void) throws -> T {
        var object = try JSONSerialization.jsonObject(with: JSONEncoder().encode(value)) as! [String:Any]
        edit(&object)
        return try JSONDecoder().decode(T.self, from: JSONSerialization.data(withJSONObject: object))
    }
    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let folder = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let root = URL(fileURLWithPath: folder, isDirectory: true)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: root.appendingPathComponent(name + ".json"), options: .atomic)
    }

    @Test func fullBathECCAndCorrelatedSolventReachOnePhysicalClosure() throws {
        let request = VivoECCSolventClosure.hydrogenControl()
        let result = try VivoECCSolventClosure.run(request)
        let reference = try VivoCorrelatedSolvation.solve(request.molecule)
        #expect(result.ecc.converged)
        #expect(result.history.count >= 2)
        #expect(result.aggregateECCFrameEvaluations <= request.configuration.maximumAggregateECCFrameEvaluations)
        #expect(result.state.alphaElectrons == 1 && result.state.betaElectrons == 1)
        #expect(result.state.orbitalCount == 2)
        #expect(result.globalSubspace.rows == 4 && result.globalSubspace.columns == 4)
        #expect(abs(result.occupations.reduce(0,+) - 2) < 1e-10)
        #expect(result.densityResidual <= request.molecule.configuration.densityTolerance)
        #expect(result.potentialResidualHartree <= request.molecule.configuration.potentialToleranceHartree)
        #expect(result.projectorResidual <= request.configuration.projectorTolerance)
        #expect(result.projectedResidualHartree <= request.configuration.projectedStationarityToleranceHartree)
        #expect(result.externalResidualHartree < 1e-8)
        #expect(abs(result.energyHartree - reference.energyHartree) < 1e-8)
        #expect(abs(result.gasEnergyHartree + result.equilibriumField.polarizationEnergyHartree - result.energyHartree) < 1e-12)
        #expect(try result.densityAO.adding(reference.totalDensityAO, scale: -1).frobeniusNorm < 1e-7)
        #expect(result.method.contains("one N-representable solvent density"))
        #expect(result.method.contains("not a proof of complete ECC-DMET/PCM"))
        try VivoECCSolventClosure.validate(result, request: request)
        try record(result, "ecc-solvent-full-bath-control")
    }

    @Test func closureCannotInventInactiveOccupationsOrStoredConvergence() throws {
        let request = VivoECCSolventClosure.hydrogenControl()
        let wrongCore = VivoECCSolventClosureRequest(molecule: request.molecule, embedding: request.embedding,
            inactiveDoublyOccupiedColumns: [[1]], globalSpace: request.globalSpace, configuration: request.configuration)
        rejects { _ = try VivoECCSolventClosure.run(wrongCore) }
        let tiny = VivoECCSolventClosureRequest(molecule: request.molecule, embedding: request.embedding,
            inactiveDoublyOccupiedColumns: request.inactiveDoublyOccupiedColumns,
            globalSpace: .init(maximumDimension: 1), configuration: request.configuration)
        rejects { _ = try VivoECCSolventClosure.run(tiny) }
        let result = try VivoECCSolventClosure.run(request)
        let forged = try mutate(result) { object in
            object["energyHartree"] = (object["energyHartree"] as! Double) - 0.1
            object["projectorResidual"] = 0.0
        }
        rejects { try VivoECCSolventClosure.validate(forged, request: request) }
        let insufficient = VivoECCSolventClosureRequest(molecule: request.molecule, embedding: request.embedding,
            inactiveDoublyOccupiedColumns: request.inactiveDoublyOccupiedColumns,
            globalSpace: request.globalSpace,
            configuration: .init(maximumOuterIterations: 2, projectorTolerance: 1e-14,
                                 projectedStationarityToleranceHartree: 1e-12,
                                 maximumAggregateECCFrameEvaluations: 4096))
        rejects { _ = try VivoECCSolventClosure.run(insufficient) }
    }
}
