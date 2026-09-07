import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PlatformMechanismWorkflowTests {
    private func identity() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data("platform-mechanism-test-implementation".utf8))
    }
    private func request() -> VivoMechanismCandidateRequest {
        let carbon = VivoElement(atomicNumber: 6, symbol: "C")
        let atoms = (0..<4).map { VivoMolecularAtom(index: UInt32($0), name: "C\($0)", element: carbon) }
        let structure = VivoMolecularStructure(identifier: "platform-mechanism",
            atoms: atoms,
            bonds: [.init(atomA:0,atomB:1,order:.single),.init(atomA:2,atomB:3,order:.single)],
            conformers: [.init(positionsNM:[.zero,.init(0.15,0,0),.init(0.30,0,0),.init(0.45,0,0)])])
        return .init(identifier: "platform-mechanism", structure: structure,
            editableBonds: [
                .init(atomA:1,atomB:2,allowedTargets:[.single]),
                .init(atomA:2,atomB:3,allowedTargets:[.absent])
            ], maximumEditsPerCandidate: 2)
    }

    @Test func registeredOperationReconstructsAndRejectsTampering() async throws {
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let definition = try registry.definition("vivo.platform.mechanism-candidates")
        #expect(definition.summary.contains("hypotheses"))
        #expect(definition.validationScope.contains("scientific scope"))
        let operation = definition.operation
        let req = request()
        let inputs = ["request": try VivoCanonicalJSON.encode(req)]
        let outputs = try await operation.execute(.object([:]), inputs, .init())
        try operation.validateOutputs(.object([:]), inputs, outputs, .init())
        let result = try VivoCanonicalJSON.decode(VivoMechanismCandidateResult.self,
            from: try #require(outputs["candidates"]))
        #expect(result.candidates.count == 3)
        #expect(result.candidates.allSatisfy { $0.interpretation.contains("graph hypothesis only") })
        let direct = try VivoMechanismCandidateEnumeration.enumerate(req)
        #expect(result == direct)

        var tampered = outputs
        let replacement = VivoMechanismCandidateResult(schema: result.schema,
            requestFingerprint: result.requestFingerprint,
            sourceStructureFingerprint: result.sourceStructureFingerprint,
            candidates: [], rejectedByValence: result.rejectedByValence,
            interpretation: result.interpretation,
            evidenceFingerprint: result.evidenceFingerprint)
        tampered["candidates"] = try VivoCanonicalJSON.encode(replacement)
        #expect(throws: (any Error).self) {
            try operation.validateOutputs(.object([:]), inputs, tampered, .init())
        }
    }
}
