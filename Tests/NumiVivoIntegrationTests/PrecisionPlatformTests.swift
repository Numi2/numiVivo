import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PrecisionPlatformTests {
    @Test func correctionWorkflowReconstructsAndRejectsAlteredOutputs() async throws {
        let identity = try PrecisionSamplingFixtures.id("precision-operation")
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity)
        for name in ["barrier-tunnelling","global-kinetic-uncertainty","reactive-surrogate-train",
                     "ring-polymer-sample","reactive-surrogate-sample","reactive-surrogate-label",
                     "ring-polymer-convergence","reactive-surrogate-coverage"] {
            _ = try registry.definition("vivo.platform."+name)
        }
        let operation = try registry.definition("vivo.platform.barrier-tunnelling").operation
        let r = VivoBarrierTunnellingRequest(model: .asymmetricEckart,temperatureK: 300,imaginaryWavenumberPerCM: 1000,
            forwardBarrierKJPerMol: 40,reverseBarrierKJPerMol: 60,hamiltonianFingerprint: identity,stationaryPointEvidence: identity,
            energyConvention: "synthetic barrier for interface qualification only")
        let input = ["request":try VivoCanonicalJSON.encode(r)]
        let output = try await operation.execute(.object([:]),input,.init())
        try operation.validateOutputs(.object([:]),input,output,.init())
        let result = try VivoPlatformOperations.input(VivoBarrierTunnellingResult.self,"correction",output)
        #expect(abs(result.logFactor-1.111947262990766) < 1e-7)
        var forged = output; forged["correction"] = Data("null".utf8)
        #expect(throws: (any Error).self) { try operation.validateOutputs(.object([:]),input,forged,.init()) }
    }
}
