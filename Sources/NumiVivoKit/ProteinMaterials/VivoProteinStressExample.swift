import Foundation

public enum VivoProteinStressExample {
    /// A two-particle harmonic mechanics fixture, NOT a protein, force-field
    /// qualification, or reproduction of the Nature Chemistry experiment.
    public static func harmonicFixture() throws -> VivoProteinStressRequest {
        let identity = try VivoCanonicalJSON.fingerprint(Data("numivivo/protein-stress/harmonic-fixture/v1".utf8))
        let particles = (0..<2).map { VivoClassicalParticle(index: UInt32($0), atomIndex: UInt32($0),
            typeIdentifier: "fixture-C", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) }
        let system = VivoClassicalSystem(identifier: "protein-stress-harmonic-fixture", structureFingerprint: identity,
            particles: particles, bonds: [.init(a: 0, b: 1, lengthNM: 0.25, forceConstant: 100)])
        let configuration = VivoMDConfiguration(timeStepPS: 0.0005, cutoffNM: 1, electrostatics: .cutoff,
            targetTemperatureK: 300, neighborListEnabled: false, randomSeed: 42)
        let checkpoint = VivoMDCheckpoint(systemFingerprint: try system.fingerprint(), configurationFingerprint: try configuration.fingerprint(),
            acceptedStep: 0, timePS: 0, positionsNM: [.zero, .init(0.25, 0, 0)],
            velocitiesNMPerPS: [.zero, .zero], periodicCell: nil)
        return .init(system: system, sourceConfiguration: configuration, sourceCheckpoint: checkpoint,
            sourceDescription: "synthetic two-particle harmonic fixture; not an equilibrated protein", replicaID: "fixture-42", randomSeed: 42,
            pull: .init(reference: .init(particles: [0]), moving: .init(particles: [1]), stiffnessKJPerMolNM2: 50),
            stages: [.init(name: "hold", steps: 4, temperatureK: 300, pullReferenceNM: 0.25),
                     .init(name: "pull", steps: 4, temperatureK: 320, pullReferenceNM: 0.30)],
            sampleEvery: 2, selection: [0, 1], nativeContacts: [.init(a: 0, b: 1, referenceDistanceNM: 0.25)])
    }
}
