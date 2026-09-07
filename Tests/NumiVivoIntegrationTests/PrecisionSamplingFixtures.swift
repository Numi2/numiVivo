import Foundation
import Testing
@testable import NumiVivoKit

enum PrecisionSamplingFixtures {
    static func id(_ text: String) throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(Data(text.utf8)) }
    static func harmonic(k: Double = 200_000, mass: Double = 1, name: String = "harmonic") throws -> VivoNuclearPotential {
        let definition = try VivoNuclearPotentialDefinition(hamiltonianFingerprint: id(name), atomIndices: [0], particleIndices: [0], massesDa: [mass])
        return try .init(definition: definition) { q in
            try .init(definition: definition, positionsNM: q, energyKJPerMol: 0.5*k*q[0].squaredNorm,
                      forcesKJPerMolNM: [q[0] * -k], normalizedConvergenceResidual: 0)
        }
    }
    static func pairPositions(_ r: Double) -> [VivoVector3D] { [.init(-r/2,0,0), .init(r/2,0,0)] }
    static func pairPotential(authority: Bool) throws -> VivoNuclearPotential {
        let d = try VivoNuclearPotentialDefinition(hamiltonianFingerprint: id(authority ? "pair-target" : "pair-baseline"),
            atomIndices: [5,8], particleIndices: [2,0], massesDa: [12,1])
        return try .init(definition: d) { q in
            var energy = q.reduce(0) { $0+0.1*$1.squaredNorm }, f = q.map { $0 * -0.2 }
            if authority {
                let dr = q[1]-q[0], r = dr.norm, x = r-1
                guard r > 1e-10 else { throw VivoChemistryError.invalid("overlapping pair fixture") }
                energy += 20*x*x+2*x*x*x
                let force = dr*((40*x+6*x*x)/r)
                f[0] = f[0]+force; f[1] = f[1]-force
            }
            return try .init(definition: d, positionsNM: q, energyKJPerMol: energy, forcesKJPerMolNM: f, normalizedConvergenceResidual: 0)
        }
    }
    static func trainedPair() async throws -> (VivoReactiveSurrogateModel,VivoNuclearPotential,VivoNuclearPotential,[VivoReactiveTrainingLabel]) {
        let a = try pairPotential(authority: true), b = try pairPotential(authority: false)
        var labels: [VivoReactiveTrainingLabel] = []
        for i in 0...20 {
            let q = pairPositions(0.8+0.02*Double(i))
            labels.append(try await .init(identifier: "train-\(i)", sourceGroup: "train-\(i)", authority: a.checked(q), baseline: b.checked(q)))
        }
        for (i,r) in [0.85,1.05,1.15].enumerated() {
            let q = pairPositions(r)
            labels.append(try await .init(identifier: "held-\(i)", sourceGroup: "held", authority: a.checked(q), baseline: b.checked(q)))
        }
        let cfg = VivoReactiveSurrogateConfiguration(features: [.init(atomA: 0,atomB: 1,scaleNM: 0.2)],
            maximumCenters: 12, committeeSize: 3, ridge: 1e-8,
            energyFitScaleKJPerMol: 0.1, forceFitScaleKJPerMolNM: 1,
            maximumHeldOutEnergyErrorKJPerMol: 0.01, maximumHeldOutForceErrorKJPerMolNM: 0.1,
            maximumEnergyDisagreementKJPerMol: 0.1, maximumForceDisagreementKJPerMolNM: 1,
            maximumNormalizedDistance: 1, seed: 72)
        let model = try VivoReactiveDeltaSurrogate.train(authority: a.definition, baseline: b.definition,
            labels: labels, heldOutGroups: ["held"], configuration: cfg)
        return (model,a,b,labels)
    }
}
