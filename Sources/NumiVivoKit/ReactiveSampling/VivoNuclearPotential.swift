import Foundation

/// One independent physical-atom coordinate layout shared by nuclear sampling
/// and reactive surrogates. Virtual sites remain owned by the force authority.
public struct VivoNuclearPotentialDefinition: Codable, Sendable, Equatable {
    public enum CoordinateEvaluation: String, Codable, Sendable {
        case continuousFP64
        /// Explicit numerical potential U_FP32(round(q)); NEVER presented as a
        /// continuous FP64 energy. Outer Monte Carlo uses this same definition.
        case projectedFP32
    }
    public let hamiltonianFingerprint: VivoFingerprint
    public let atomIndices: [UInt32]
    public let particleIndices: [UInt32]
    public let massesDa: [Double]
    public let periodicCell: VivoPeriodicCell?
    /// Independent-atom slots translated together when choosing a periodic chart.
    public let periodicMoleculeGroups: [[Int]]
    public let coordinateEvaluation: CoordinateEvaluation
    public init(hamiltonianFingerprint: VivoFingerprint, atomIndices: [UInt32], particleIndices: [UInt32],
                massesDa: [Double], periodicCell: VivoPeriodicCell? = nil,
                periodicMoleculeGroups: [[Int]] = [], coordinateEvaluation: CoordinateEvaluation = .continuousFP64) throws {
        self.hamiltonianFingerprint = hamiltonianFingerprint; self.atomIndices = atomIndices
        self.particleIndices = particleIndices; self.massesDa = massesDa; self.periodicCell = periodicCell
        self.coordinateEvaluation = coordinateEvaluation; self.periodicMoleculeGroups = periodicMoleculeGroups
        try validate()
    }
    public func validate() throws {
        let n = atomIndices.count
        guard (1...100_000).contains(n), particleIndices.count == n, massesDa.count == n,
              Set(atomIndices).count == n, Set(particleIndices).count == n,
              massesDa.allSatisfy({ $0.isFinite && $0 > 0 && $0 <= 1e6 }),
              periodicCell?.isValid != false,
              periodicCell == nil ? periodicMoleculeGroups.isEmpty : (periodicMoleculeGroups.allSatisfy({ !$0.isEmpty })
                && periodicMoleculeGroups.flatMap({ $0 }).sorted() == Array(0..<n)) else {
            throw VivoChemistryError.invalid("nuclear potential atom/mass/mapping contract")
        }
    }
    /// Translate whole molecules, and ALL beads of each molecule, by one lattice
    /// vector. Internal geometry and ring springs survive; no atomwise wrapping.
    /// This represents zero-winding paths on a periodic molecular coordinate chart.
    public func canonicalBeads(_ beads: [[VivoVector3D]]) throws -> [[VivoVector3D]] {
        try validate()
        guard !beads.isEmpty, beads.allSatisfy({ $0.count == massesDa.count && $0.allSatisfy(\.isFinite) }) else {
            throw VivoChemistryError.invalid("periodic nuclear chart shape")
        }
        guard let cell = periodicCell else { return beads }
        let det = cell.a.dot(cell.b.cross(cell.c))
        var result = beads
        for group in periodicMoleculeGroups {
            let anchor = beads.reduce(VivoVector3D.zero) { $0+$1[group[0]] }/Double(beads.count)
            let fractional = [cell.b.cross(cell.c).dot(anchor)/det,cell.c.cross(cell.a).dot(anchor)/det,cell.a.cross(cell.b).dot(anchor)/det]
            guard fractional.allSatisfy({ $0.isFinite && abs($0) < 0x1p40 }) else { throw VivoChemistryError.resourceLimit("periodic nuclear image range") }
            let shift = cell.a*floor(fractional[0])+cell.b*floor(fractional[1])+cell.c*floor(fractional[2])
            for b in beads.indices { for i in group { result[b][i] = beads[b][i]-shift } }
        }
        return result
    }
    public func canonicalPositions(_ positions: [VivoVector3D]) throws -> [VivoVector3D] { try canonicalBeads([positions])[0] }
    public func fingerprint() throws -> VivoFingerprint {
        try validate(); return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}
public struct VivoNuclearPotentialEvaluation: Codable, Sendable, Equatable {
    public let definitionFingerprint: VivoFingerprint
    public let requestedPositionsNM: [VivoVector3D]
    public let energyKJPerMol: Double
    public let forcesKJPerMolNM: [VivoVector3D]
    public let normalizedConvergenceResidual: Double
    public init(definition: VivoNuclearPotentialDefinition, positionsNM: [VivoVector3D],
                energyKJPerMol: Double, forcesKJPerMolNM: [VivoVector3D],
                normalizedConvergenceResidual: Double) throws {
        definitionFingerprint = try definition.fingerprint(); requestedPositionsNM = positionsNM
        self.energyKJPerMol = energyKJPerMol; self.forcesKJPerMolNM = forcesKJPerMolNM
        self.normalizedConvergenceResidual = normalizedConvergenceResidual
        try validate(definition: definition, positionsNM: positionsNM)
    }
    public func validate(definition: VivoNuclearPotentialDefinition, positionsNM: [VivoVector3D]) throws {
        guard definitionFingerprint == (try definition.fingerprint()), requestedPositionsNM == positionsNM,
              positionsNM.count == definition.atomIndices.count, positionsNM.allSatisfy(\.isFinite),
              forcesKJPerMolNM.count == positionsNM.count, forcesKJPerMolNM.allSatisfy(\.isFinite), energyKJPerMol.isFinite,
              normalizedConvergenceResidual.isFinite, normalizedConvergenceResidual >= 0,
              normalizedConvergenceResidual <= 1 else {
            throw VivoChemistryError.convergence("incomplete, stale or unconverged nuclear potential evaluation")
        }
    }
}
public struct VivoNuclearPotential: Sendable {
    public let definition: VivoNuclearPotentialDefinition
    public let evaluate: @Sendable ([VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation
    public init(definition: VivoNuclearPotentialDefinition,
                evaluate: @escaping @Sendable ([VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation) throws {
        try definition.validate(); self.definition = definition; self.evaluate = evaluate
    }
    public func checked(_ positions: [VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation {
        guard positions.count == definition.atomIndices.count, positions.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("nuclear sampler geometry")
        }
        try Task.checkCancellation()
        let result = try await evaluate(positions)
        try result.validate(definition: definition, positionsNM: positions)
        return result
    }
}
