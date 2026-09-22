import Foundation

public enum VivoProteinStressError: Error, Sendable, CustomStringConvertible {
    case invalid(String)
    public var description: String {
        switch self { case .invalid(let reason): return "Protein stress test: " + reason }
    }
}

/// Physical-particle indices, not structure-atom indices. The compiler supplies
/// actual particle masses; weights are normalized here and in the force gradient.
public struct VivoProteinPullGroup: Codable, Sendable, Equatable {
    public let particles: [UInt32]
    public init(particles: [UInt32]) { self.particles = particles }
}

public struct VivoProteinPullDefinition: Codable, Sendable, Equatable {
    public let reference: VivoProteinPullGroup
    public let moving: VivoProteinPullGroup
    /// nil is the radial COM distance; otherwise this is a fixed laboratory axis.
    public let projectionAxis: VivoVector3D?
    public let stiffnessKJPerMolNM2: Double
    public init(reference: VivoProteinPullGroup, moving: VivoProteinPullGroup,
                projectionAxis: VivoVector3D? = nil, stiffnessKJPerMolNM2: Double) {
        self.reference = reference; self.moving = moving
        self.projectionAxis = projectionAxis; self.stiffnessKJPerMolNM2 = stiffnessKJPerMolNM2
    }
    public func validate(massesDa: [Double]) throws {
        guard !massesDa.isEmpty, massesDa.count <= Int(UInt32.max),
              massesDa.allSatisfy({ $0.isFinite && $0 >= 0 }),
              stiffnessKJPerMolNM2.isFinite, stiffnessKJPerMolNM2 > 0 else {
            throw VivoProteinStressError.invalid("invalid masses or spring stiffness")
        }
        for group in [reference, moving] {
            guard !group.particles.isEmpty, Set(group.particles).count == group.particles.count,
                  group.particles.allSatisfy({ Int($0) < massesDa.count && massesDa[Int($0)] > 0 }) else {
                throw VivoProteinStressError.invalid("pull groups must contain unique massive physical particles")
            }
            let total = group.particles.reduce(0.0) { $0 + massesDa[Int($1)] }
            guard total.isFinite, total > 0 else { throw VivoProteinStressError.invalid("pull-group mass overflow") }
        }
        guard Set(reference.particles).isDisjoint(with: moving.particles) else {
            throw VivoProteinStressError.invalid("pull groups overlap")
        }
        if let axis = projectionAxis {
            guard axis.isFinite, axis.norm.isFinite, axis.norm > 1e-12 else {
                throw VivoProteinStressError.invalid("projection axis is zero or nonfinite")
            }
        }
    }
}

public struct VivoProteinPullEvaluation: Codable, Sendable, Equatable {
    public let coordinateNM: Double
    public let referenceNM: Double
    public let energyKJPerMol: Double
    /// Positive means the restraint pulls the moving group toward larger q.
    public let tensileForcePN: Double
    public let forcesKJPerMolNM: [VivoVector3D]
}

public enum VivoProteinPullMath {
    /// Exact SI Avogadro constant. 1 kJ mol^-1 nm^-1 = 1.660539... pN.
    public static let forceToPN = 1e24 / 6.02214076e23

    /// Coordinates must already describe one whole molecule. Never independently
    /// minimum-image the two COMs: doing so changes the measured extension.
    public static func evaluate(_ definition: VivoProteinPullDefinition,
                                wholePositionsNM: [VivoVector3D], massesDa: [Double],
                                referenceNM: Double) throws -> VivoProteinPullEvaluation {
        try definition.validate(massesDa: massesDa)
        guard wholePositionsNM.count == massesDa.count, wholePositionsNM.allSatisfy(\.isFinite),
              referenceNM.isFinite else { throw VivoProteinStressError.invalid("pull geometry or reference") }
        func center(_ group: VivoProteinPullGroup) -> VivoVector3D {
            let total = group.particles.reduce(0.0) { $0 + massesDa[Int($1)] }
            // Normalize first to avoid multiplying large coordinates by masses.
            return group.particles.reduce(.zero) { $0 + wholePositionsNM[Int($1)] * (massesDa[Int($1)] / total) }
        }
        let delta = center(definition.moving) - center(definition.reference)
        let direction: VivoVector3D, q: Double
        if let axis = definition.projectionAxis {
            direction = axis / axis.norm; q = delta.dot(direction)
        } else {
            guard delta.norm.isFinite, delta.norm > 1e-12, referenceNM >= 0 else {
                throw VivoProteinStressError.invalid("radial pulling requires distinct COMs and a nonnegative reference")
            }
            direction = delta / delta.norm; q = delta.norm
        }
        let displacement = q - referenceNM
        let force = -definition.stiffnessKJPerMolNM2 * displacement
        let energy = 0.5 * definition.stiffnessKJPerMolNM2 * displacement * displacement
        var forces = [VivoVector3D](repeating: .zero, count: massesDa.count)
        for (group, sign) in [(definition.reference, -1.0), (definition.moving, 1.0)] {
            let total = group.particles.reduce(0.0) { $0 + massesDa[Int($1)] }
            for index in group.particles {
                forces[Int(index)] = direction * (sign * force * massesDa[Int(index)] / total)
            }
        }
        guard q.isFinite, energy.isFinite, (force * forceToPN).isFinite,
              forces.allSatisfy(\.isFinite) else { throw VivoProteinStressError.invalid("pull evaluation overflow") }
        return .init(coordinateNM: q, referenceNM: referenceNM, energyKJPerMol: energy,
                     tensileForcePN: force * forceToPN, forcesKJPerMolNM: forces)
    }

    /// Exact protocol work for an instantaneous parameter switch at fixed x.
    /// This is not displacement work, thermostat heat, or a folding free energy.
    public static func switchingWork(coordinateNM: Double, from: Double, to: Double,
                                     stiffnessKJPerMolNM2: Double) throws -> Double {
        guard [coordinateNM, from, to, stiffnessKJPerMolNM2].allSatisfy(\.isFinite),
              stiffnessKJPerMolNM2 > 0 else { throw VivoProteinStressError.invalid("switching-work inputs") }
        let result = stiffnessKJPerMolNM2 * ((from / 2 + to / 2) - coordinateNM) * (to - from)
        guard result.isFinite else { throw VivoProteinStressError.invalid("switching-work overflow") }
        return result
    }
}

/// Explicit D-H-A identity. A donor's covalent hydrogen is not itself the
/// noncovalent hydrogen bond being measured. No donor/acceptor chemistry is guessed.
public struct VivoProteinHydrogenBond: Codable, Sendable, Equatable {
    public let donor: UInt32
    public let hydrogen: UInt32
    public let acceptor: UInt32
    public init(donor: UInt32, hydrogen: UInt32, acceptor: UInt32) {
        self.donor = donor; self.hydrogen = hydrogen; self.acceptor = acceptor
    }
}

public struct VivoProteinHydrogenBondCriteria: Codable, Sendable, Equatable {
    public let maximumDonorAcceptorNM: Double
    public let maximumHydrogenAcceptorNM: Double
    /// D-H-A angle; a straight bond is 180 degrees, not zero degrees.
    public let minimumDHAAngleDegrees: Double
    public init(maximumDonorAcceptorNM: Double = 0.35, maximumHydrogenAcceptorNM: Double = 0.25,
                minimumDHAAngleDegrees: Double = 150) {
        self.maximumDonorAcceptorNM = maximumDonorAcceptorNM
        self.maximumHydrogenAcceptorNM = maximumHydrogenAcceptorNM
        self.minimumDHAAngleDegrees = minimumDHAAngleDegrees
    }
    public func validate() throws {
        guard maximumDonorAcceptorNM.isFinite, maximumDonorAcceptorNM > 0,
              maximumHydrogenAcceptorNM.isFinite, maximumHydrogenAcceptorNM > 0,
              minimumDHAAngleDegrees.isFinite, minimumDHAAngleDegrees > 0,
              minimumDHAAngleDegrees <= 180 else {
            throw VivoProteinStressError.invalid("hydrogen-bond geometric criteria")
        }
    }
}

public struct VivoProteinHydrogenBondGeometry: Codable, Sendable, Equatable {
    public let present: Bool
    public let donorAcceptorNM: Double
    public let hydrogenAcceptorNM: Double
    public let dhaAngleDegrees: Double
}

public enum VivoProteinHydrogenBondAnalysis {
    public static func evaluate(_ bonds: [VivoProteinHydrogenBond], positionsNM: [VivoVector3D],
                                criteria: VivoProteinHydrogenBondCriteria = .init()) throws
    -> [VivoProteinHydrogenBondGeometry] {
        try criteria.validate()
        guard bonds.count <= 1_000_000, positionsNM.allSatisfy(\.isFinite) else {
            throw VivoProteinStressError.invalid("hydrogen-bond analysis size or coordinates")
        }
        var seen = Set<[UInt32]>()
        return try bonds.map { bond in
            let indices = [bond.donor, bond.hydrogen, bond.acceptor]
            guard Set(indices).count == 3, indices.allSatisfy({ Int($0) < positionsNM.count }),
                  seen.insert(indices).inserted else {
                throw VivoProteinStressError.invalid("duplicate or invalid D-H-A identity")
            }
            let d = positionsNM[Int(bond.donor)], h = positionsNM[Int(bond.hydrogen)]
            let a = positionsNM[Int(bond.acceptor)], dh = d - h, ah = a - h
            let da = (d - a).norm, ha = ah.norm, donorHydrogen = dh.norm
            guard [da, ha, donorHydrogen].allSatisfy(\.isFinite), ha > 1e-12, donorHydrogen > 1e-12 else {
                throw VivoProteinStressError.invalid("degenerate D-H-A geometry")
            }
            let cosine = max(-1, min(1, (dh / donorHydrogen).dot(ah / ha)))
            let angle = acos(cosine) * 180 / Double.pi
            return .init(present: da <= criteria.maximumDonorAcceptorNM && ha <= criteria.maximumHydrogenAcceptorNM
                            && angle >= criteria.minimumDHAAngleDegrees,
                         donorAcceptorNM: da, hydrogenAcceptorNM: ha, dhaAngleDegrees: angle)
        }
    }

    /// Occupancy is time-weighted with a left-constant interpolation convention.
    /// First-loss times are observation times, not resolved bond lifetimes.
    public static func persistence(timesPS: [Double], present: [[Bool]]) throws -> [VivoProteinHydrogenBondPersistence] {
        guard timesPS.count >= 2, present.count == timesPS.count,
              timesPS.allSatisfy({ $0.isFinite && $0 >= 0 }),
              zip(timesPS, timesPS.dropFirst()).allSatisfy({ $0 < $1 }),
              let first = present.first, present.allSatisfy({ $0.count == first.count }) else {
            throw VivoProteinStressError.invalid("persistence requires aligned observations at increasing times")
        }
        let duration = timesPS.last! - timesPS[0]
        return first.indices.map { bond in
            var occupied = 0.0, firstLoss: Double?, reformations = 0
            for i in 0..<(timesPS.count - 1) {
                if present[i][bond] { occupied += timesPS[i + 1] - timesPS[i] }
                if present[i][bond] && !present[i + 1][bond] && firstLoss == nil { firstLoss = timesPS[i + 1] }
                if !present[i][bond] && present[i + 1][bond] { reformations += 1 }
            }
            return .init(initiallyPresent: first[bond], timeWeightedOccupancy: occupied / duration,
                         firstObservedLossTimePS: firstLoss, observedReformations: reformations)
        }
    }
}

public struct VivoProteinHydrogenBondPersistence: Codable, Sendable, Equatable {
    public let initiallyPresent: Bool
    public let timeWeightedOccupancy: Double
    public let firstObservedLossTimePS: Double?
    public let observedReformations: Int
}

/// Caller-selected reference contacts. These need not be hydrogen bonds.
public struct VivoProteinNativeContact: Codable, Sendable, Equatable {
    public let a: UInt32
    public let b: UInt32
    public let referenceDistanceNM: Double
    public init(a: UInt32, b: UInt32, referenceDistanceNM: Double) {
        self.a = a; self.b = b; self.referenceDistanceNM = referenceDistanceNM
    }
}

public struct VivoProteinStructuralRetention: Codable, Sendable, Equatable {
    public let retainedContactFraction: Double?
    public let contactDistanceRMSErrorNM: Double?
    public let radiusOfGyrationNM: Double
}

public enum VivoProteinStructuralAnalysis {
    public static func evaluate(positionsNM: [VivoVector3D], massesDa: [Double],
                                selection: [UInt32], contacts: [VivoProteinNativeContact],
                                maximumDistanceRatio: Double = 1.2) throws -> VivoProteinStructuralRetention {
        guard positionsNM.count == massesDa.count, positionsNM.allSatisfy(\.isFinite),
              !selection.isEmpty, Set(selection).count == selection.count,
              selection.allSatisfy({ Int($0) < massesDa.count && massesDa[Int($0)].isFinite && massesDa[Int($0)] > 0 }),
              contacts.count <= 1_000_000, maximumDistanceRatio.isFinite, maximumDistanceRatio >= 1 else {
            throw VivoProteinStressError.invalid("structural-retention selection or criteria")
        }
        let members = Set(selection)
        let totalMass = selection.reduce(0.0) { $0 + massesDa[Int($1)] }
        guard totalMass.isFinite else { throw VivoProteinStressError.invalid("selection mass overflow") }
        let center = selection.reduce(VivoVector3D.zero) { $0 + positionsNM[Int($1)] * (massesDa[Int($1)] / totalMass) }
        let rg2 = selection.reduce(0.0) { $0 + (positionsNM[Int($1)] - center).squaredNorm * (massesDa[Int($1)] / totalMass) }
        var retained = 0, error2 = 0.0, seen = Set<UInt64>()
        for contact in contacts {
            let key = UInt64(min(contact.a, contact.b)) << 32 | UInt64(max(contact.a, contact.b))
            guard contact.a != contact.b, members.contains(contact.a), members.contains(contact.b),
                  contact.referenceDistanceNM.isFinite, contact.referenceDistanceNM > 0,
                  seen.insert(key).inserted else { throw VivoProteinStressError.invalid("invalid or duplicate reference contact") }
            let distance = (positionsNM[Int(contact.a)] - positionsNM[Int(contact.b)]).norm
            if distance / contact.referenceDistanceNM <= maximumDistanceRatio { retained += 1 }
            error2 += pow(distance - contact.referenceDistanceNM, 2)
        }
        guard rg2.isFinite, error2.isFinite else { throw VivoProteinStressError.invalid("structural metric overflow") }
        return .init(retainedContactFraction: contacts.isEmpty ? nil : Double(retained) / Double(contacts.count),
                     contactDistanceRMSErrorNM: contacts.isEmpty ? nil : sqrt(error2 / Double(contacts.count)),
                     radiusOfGyrationNM: sqrt(rg2))
    }
}

public struct VivoProteinReplicaStatistics: Codable, Sendable, Equatable {
    public let count: Int
    public let mean: Double
    public let sampleStandardDeviation: Double?
    public let standardErrorOfMean: Double?
    public let minimum: Double
    public let maximum: Double
    /// Descriptive between-replica spread; not a calibrated physical uncertainty.
    public static func calculate(_ values: [Double]) throws -> Self {
        guard !values.isEmpty, values.allSatisfy(\.isFinite) else {
            throw VivoProteinStressError.invalid("replica statistics require finite observations")
        }
        var mean = 0.0, m2 = 0.0
        for (i, x) in values.enumerated() { let d = x - mean; mean += d / Double(i + 1); m2 += d * (x - mean) }
        guard mean.isFinite, m2.isFinite else { throw VivoProteinStressError.invalid("replica statistics overflow") }
        let sd = values.count > 1 ? sqrt(max(0, m2) / Double(values.count - 1)) : nil
        return .init(count: values.count, mean: mean, sampleStandardDeviation: sd,
                     standardErrorOfMean: sd.map { $0 / sqrt(Double(values.count)) },
                     minimum: values.min()!, maximum: values.max()!)
    }
}
