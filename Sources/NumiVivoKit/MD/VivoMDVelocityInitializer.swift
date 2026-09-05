import Foundation

public struct VivoMDVelocityInitializationResult: Codable, Sendable, Equatable {
    public let recipe: VivoMDVelocityInitialization
    public let velocitiesNMPerPS: [VivoVector3D]
    public let projectionIterations: UInt32
    public let maximumNormalizedConstraintResidual: Double
    public let removesCenterOfMassMomentum: Bool
}

/// One-time preparation, not a thermostat. Samples a Gaussian in mass-weighted
/// coordinates, then projects onto the tangent space of the fixed-position
/// distance constraints. It does not rescale kinetic energy to force an exact
/// temperature, which would change the Maxwell kinetic-energy distribution.
public enum VivoMDVelocityInitializer {
    public static func make(system: VivoClassicalSystem,
                            positionsNM: [VivoVector3D], periodicCell: VivoPeriodicCell?,
                            recipe: VivoMDVelocityInitialization,
                            tolerance: Double = 1e-8,
                            maximumProjectionIterations: UInt32 = 2048) throws -> VivoMDVelocityInitializationResult {
        try VivoClassicalSystemValidator.validate(system)
        try recipe.validate()
        guard positionsNM.count == system.particles.count,
              positionsNM.allSatisfy(\.isFinite), periodicCell?.isValid != false,
              tolerance.isFinite, tolerance > 0, maximumProjectionIterations > 0 else {
            throw VivoArtifactValidationError.invalid("invalid thermal velocity preparation inputs")
        }
        var random = GaussianStream(seed: recipe.seed)
        let thermalEnergy = 0.00831446261815324 * recipe.temperatureK
        guard thermalEnergy.isFinite, thermalEnergy > 0 else {
            throw VivoArtifactValidationError.invalid("thermal energy is not representable")
        }
        var velocities = [VivoVector3D](repeating: .zero, count: system.particles.count)
        var inverseMass = [Double](repeating: 0, count: system.particles.count)
        for particle in system.particles {
            guard particle.role == .atom || particle.role == .virtualSite else {
                throw VivoArtifactValidationError.incompatible("thermal initialization does not support Drude dynamics")
            }
            if particle.massDa == 0 { continue }
            let index = Int(particle.index)
            inverseMass[index] = 1 / particle.massDa
            let sigma = sqrt(thermalEnergy * inverseMass[index])
            guard sigma.isFinite else { throw VivoArtifactValidationError.invalid("thermal velocity variance overflow") }
            velocities[index] = .init(sigma * random.normal(), sigma * random.normal(), sigma * random.normal())
        }
        if recipe.removeCenterOfMass {
            let mass = system.particles.reduce(0.0) { $0 + $1.massDa }
            guard mass.isFinite, mass > 0 else { throw VivoArtifactValidationError.invalid("total physical mass is invalid") }
            var momentum = VivoVector3D.zero
            for particle in system.particles where particle.massDa > 0 {
                momentum = momentum + velocities[Int(particle.index)] * particle.massDa
            }
            let drift = momentum / mass
            for particle in system.particles where particle.massDa > 0 {
                velocities[Int(particle.index)] = velocities[Int(particle.index)] - drift
            }
        }

        struct ConstraintRow {
            let a: Int
            let b: Int
            let direction: VivoVector3D
            let denominator: Double
            let residualScale: Double
        }
        var rows: [ConstraintRow] = []
        rows.reserveCapacity(system.constraints.count)
        for constraint in system.constraints {
            let a = Int(constraint.a), b = Int(constraint.b)
            guard inverseMass[a] > 0, inverseMass[b] > 0 else {
                throw VivoArtifactValidationError.incompatible("distance constraints must join massive physical particles")
            }
            let d = try VivoMDPreparationGeometry.minimumImage(positionsNM[a] - positionsNM[b], cell: periodicCell)
            let r2 = d.squaredNorm
            let sumInverseMass = inverseMass[a] + inverseMass[b]
            let denominator = r2 * sumInverseMass
            let residualScale = sqrt(r2 * thermalEnergy * sumInverseMass)
            guard r2.isFinite, r2 > 1e-20, denominator.isFinite, denominator > 0,
                  residualScale.isFinite, residualScale > 0 else {
                throw VivoArtifactValidationError.invalid("degenerate constraint in thermal velocity projection")
            }
            rows.append(.init(a: a, b: b, direction: d, denominator: denominator, residualScale: residualScale))
        }
        var residual = Double.infinity
        var iterations: UInt32 = 0
        if rows.isEmpty { residual = 0 }
        while !rows.isEmpty, iterations < maximumProjectionIterations {
            for row in rows {
                let multiplier = row.direction.dot(velocities[row.a] - velocities[row.b]) / row.denominator
                velocities[row.a] = velocities[row.a] - row.direction * (inverseMass[row.a] * multiplier)
                velocities[row.b] = velocities[row.b] + row.direction * (inverseMass[row.b] * multiplier)
            }
            iterations += 1
            residual = 0
            for row in rows {
                residual = max(residual, abs(row.direction.dot(velocities[row.a] - velocities[row.b])) / row.residualScale)
            }
            if residual <= tolerance { break }
        }
        guard residual.isFinite, residual <= tolerance else {
            throw VivoArtifactValidationError.invalid("thermal velocity constraint projection exhausted its iteration budget")
        }
        let sites = system.linearVirtualSites ?? []
        let expectedSites = Set(system.particles.filter { $0.role == .virtualSite }.map(\.index))
        guard Set(sites.map(\.siteParticle)) == expectedSites else {
            throw VivoArtifactValidationError.unresolved("thermal initialization requires resolved virtual sites")
        }
        for site in sites {
            var velocity = VivoVector3D.zero
            for (parent, weight) in zip(site.parentParticles, site.weights) {
                velocity = velocity + velocities[Int(parent)] * weight
            }
            velocities[Int(site.siteParticle)] = velocity
        }
        guard velocities.allSatisfy(\.isFinite) else {
            throw VivoArtifactValidationError.invalid("initialized velocities contain nonfinite values")
        }
        return .init(recipe: recipe, velocitiesNMPerPS: velocities,
                     projectionIterations: iterations,
                     maximumNormalizedConstraintResidual: residual,
                     removesCenterOfMassMomentum: recipe.removeCenterOfMass)
    }

    private struct GaussianStream {
        var state: UInt64
        var spare: Double?
        init(seed: UInt64) { state = seed ^ 0x56454c4f43495459 }
        mutating func uniform() -> Double {
            state &+= 0x9e3779b97f4a7c15
            var value = state
            value = (value ^ (value >> 30)) &* 0xbf58476d1ce4e5b9
            value = (value ^ (value >> 27)) &* 0x94d049bb133111eb
            value ^= value >> 31
            // A 52-bit midpoint grid avoids rounding the largest draw to 1.
            return (Double(value >> 12) + 0.5) * 0x1p-52
        }
        mutating func normal() -> Double {
            if let value = spare { spare = nil; return value }
            let radius = sqrt(-2 * log(uniform()))
            let angle = 2 * Double.pi * uniform()
            spare = radius * sin(angle)
            return radius * cos(angle)
        }
    }
}

/// Exact bounded closest-image search for preparation-time geometry. Fractional
/// rounding alone is not a nearest-image algorithm for arbitrary skew cells.
/// The reciprocal-vector bounds enumerate every image that could improve the
/// initial candidate. Highly ill-conditioned cells fail rather than guessing.
enum VivoMDPreparationGeometry {
    static func minimumImage(_ displacement: VivoVector3D,
                             cell: VivoPeriodicCell?) throws -> VivoVector3D {
        guard displacement.isFinite else { throw VivoArtifactValidationError.invalid("nonfinite displacement") }
        guard let cell else { return displacement }
        guard cell.isValid else { throw VivoArtifactValidationError.invalid("invalid periodic cell") }
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        let reciprocals = [cell.b.cross(cell.c) / determinant,
                           cell.c.cross(cell.a) / determinant,
                           cell.a.cross(cell.b) / determinant]
        let fractions = reciprocals.map { $0.dot(displacement) }
        let reduced = displacement - cell.a * fractions[0].rounded()
            - cell.b * fractions[1].rounded() - cell.c * fractions[2].rounded()
        var best = reduced
        var bestSquared = best.squaredNorm
        guard bestSquared.isFinite else { throw VivoArtifactValidationError.invalid("periodic displacement overflow") }
        let radius = sqrt(bestSquared)
        var intervals: [ClosedRange<Int>] = []
        var candidateCount = 1
        for reciprocal in reciprocals {
            let center = reciprocal.dot(reduced)
            let bound = reciprocal.norm * radius + 32 * Double.ulpOfOne
            let low = ceil(center - bound), high = floor(center + bound)
            guard low.isFinite, high.isFinite, low >= -4096, high <= 4096, low <= high else {
                throw VivoArtifactValidationError.incompatible("periodic cell is too ill-conditioned for bounded image selection")
            }
            let lower = Int(low), upper = Int(high)
            let product = candidateCount.multipliedReportingOverflow(by: upper - lower + 1)
            guard !product.overflow, product.partialValue <= 4096 else {
                throw VivoArtifactValidationError.incompatible("periodic image search exceeds 4096 candidates")
            }
            candidateCount = product.partialValue
            intervals.append(lower...upper)
        }
        for i in intervals[0] {
            for j in intervals[1] {
                for k in intervals[2] {
                    let candidate = reduced - cell.a * Double(i) - cell.b * Double(j) - cell.c * Double(k)
                    let squared = candidate.squaredNorm
                    if squared < bestSquared { best = candidate; bestSquared = squared }
                }
            }
        }
        return best
    }
}
