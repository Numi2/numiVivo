import Foundation

/// Shared hot-path MD command. For electrostatics mode 2 (PME),
/// `reactionFieldK` carries the Ewald splitting parameter beta in nm^-1 and
/// `reactionFieldC` is zero. This preserves the 176-byte ABI while mode-specific
/// kernels interpret the fields explicitly.
struct VivoMDMetalCommand {
    var particleCount: UInt32
    var typeCount: UInt32
    var electrostatics: UInt32
    var periodic: UInt32
    var stepLow: UInt32
    var stepHigh: UInt32
    var seedLow: UInt32
    var seedHigh: UInt32
    var dtPS: Float
    var cutoffNM: Float
    var coulombPrefactor: Float
    var reactionFieldK: Float
    var reactionFieldC: Float
    var minimumDistanceNM: Float
    var constraintTolerance: Float
    var langevinA: Float
    var targetTemperatureK: Float
    var boltzmannKJPerMolK: Float
    var neighborCapacity: UInt32
    var neighborRadiusNM: Float
    var cellA: SIMD4<Float>
    var cellB: SIMD4<Float>
    var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

struct VivoMDMetalStatus {
    var flags: UInt32
    var firstParticle: UInt32
    var violationCount: UInt32
    var reserved: UInt32
}

enum VivoMDMetalABI {
    static func validate() throws {
        let expected: [(String, Int, Int)] = [
            ("Command", MemoryLayout<VivoMDMetalCommand>.stride, 176),
            ("Status", MemoryLayout<VivoMDMetalStatus>.stride, 16),
            ("Bond", MemoryLayout<VivoMDPackedSystem.Bond>.stride, 16),
            ("Angle", MemoryLayout<VivoMDPackedSystem.Angle>.stride, 32),
            ("Torsion", MemoryLayout<VivoMDPackedSystem.Torsion>.stride, 32),
            ("Constraint", MemoryLayout<VivoMDPackedSystem.Constraint>.stride, 16),
            ("Incidence", MemoryLayout<VivoMDPackedSystem.Incidence>.stride, 8),
            ("PairException", MemoryLayout<VivoMDPackedSystem.PairException>.stride, 32),
            ("LinearVirtualSite", MemoryLayout<VivoMDPackedSystem.LinearVirtualSite>.stride, 48),
            ("VirtualParentIncidence", MemoryLayout<VivoMDPackedSystem.VirtualParentIncidence>.stride, 8)
        ]
        for (name, actual, required) in expected where actual != required {
            throw VivoMDRuntimeError.metal("Swift/Metal ABI mismatch for \(name): \(actual) != \(required)")
        }
    }

    static func command(packed: VivoMDPackedSystem,
                        configuration: VivoMDConfiguration,
                        cell: VivoPeriodicCell?,
                        stepIndex: UInt64) throws -> VivoMDMetalCommand {
        let mode: UInt32
        switch configuration.electrostatics {
        case .cutoff: mode = 0
        case .reactionField: mode = 1
        case .pme: mode = 2
        }
        let eps = configuration.relativeDielectric
        let cutoff = configuration.cutoffNM
        var modeK: Float = 0
        var modeC: Float = 0
        if configuration.electrostatics == .reactionField {
            let external = configuration.reactionFieldDielectric
            modeK = Float((external - eps) / (2 * external + eps) / pow(cutoff, 3))
            modeC = Float(3 * external / (2 * external + eps) / cutoff)
        } else if configuration.electrostatics == .pme {
            guard let cell else {
                throw VivoMDRuntimeError.unsupported(["PME command construction requires a periodic cell"])
            }
            let plan = try VivoPMEPlan.make(cell: cell,
                                            cutoffNM: cutoff,
                                            tolerance: configuration.resolvedPMETolerance,
                                            targetGridSpacingNM: configuration.resolvedPMEGridSpacingNM)
            guard plan.ewaldBetaPerNM <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoMDRuntimeError.metal("PME Ewald beta exceeds FP32 command range")
            }
            modeK = Float(plan.ewaldBetaPerNM)
        }

        let langevinA: Float = configuration.thermostat == .langevinMiddle
            ? Float(exp(-(configuration.frictionPerPS ?? 0) * configuration.timeStepPS)) : 1
        let possible = packed.particleCount > 0 ? packed.particleCount - 1 : 0
        let capacity = max(1, min(configuration.resolvedMaximumNeighborsPerParticle,
                                  max(possible, 1)))
        let radius = cutoff + configuration.neighborSkinNM
        guard radius.isFinite, radius > cutoff,
              radius <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoMDRuntimeError.metal("neighbor-list radius is invalid")
        }
        let reciprocal = try cell.map(reciprocalCell)
        func f4(_ value: VivoVector3D?) -> SIMD4<Float> {
            guard let value else { return .zero }
            return .init(Float(value.x), Float(value.y), Float(value.z), 0)
        }
        return .init(
            particleCount: packed.particleCount, typeCount: packed.typeCount,
            electrostatics: mode, periodic: cell == nil ? 0 : 1,
            stepLow: UInt32(truncatingIfNeeded: stepIndex),
            stepHigh: UInt32(truncatingIfNeeded: stepIndex >> 32),
            seedLow: UInt32(truncatingIfNeeded: configuration.randomSeed),
            seedHigh: UInt32(truncatingIfNeeded: configuration.randomSeed >> 32),
            dtPS: Float(configuration.timeStepPS), cutoffNM: Float(cutoff),
            coulombPrefactor: Float(138.935456 / eps),
            reactionFieldK: modeK, reactionFieldC: modeC,
            minimumDistanceNM: 1e-5,
            constraintTolerance: Float(configuration.constraintTolerance),
            langevinA: langevinA,
            targetTemperatureK: Float(configuration.targetTemperatureK ?? 0),
            boltzmannKJPerMolK: 0.00831446261815324,
            neighborCapacity: capacity, neighborRadiusNM: Float(radius),
            cellA: f4(cell?.a), cellB: f4(cell?.b), cellC: f4(cell?.c),
            reciprocalA: f4(reciprocal?.0), reciprocalB: f4(reciprocal?.1),
            reciprocalC: f4(reciprocal?.2)
        )
    }

    static func validateMinimumImageCutoff(_ cutoff: Double, cell: VivoPeriodicCell) throws {
        try validateRadius(cutoff, cell: cell, label: "cutoff")
    }
    static func validateNeighborRadius(_ radius: Double, cell: VivoPeriodicCell) throws {
        try validateRadius(radius, cell: cell, label: "cutoff + neighbor skin")
    }
    private static func validateRadius(_ radius: Double,
                                       cell: VivoPeriodicCell,
                                       label: String) throws {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        let heights = [abs(determinant) / cell.b.cross(cell.c).norm,
                       abs(determinant) / cell.c.cross(cell.a).norm,
                       abs(determinant) / cell.a.cross(cell.b).norm]
        guard heights.allSatisfy({ $0.isFinite && $0 > 0 }),
              radius < 0.5 * (heights.min() ?? 0) else {
            throw VivoMDRuntimeError.unsupported([
                "\(label) must be less than half the shortest periodic face height"
            ])
        }
    }
    private static func reciprocalCell(_ cell: VivoPeriodicCell) throws
    -> (VivoVector3D, VivoVector3D, VivoVector3D) {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else {
            throw VivoMDRuntimeError.metal("singular periodic cell")
        }
        return (cell.b.cross(cell.c) / determinant,
                cell.c.cross(cell.a) / determinant,
                cell.a.cross(cell.b) / determinant)
    }
}
