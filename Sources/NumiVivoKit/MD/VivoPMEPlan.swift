import Foundation

/// Immutable reciprocal-space plan for particle-mesh Ewald electrostatics.
/// Mesh dimensions are powers of two so the Metal backend can use an in-place
/// radix-2 Stockham/Cooley-Tukey schedule without depending on platform FFT APIs.
public struct VivoPMEPlan: Codable, Sendable, Equatable {
    public var gridX: UInt32
    public var gridY: UInt32
    public var gridZ: UInt32
    public var gridPointCount: UInt64
    public var interpolationOrder: UInt32
    public var ewaldBetaPerNM: Double
    public var tolerance: Double
    public var targetGridSpacingNM: Double
    public var cellVolumeNM3: Double
    public var reciprocalA: VivoVector3D
    public var reciprocalB: VivoVector3D
    public var reciprocalC: VivoVector3D

    public static func make(cell: VivoPeriodicCell,
                            cutoffNM: Double,
                            tolerance: Double,
                            targetGridSpacingNM: Double,
                            interpolationOrder: UInt32 = 4) throws -> Self {
        guard cell.isValid, cutoffNM.isFinite, cutoffNM > 0,
              tolerance.isFinite, tolerance > 0, tolerance < 0.1,
              targetGridSpacingNM.isFinite, targetGridSpacingNM > 0,
              interpolationOrder == 4 else {
            throw VivoMDRuntimeError.metal("invalid PME planning request")
        }
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else {
            throw VivoMDRuntimeError.metal("PME periodic cell is singular")
        }
        let reciprocalA = cell.b.cross(cell.c) / determinant
        let reciprocalB = cell.c.cross(cell.a) / determinant
        let reciprocalC = cell.a.cross(cell.b) / determinant
        let beta = try solveBeta(cutoffNM: cutoffNM, tolerance: tolerance)

        let gx = try meshDimension(lengthNM: cell.a.norm, spacingNM: targetGridSpacingNM)
        let gy = try meshDimension(lengthNM: cell.b.norm, spacingNM: targetGridSpacingNM)
        let gz = try meshDimension(lengthNM: cell.c.norm, spacingNM: targetGridSpacingNM)
        let xy = UInt64(gx) * UInt64(gy)
        let count = xy.multipliedReportingOverflow(by: UInt64(gz))
        guard !count.overflow, count.partialValue > 0 else {
            throw VivoMDRuntimeError.metal("PME grid point count overflow")
        }
        return .init(gridX: gx, gridY: gy, gridZ: gz,
                     gridPointCount: count.partialValue,
                     interpolationOrder: interpolationOrder,
                     ewaldBetaPerNM: beta,
                     tolerance: tolerance,
                     targetGridSpacingNM: targetGridSpacingNM,
                     cellVolumeNM3: abs(determinant),
                     reciprocalA: reciprocalA,
                     reciprocalB: reciprocalB,
                     reciprocalC: reciprocalC)
    }

    private static func solveBeta(cutoffNM: Double, tolerance: Double) throws -> Double {
        // Real-space truncation criterion erfc(beta*r_c)/r_c <= tolerance.
        // Solve monotonically by bisection rather than relying on a fragile closed
        // approximation. This is planning-time work, not a per-step cost.
        func residual(_ beta: Double) -> Double {
            Foundation.erfc(beta * cutoffNM) / cutoffNM
        }
        var low = 0.0
        var high = 1.0 / cutoffNM
        var guardIterations = 0
        while residual(high) > tolerance {
            high *= 2
            guardIterations += 1
            guard high.isFinite, guardIterations < 128 else {
                throw VivoMDRuntimeError.metal("PME Ewald-beta search failed to bracket the requested tolerance")
            }
        }
        for _ in 0..<80 {
            let mid = 0.5 * (low + high)
            if residual(mid) > tolerance { low = mid } else { high = mid }
        }
        guard high.isFinite, high > 0 else {
            throw VivoMDRuntimeError.metal("PME Ewald beta is invalid")
        }
        return high
    }

    private static func meshDimension(lengthNM: Double, spacingNM: Double) throws -> UInt32 {
        guard lengthNM.isFinite, lengthNM > 0 else {
            throw VivoMDRuntimeError.metal("PME lattice-vector length is invalid")
        }
        let requested = max(4, Int(ceil(lengthNM / spacingNM)))
        var value = 1
        while value < requested {
            let doubled = value.multipliedReportingOverflow(by: 2)
            guard !doubled.overflow, doubled.partialValue <= Int(UInt32.max) else {
                throw VivoMDRuntimeError.metal("PME mesh dimension exceeds UInt32")
            }
            value = doubled.partialValue
        }
        return UInt32(value)
    }
}
