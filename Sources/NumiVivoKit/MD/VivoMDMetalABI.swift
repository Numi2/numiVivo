import Foundation

struct VivoMDMetalCommand {
    var particleCount: UInt32
    var typeCount: UInt32
    var electrostatics: UInt32
    var periodic: UInt32
    var dtPS: Float
    var cutoffNM: Float
    var coulombPrefactor: Float
    var reactionFieldK: Float
    var reactionFieldC: Float
    var minimumDistanceNM: Float
    var constraintTolerance: Float
    var reserved0: Float = 0
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
            ("Command", MemoryLayout<VivoMDMetalCommand>.stride, 144),
            ("Status", MemoryLayout<VivoMDMetalStatus>.stride, 16),
            ("Bond", MemoryLayout<VivoMDPackedSystem.Bond>.stride, 16),
            ("Angle", MemoryLayout<VivoMDPackedSystem.Angle>.stride, 32),
            ("Torsion", MemoryLayout<VivoMDPackedSystem.Torsion>.stride, 32),
            ("Constraint", MemoryLayout<VivoMDPackedSystem.Constraint>.stride, 16),
            ("Incidence", MemoryLayout<VivoMDPackedSystem.Incidence>.stride, 8),
            ("PairException", MemoryLayout<VivoMDPackedSystem.PairException>.stride, 32)
        ]
        for (name, actual, required) in expected where actual != required {
            throw VivoMDRuntimeError.metal("Swift/Metal ABI mismatch for \(name): \(actual) != \(required)")
        }
    }

    static func command(packed: VivoMDPackedSystem,
                        configuration: VivoMDConfiguration,
                        cell: VivoPeriodicCell?) throws -> VivoMDMetalCommand {
        let mode: UInt32
        switch configuration.electrostatics {
        case .cutoff: mode = 0
        case .reactionField: mode = 1
        case .pme: mode = 2
        }
        let eps = configuration.relativeDielectric, cutoff = configuration.cutoffNM
        var krf: Float = 0, crf: Float = 0
        if configuration.electrostatics == .reactionField {
            let external = configuration.reactionFieldDielectric
            krf = Float((external - eps) / (2 * external + eps) / pow(cutoff, 3))
            crf = Float(3 * external / (2 * external + eps) / cutoff)
        }
        let reciprocal = try cell.map(reciprocalCell)
        func f4(_ value: VivoVector3D?) -> SIMD4<Float> {
            guard let value else { return .zero }
            return .init(Float(value.x), Float(value.y), Float(value.z), 0)
        }
        return .init(particleCount: packed.particleCount, typeCount: packed.typeCount,
                     electrostatics: mode, periodic: cell == nil ? 0 : 1,
                     dtPS: Float(configuration.timeStepPS), cutoffNM: Float(cutoff),
                     coulombPrefactor: Float(138.935456 / eps), reactionFieldK: krf,
                     reactionFieldC: crf, minimumDistanceNM: 1e-5,
                     constraintTolerance: Float(configuration.constraintTolerance),
                     cellA: f4(cell?.a), cellB: f4(cell?.b), cellC: f4(cell?.c),
                     reciprocalA: f4(reciprocal?.0), reciprocalB: f4(reciprocal?.1),
                     reciprocalC: f4(reciprocal?.2))
    }

    static func validateMinimumImageCutoff(_ cutoff: Double, cell: VivoPeriodicCell) throws {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        let heights = [abs(determinant) / cell.b.cross(cell.c).norm,
                       abs(determinant) / cell.c.cross(cell.a).norm,
                       abs(determinant) / cell.a.cross(cell.b).norm]
        guard heights.allSatisfy({ $0.isFinite && $0 > 0 }),
              cutoff < 0.5 * (heights.min() ?? 0) else {
            throw VivoMDRuntimeError.unsupported([
                "cutoff must be less than half the shortest periodic face height for minimum-image execution"
            ])
        }
    }

    private static func reciprocalCell(_ cell: VivoPeriodicCell) throws -> (VivoVector3D,VivoVector3D,VivoVector3D) {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else {
            throw VivoMDRuntimeError.metal("singular periodic cell")
        }
        return (cell.b.cross(cell.c) / determinant,
                cell.c.cross(cell.a) / determinant,
                cell.a.cross(cell.b) / determinant)
    }
}
