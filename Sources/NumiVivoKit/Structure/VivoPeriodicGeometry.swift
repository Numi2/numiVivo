import Foundation

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
