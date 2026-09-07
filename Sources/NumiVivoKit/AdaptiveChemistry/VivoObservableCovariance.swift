import Foundation

public struct VivoObservableCovarianceProjection: Codable, Sendable, Equatable {
    public let variance: Double
    public let standardDeviation: Double
    public let signedParameterContributions: [Double]
    public let minimumCovarianceEigenvalue: Double
}
public enum VivoObservableCovariance {
    /// Local linear propagation g^T Cov g. Off-diagonal terms are retained; a
    /// missing joint covariance is never replaced by independent marginals.
    public static func project(derivatives: [Double],covariance: VivoQMMatrix) throws -> VivoObservableCovarianceProjection {
        let n = derivatives.count
        guard (1...256).contains(n), covariance.rows == n, covariance.columns == n,
              derivatives.allSatisfy(\.isFinite), covariance.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("observable covariance dimensions")
        }
        let scale = covariance.values.map(abs).max()!
        if scale == 0 { return .init(variance: 0,standardDeviation: 0,
                                    signedParameterContributions: Array(repeating: 0,count: n),minimumCovarianceEigenvalue: 0) }
        for i in 0..<n { for j in 0..<n where abs(covariance[i,j]-covariance[j,i]) > 1e-12*scale {
            throw VivoChemistryError.invalid("observable covariance is not symmetric")
        } }
        let spectrum = try VivoQMDenseAlgebra.symmetricEigen(covariance.scaled(1/scale))
        guard spectrum.values.first! >= -1e-12 else { throw VivoChemistryError.invalid("observable covariance is not positive semidefinite") }
        let contributions = (0..<n).map { i in derivatives[i]*(0..<n).reduce(0.0) { $0+covariance[i,$1]*derivatives[$1] } }
        let variance = contributions.reduce(0,+)
        guard variance.isFinite, variance >= -1e-12*scale*max(1,derivatives.reduce(0) { $0+$1*$1 }) else {
            throw VivoChemistryError.convergence("observable covariance projection overflow")
        }
        return .init(variance: max(0,variance),standardDeviation: sqrt(max(0,variance)),
                     signedParameterContributions: contributions,minimumCovarianceEigenvalue: scale*spectrum.values.first!)
    }
}
