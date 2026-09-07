import Foundation

public struct VivoObservableCovarianceResult: Codable, Sendable, Equatable {
    public let variance: Double
    public let standardDeviation: Double
    public let signedParameterContributions: [Double]
    public let minimumCovarianceEigenvalue: Double
}
/// Local linear uncertainty propagation. Off-diagonal covariance is retained;
/// marginal standard deviations do not supply a missing joint covariance.
public enum VivoObservableCovariance {
    public static func project(derivatives: [Double], covariance: VivoQMMatrix) throws -> VivoObservableCovarianceResult {
        let n = derivatives.count
        guard (1...256).contains(n), covariance.rows == n, covariance.columns == n,
              derivatives.allSatisfy(\.isFinite), covariance.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("observable covariance dimensions")
        }
        let scale = covariance.values.map(abs).max()!
        if scale == 0 { return .init(variance: 0,standardDeviation: 0,
                                    signedParameterContributions: Array(repeating: 0,count: n),minimumCovarianceEigenvalue: 0) }
        for i in 0..<n { for j in 0..<n where abs(covariance[i,j]-covariance[j,i]) > 1e-12*scale {
            throw VivoChemistryError.invalid("nonsymmetric parameter covariance")
        } }
        let normalized = try VivoQMMatrix(rows: n,columns: n,values: covariance.values.map { $0/scale })
        let spectrum = try VivoQMDenseAlgebra.symmetricEigen(normalized)
        guard spectrum.values.first! >= -1e-12 else { throw VivoChemistryError.invalid("parameter covariance is not positive semidefinite") }
        let contributions = (0..<n).map { i in derivatives[i]*(0..<n).reduce(0.0) { $0+covariance[i,$1]*derivatives[$1] } }
        let variance = contributions.reduce(0,+)
        let derivativeNormSquared = derivatives.reduce(0) { $0+$1*$1 }
        guard derivativeNormSquared.isFinite, contributions.allSatisfy(\.isFinite), variance.isFinite,
              variance >= -1e-12*scale*max(1,derivativeNormSquared) else {
            throw VivoChemistryError.convergence("observable covariance contraction")
        }
        return .init(variance: max(0,variance),standardDeviation: sqrt(max(0,variance)),
                     signedParameterContributions: contributions,minimumCovarianceEigenvalue: scale*spectrum.values.first!)
    }
}
