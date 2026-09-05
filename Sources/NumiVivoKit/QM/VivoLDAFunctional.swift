import Foundation

public struct VivoLDAEvaluation: Codable, Sendable, Equatable {
    public let densityPerBohr3: Double
    public let rs: Double?
    public let exchangeEnergyPerElectronHartree: Double
    public let correlationEnergyPerElectronHartree: Double
    public let exchangePotentialHartree: Double
    public let correlationPotentialHartree: Double
    public var xcEnergyPerElectronHartree: Double { exchangeEnergyPerElectronHartree + correlationEnergyPerElectronHartree }
    public var xcPotentialHartree: Double { exchangePotentialHartree + correlationPotentialHartree }
}

/// Spin-unpolarized Dirac exchange + Perdew-Zunger (1981) correlation.
/// Density is total electron density in a0^-3 and energies are Hartree.
public enum VivoPZ81LDA {
    private static let A = 0.0311
    private static let B = -0.0480
    private static let C = 0.0020
    private static let D = -0.0116
    private static let gamma = -0.1423
    private static let beta1 = 1.0529
    private static let beta2 = 0.3334
    private static let exchangeCoefficient = -0.75 * pow(3 / Double.pi, 1.0 / 3.0)

    public static func evaluate(density rho: Double) throws -> VivoLDAEvaluation {
        guard rho.isFinite, rho >= 0 else { throw VivoChemistryError.invalid("LDA density") }
        if rho <= 1e-24 {
            return .init(densityPerBohr3: rho, rs: nil,
                         exchangeEnergyPerElectronHartree: 0, correlationEnergyPerElectronHartree: 0,
                         exchangePotentialHartree: 0, correlationPotentialHartree: 0)
        }
        let rs = pow(3 / (4 * Double.pi * rho), 1.0 / 3.0)
        guard rs.isFinite, rs > 0 else { throw VivoChemistryError.invalid("LDA density parameter") }
        let exchange = exchangeCoefficient * pow(rho, 1.0 / 3.0)
        let exchangePotential = (4.0 / 3.0) * exchange
        let correlation: Double, correlationPotential: Double
        if rs < 1 {
            let logRS = log(rs)
            correlation = A * logRS + B + C * rs * logRS + D * rs
            correlationPotential = A * logRS + B - A / 3
                + (2.0 / 3.0) * C * rs * logRS + (2 * D - C) * rs / 3
        } else {
            let root = sqrt(rs)
            let denominator = 1 + beta1 * root + beta2 * rs
            correlation = gamma / denominator
            correlationPotential = gamma * (1 + (7.0 / 6.0) * beta1 * root + (4.0 / 3.0) * beta2 * rs)
                / (denominator * denominator)
        }
        guard [exchange, exchangePotential, correlation, correlationPotential].allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite LDA functional")
        }
        return .init(densityPerBohr3: rho, rs: rs,
                     exchangeEnergyPerElectronHartree: exchange,
                     correlationEnergyPerElectronHartree: correlation,
                     exchangePotentialHartree: exchangePotential,
                     correlationPotentialHartree: correlationPotential)
    }
}

public struct VivoDFTGridConfiguration: Codable, Sendable, Equatable {
    public var radialPoints: Int
    public var angularPoints: Int
    public var radialScaleBohr: Double
    public var partitionSharpnessPerBohr: Double
    public var maximumGridPoints: Int
    public init(radialPoints: Int = 64, angularPoints: Int = 302, radialScaleBohr: Double = 1,
                partitionSharpnessPerBohr: Double = 2, maximumGridPoints: Int = 2_000_000) {
        self.radialPoints = radialPoints; self.angularPoints = angularPoints
        self.radialScaleBohr = radialScaleBohr; self.partitionSharpnessPerBohr = partitionSharpnessPerBohr
        self.maximumGridPoints = maximumGridPoints
    }
    public func validate(atomCount: Int) throws {
        let (perAtom, overflow1) = radialPoints.multipliedReportingOverflow(by: angularPoints)
        let (total, overflow2) = perAtom.multipliedReportingOverflow(by: atomCount)
        guard atomCount > 0, radialPoints >= 8, radialPoints <= 1024,
              angularPoints >= 26, angularPoints <= 4096,
              !overflow1, !overflow2, total <= maximumGridPoints,
              radialScaleBohr.isFinite, radialScaleBohr > 0,
              partitionSharpnessPerBohr.isFinite, partitionSharpnessPerBohr > 0,
              maximumGridPoints > 0 else {
            throw VivoChemistryError.resourceLimit("DFT quadrature configuration")
        }
    }
}

struct VivoDFTGridPoint: Sendable {
    let position: SIMD3<Double>
    let weightBohr3: Double
}

enum VivoDFTQuadrature {
    static func gaussLegendre(_ n: Int) throws -> (nodes: [Double], weights: [Double]) {
        guard n > 0, n <= 4096 else { throw VivoChemistryError.resourceLimit("Gauss-Legendre order") }
        var nodes = [Double](repeating: 0, count: n), weights = nodes
        let half = (n + 1) / 2
        for i in 0..<half {
            var z = cos(Double.pi * (Double(i) + 0.75) / (Double(n) + 0.5))
            var derivative = 0.0, converged = false
            for _ in 0..<64 {
                var p0 = 1.0, p1 = z
                if n >= 2 {
                    for k in 2...n {
                        let p = (Double(2*k-1) * z * p1 - Double(k-1) * p0) / Double(k)
                        p0 = p1; p1 = p
                    }
                }
                let pn = n == 1 ? z : p1, pn1 = n == 1 ? 1 : p0
                derivative = Double(n) * (z * pn - pn1) / (z*z - 1)
                let next = z - pn / derivative
                if abs(next - z) < 2e-15 { z = next; converged = true; break }
                z = next
            }
            guard converged else { throw VivoChemistryError.convergence("Gauss-Legendre root") }
            var p0 = 1.0, p1 = z
            if n >= 2 { for k in 2...n { let p=(Double(2*k-1)*z*p1-Double(k-1)*p0)/Double(k); p0=p1; p1=p } }
            let pn = n == 1 ? z : p1, pn1 = n == 1 ? 1 : p0
            derivative = Double(n) * (z * pn - pn1) / (z*z - 1)
            let weight = 2 / ((1-z*z) * derivative * derivative)
            nodes[i] = -z; nodes[n-1-i] = z; weights[i] = weight; weights[n-1-i] = weight
        }
        return (nodes, weights)
    }

    static func fibonacciSphere(_ n: Int) -> [(SIMD3<Double>, Double)] {
        let golden = Double.pi * (3 - sqrt(5.0)), weight = 4 * Double.pi / Double(n)
        return (0..<n).map { index in
            let z = 1 - 2 * (Double(index) + 0.5) / Double(n)
            let radius = sqrt(max(0, 1-z*z)), phi = golden * Double(index)
            return (SIMD3<Double>(radius*cos(phi), radius*sin(phi), z), weight)
        }
    }

    static func build(system: VivoElectronicSystem, configuration cfg: VivoDFTGridConfiguration) throws -> [VivoDFTGridPoint] {
        try cfg.validate(atomCount: system.nuclei.count)
        let radial = try gaussLegendre(cfg.radialPoints), angular = fibonacciSphere(cfg.angularPoints)
        var result: [VivoDFTGridPoint] = []
        result.reserveCapacity(system.nuclei.count * cfg.radialPoints * cfg.angularPoints)
        for (atomIndex, nucleus) in system.nuclei.enumerated() {
            for i in radial.nodes.indices {
                let x = 0.5 * (radial.nodes[i] + 1), wx = 0.5 * radial.weights[i]
                let denominator = 1-x, radius = cfg.radialScaleBohr * x / denominator
                let radialWeight = wx * radius * radius * cfg.radialScaleBohr / (denominator * denominator)
                for (direction, angularWeight) in angular {
                    let point = nucleus.positionBohr + radius * direction
                    var distances = [Double](repeating: 0, count: system.nuclei.count), minimum = Double.infinity
                    for (j, other) in system.nuclei.enumerated() {
                        distances[j] = vivoQMNorm(point-other.positionBohr); minimum = min(minimum, distances[j])
                    }
                    var normalization = 0.0
                    for j in distances.indices {
                        distances[j] = exp(-cfg.partitionSharpnessPerBohr * (distances[j]-minimum)); normalization += distances[j]
                    }
                    let weight = radialWeight * angularWeight * distances[atomIndex] / normalization
                    guard weight.isFinite, weight >= 0 else { throw VivoChemistryError.convergence("DFT grid weight") }
                    result.append(.init(position: point, weightBohr3: weight))
                }
            }
        }
        return result
    }
}
