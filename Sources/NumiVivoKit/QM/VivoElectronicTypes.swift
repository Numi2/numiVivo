import Foundation

public enum VivoChemistryError: Error, Sendable, CustomStringConvertible {
    case invalid(String), unsupported(String), resourceLimit(String), convergence(String)
    public var description: String {
        switch self {
        case .invalid(let s): return "Invalid chemistry input: \(s)"
        case .unsupported(let s): return "Unsupported chemistry capability: \(s)"
        case .resourceLimit(let s): return "Chemistry resource limit: \(s)"
        case .convergence(let s): return "Chemistry convergence failure: \(s)"
        }
    }
}

public enum VivoAtomicUnits {
    public static let bohrInNM = 0.0529177210544
    public static let hartreeInKJPerMol = 2625.4996394799
    public static let boltzmannJPerK = 1.380649e-23
    public static let planckJS = 6.62607015e-34
    public static let gasConstantJPerMolK = 8.31446261815324
}

@inline(__always) func vivoQMDot(_ a: SIMD3<Double>, _ b: SIMD3<Double>) -> Double { a.x * b.x + a.y * b.y + a.z * b.z }
@inline(__always) func vivoQMNorm(_ a: SIMD3<Double>) -> Double { sqrt(vivoQMDot(a, a)) }
@inline(__always) func vivoQMFinite(_ a: SIMD3<Double>) -> Bool { a.x.isFinite && a.y.isFinite && a.z.isFinite }

public struct VivoQMNucleus: Codable, Sendable, Equatable {
    public var atomicNumber: Int
    public var positionBohr: SIMD3<Double>
    public var structureAtomIndex: UInt32?
    public init(atomicNumber: Int, positionBohr: SIMD3<Double>, structureAtomIndex: UInt32? = nil) { self.atomicNumber = atomicNumber; self.positionBohr = positionBohr; self.structureAtomIndex = structureAtomIndex }
}
public struct VivoQMPointCharge: Codable, Sendable, Equatable {
    public var chargeE: Double
    public var positionBohr: SIMD3<Double>
    public var classicalParticleIndex: UInt32?
    public init(chargeE: Double, positionBohr: SIMD3<Double>, classicalParticleIndex: UInt32? = nil) { self.chargeE = chargeE; self.positionBohr = positionBohr; self.classicalParticleIndex = classicalParticleIndex }
}
public struct VivoElectronicSystem: Codable, Sendable, Equatable {
    public var nuclei: [VivoQMNucleus]
    public var pointCharges: [VivoQMPointCharge]
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public init(nuclei: [VivoQMNucleus], alphaElectrons: Int, betaElectrons: Int, pointCharges: [VivoQMPointCharge] = []) { self.nuclei = nuclei; self.alphaElectrons = alphaElectrons; self.betaElectrons = betaElectrons; self.pointCharges = pointCharges }
    public func validate() throws {
        guard !nuclei.isEmpty, nuclei.count <= 100_000, alphaElectrons >= 0, betaElectrons >= 0,
              nuclei.allSatisfy({ (1...118).contains($0.atomicNumber) && vivoQMFinite($0.positionBohr) }),
              pointCharges.allSatisfy({ $0.chargeE.isFinite && vivoQMFinite($0.positionBohr) }) else { throw VivoChemistryError.invalid("nuclei, electron counts or embedding charges") }
    }
}

public struct VivoGaussianPrimitive: Codable, Sendable, Equatable { public var exponent: Double; public var coefficient: Double; public init(exponent: Double, coefficient: Double) { self.exponent = exponent; self.coefficient = coefficient } }
public struct VivoGaussianShell: Codable, Sendable, Equatable { public var nucleusIndex: Int; public var angularMomentum: Int; public var primitives: [VivoGaussianPrimitive]; public init(nucleusIndex: Int, angularMomentum: Int, primitives: [VivoGaussianPrimitive]) { self.nucleusIndex = nucleusIndex; self.angularMomentum = angularMomentum; self.primitives = primitives } }
public struct VivoGaussianBasis: Codable, Sendable, Equatable {
    public var identifier: String; public let representation: String; public var shells: [VivoGaussianShell]; public var source: String
    public init(identifier: String, shells: [VivoGaussianShell], source: String) { self.identifier = identifier; representation = "normalized-cartesian"; self.shells = shells; self.source = source }
    public func validate(nucleusCount: Int) throws {
        guard !identifier.isEmpty, !source.isEmpty, representation == "normalized-cartesian", !shells.isEmpty else { throw VivoChemistryError.invalid("basis identity, source or representation") }
        for s in shells { guard s.nucleusIndex >= 0, s.nucleusIndex < nucleusCount, (0...3).contains(s.angularMomentum), !s.primitives.isEmpty, s.primitives.allSatisfy({ $0.exponent.isFinite && $0.exponent > 0 && $0.coefficient.isFinite }) else { throw VivoChemistryError.invalid("basis shell") } }
    }
    public static func hydrogenSTO3G(nucleusIndices: [Int]) -> Self {
        let e = [3.425250914, 0.6239137298, 0.1688554040], c = [0.1543289673, 0.5353281423, 0.4446345422]
        return .init(identifier: "STO-3G-H-v1", shells: nucleusIndices.map { i in .init(nucleusIndex: i, angularMomentum: 0, primitives: zip(e,c).map { .init(exponent: $0.0, coefficient: $0.1) }) }, source: "Basis Set Exchange STO-3G version 1")
    }
}

public struct VivoChemistryBudget: Codable, Sendable, Equatable {
    public var maximumBytes: Int; public var maximumBasisFunctions: Int; public var maximumDeterminants: Int; public var maximumOperatorApplications: Int
    public init(maximumBytes: Int = 256 * 1024 * 1024, maximumBasisFunctions: Int = 64, maximumDeterminants: Int = 512, maximumOperatorApplications: Int = 100_000_000) { self.maximumBytes = maximumBytes; self.maximumBasisFunctions = maximumBasisFunctions; self.maximumDeterminants = maximumDeterminants; self.maximumOperatorApplications = maximumOperatorApplications }
    public func validate() throws { guard maximumBytes > 0, maximumBasisFunctions > 0, maximumDeterminants > 0, maximumOperatorApplications > 0 else { throw VivoChemistryError.invalid("nonpositive resource budget") } }
}
