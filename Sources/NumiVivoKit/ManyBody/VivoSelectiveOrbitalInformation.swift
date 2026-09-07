import Foundation

public struct VivoInformationOrbitalPair: Codable, Sendable, Equatable, Hashable {
    public let first: Int
    public let second: Int
    public init(_ first: Int, _ second: Int) {
        self.first = min(first,second); self.second = max(first,second)
    }
}
public struct VivoOrbitalInformationSelection: Codable, Sendable, Equatable {
    public let orbitals: [Int]
    public let pairs: [VivoInformationOrbitalPair]
    public init(orbitals: [Int], pairs: [VivoInformationOrbitalPair] = []) {
        self.orbitals = orbitals; self.pairs = pairs
    }
    public func validate(orbitalCount n: Int) throws {
        guard (1...31).contains(n), orbitals.count <= n, pairs.count <= n*(n-1)/2,
              !orbitals.isEmpty || !pairs.isEmpty,
              Set(orbitals).count == orbitals.count, Set(pairs).count == pairs.count,
              orbitals.allSatisfy({ (0..<n).contains($0) }),
              pairs.allSatisfy({ $0.first >= 0 && $0.second < n && $0.first < $0.second }) else {
            throw VivoChemistryError.invalid("selective orbital-information indices, duplicates or capacities")
        }
    }
    public static func all(orbitalCount n: Int) throws -> Self {
        guard (1...31).contains(n) else { throw VivoChemistryError.invalid("orbital-information dimension") }
        return .init(orbitals: Array(0..<n),pairs: (0..<n).flatMap { i in
            ((i+1)..<n).map { VivoInformationOrbitalPair(i,$0) }
        })
    }
}
public struct VivoSingleOrbitalInformation: Codable, Sendable, Equatable {
    public let orbital: Int
    public let occupation: Double
    public let doubleOccupation: Double
    public let entropyNats: Double
}
public struct VivoPairOrbitalInformation: Codable, Sendable, Equatable {
    public let pair: VivoInformationOrbitalPair
    public let jointEntropyNats: Double
    /// I(i:j) = S(i)+S(j)-S(ij), natural logarithm, no factor of one half.
    /// Total orbital correlation, not a measure of purely quantum entanglement.
    public let mutualInformationNats: Double
}
public struct VivoSelectiveOrbitalInformationResult: Codable, Sendable, Equatable {
    public let schema: String
    public let orbitalCount: Int
    public let selection: VivoOrbitalInformationSelection
    public let singles: [VivoSingleOrbitalInformation]
    public let pairs: [VivoPairOrbitalInformation]
    public let completeSingleOrbitalCoverage: Bool
    public let completePairCoverage: Bool
    public let marginalEvaluations: Int
    /// Conservative work reservation, not a hardware FLOP measurement.
    public let reservedPrimitiveWork: Int
}
public struct VivoOrbitalBoundaryInformation: Codable, Sendable, Equatable {
    public let inside: [Int]
    public let outside: [Int]
    public let observedMutualInformationSumNats: Double
    public let measuredPairs: Int
    public let possiblePairs: Int
    public var complete: Bool { measuredPairs == possiblePairs }
}
/// Bounded queries around the existing exact CI marginal implementation. This
/// does not infer two-spatial-orbital density operators from ordinary 1/2-RDMs.
/// Missing pairs are unexamined, not zero. Correlation cuts are diagnostics,
/// never calibrated energy/rate error bounds.
public enum VivoSelectiveOrbitalInformation {
    public static func analyze(_ state: VivoCIState, selection: VivoOrbitalInformationSelection,
                               budget: VivoChemistryBudget = .init()) throws -> VivoSelectiveOrbitalInformationResult {
        try state.validate(budget: budget); try selection.validate(orbitalCount: state.orbitalCount)
        let singles = Set(selection.orbitals + selection.pairs.flatMap { [$0.first,$0.second] }).sorted()
        let pairs = selection.pairs.sorted { $0.first == $1.first ? $0.second < $1.second : $0.first < $1.first }
        let evaluations = singles.count+pairs.count
        // Each pair marginal has at most D environment groups of 16 amplitudes,
        // a 16x16 density and bounded 16x16 eigensolver work. Reserve hashing,
        // amplitude copying, accumulation and eigen sweeps conservatively.
        _ = try budget.elements([state.determinants.count,128],simultaneousArrays: 4)
        let perMarginal = state.determinants.count.multipliedReportingOverflow(by: 1024)
        let withEigen = perMarginal.partialValue.addingReportingOverflow(512_000)
        let work = withEigen.partialValue.multipliedReportingOverflow(by: evaluations)
        guard !perMarginal.overflow, !withEigen.overflow, !work.overflow,
              work.partialValue <= budget.maximumOperatorApplications else {
            throw VivoChemistryError.resourceLimit("selective orbital-information aggregate work reservation")
        }
        var information: [VivoSingleOrbitalInformation] = [], entropies: [Int:Double] = [:]
        for i in singles {
            let rho = try VivoCIDensityMatrices.orbitalMarginal(state,orbitals: [i],budget: budget)
            let entropy = try VivoCIDensityMatrices.entropy(rho)
            let occupation = rho[1,1]+rho[2,2]+2*rho[3,3]
            guard occupation.isFinite, occupation >= -1e-10, occupation <= 2+1e-10 else {
                throw VivoChemistryError.invalid("orbital marginal occupation")
            }
            information.append(.init(orbital: i,occupation: occupation,doubleOccupation: rho[3,3],entropyNats: entropy))
            entropies[i] = entropy
        }
        var pairInformation: [VivoPairOrbitalInformation] = []
        for pair in pairs {
            let rho = try VivoCIDensityMatrices.orbitalMarginal(state,orbitals: [pair.first,pair.second],budget: budget)
            let joint = try VivoCIDensityMatrices.entropy(rho)
            let value = entropies[pair.first]!+entropies[pair.second]!-joint
            guard value.isFinite, value >= -1e-9 else { throw VivoChemistryError.invalid("negative orbital mutual information") }
            pairInformation.append(.init(pair: pair,jointEntropyNats: joint,mutualInformationNats: max(0,value)))
        }
        return .init(schema: "numivivo.org/selective-orbital-information/v1",orbitalCount: state.orbitalCount,
            selection: selection,singles: information,pairs: pairInformation,
            completeSingleOrbitalCoverage: singles.count == state.orbitalCount,
            completePairCoverage: pairs.count == state.orbitalCount*(state.orbitalCount-1)/2,
            marginalEvaluations: evaluations,reservedPrimitiveWork: work.partialValue)
    }
    public static func boundary(_ report: VivoSelectiveOrbitalInformationResult, inside: [Int], outside: [Int]) throws -> VivoOrbitalBoundaryInformation {
        let n = report.orbitalCount
        guard report.schema == "numivivo.org/selective-orbital-information/v1",
              Set(inside).count == inside.count, Set(outside).count == outside.count,
              !inside.isEmpty, !outside.isEmpty, Set(inside).isDisjoint(with: Set(outside)),
              (inside+outside).allSatisfy({ (0..<n).contains($0) }),
              Set(report.pairs.map(\.pair)).count == report.pairs.count,
              Set(report.pairs.map(\.pair)) == Set(report.selection.pairs) else {
            throw VivoChemistryError.invalid("orbital-information boundary")
        }
        try report.selection.validate(orbitalCount: n)
        let a = Set(inside), b = Set(outside)
        var sum = 0.0, count = 0
        for edge in report.pairs {
            guard edge.mutualInformationNats.isFinite, edge.mutualInformationNats >= 0,
                  edge.pair.first >= 0, edge.pair.first < edge.pair.second, edge.pair.second < n else {
                throw VivoChemistryError.invalid("orbital-information edge")
            }
            if (a.contains(edge.pair.first) && b.contains(edge.pair.second)) ||
               (b.contains(edge.pair.first) && a.contains(edge.pair.second)) {
                sum += edge.mutualInformationNats; count += 1
            }
        }
        return .init(inside: inside.sorted(),outside: outside.sorted(),observedMutualInformationSumNats: sum,
                     measuredPairs: count,possiblePairs: inside.count*outside.count)
    }
    public static func validate(_ result: VivoSelectiveOrbitalInformationResult, state: VivoCIState,
                                selection: VivoOrbitalInformationSelection,
                                budget: VivoChemistryBudget = .init()) throws {
        guard result == (try analyze(state,selection: selection,budget: budget)) else {
            throw VivoChemistryError.invalid("orbital information differs from state reconstruction")
        }
    }
}
