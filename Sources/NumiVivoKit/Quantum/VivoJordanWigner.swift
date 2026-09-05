import Foundation

public struct VivoComplex64: Codable, Sendable, Equatable {
    public var real: Double
    public var imaginary: Double
    public init(_ real: Double, _ imaginary: Double = 0) { self.real = real; self.imaginary = imaginary }
    public var conjugate: Self { .init(real,-imaginary) }
    public var squaredNorm: Double { real*real+imaginary*imaginary }
    public var isFinite: Bool { real.isFinite && imaginary.isFinite }
    public static func +(a: Self,b: Self) -> Self { .init(a.real+b.real,a.imaginary+b.imaginary) }
    public static func *(a: Self,b: Self) -> Self { .init(a.real*b.real-a.imaginary*b.imaginary,a.real*b.imaginary+a.imaginary*b.real) }
    public static func *(a: Self,b: Double) -> Self { .init(a.real*b,a.imaginary*b) }
}

/// P = tensor product of I/X/Y/Z. (x,z)=(1,1) means Hermitian Y,
/// not XZ. Qubit zero is the least-significant occupation bit.
public struct VivoPauliWord: Codable, Sendable, Hashable {
    public let xMask: UInt64
    public let zMask: UInt64
    public init(xMask: UInt64, zMask: UInt64) { self.xMask = xMask; self.zMask = zMask }
    public func qubitWiseCommutes(with other: Self) -> Bool {
        let common = (xMask|zMask)&(other.xMask|other.zMask)
        return ((xMask^other.xMask)|(zMask^other.zMask))&common == 0
    }
    public func multiplied(by b: Self) -> (Self,VivoComplex64) {
        let result = Self(xMask:xMask^b.xMask,zMask:zMask^b.zMask)
        // P(x,z)=i^(popcount(x&z)) X^x Z^z.
        let exponent = (xMask&zMask).nonzeroBitCount+(b.xMask&b.zMask).nonzeroBitCount
            - (result.xMask&result.zMask).nonzeroBitCount+2*(zMask&b.xMask).nonzeroBitCount
        let phase = ((exponent%4)+4)%4
        return (result,[VivoComplex64(1),.init(0,1),.init(-1),.init(0,-1)][phase])
    }
}
public struct VivoPauliTerm: Codable, Sendable {
    public let word: VivoPauliWord
    public let coefficient: Double
}
public struct VivoQubitHamiltonian: Codable, Sendable {
    public let qubitCount: Int
    public let terms: [VivoPauliTerm]
    public let discardedCoefficientL1: Double
    public func validate() throws {
        guard (1...24).contains(qubitCount), terms.count <= 1_000_000,
              Set(terms.map(\.word)).count == terms.count,
              discardedCoefficientL1.isFinite, discardedCoefficientL1 >= 0,
              terms.allSatisfy({ $0.coefficient.isFinite && (($0.word.xMask|$0.word.zMask)>>qubitCount) == 0 }) else {
            throw VivoQMError.invalid("Pauli Hamiltonian storage or qubit range")
        }
    }
    /// Deterministic greedy QWC grouping, not general commuting-group diagonalization.
    public func qubitWiseCommutingGroups() throws -> [[Int]] {
        try validate()
        guard terms.count <= 20_000 else { throw VivoQMError.capacity("quadratic reference measurement grouping") }
        var groups: [[Int]] = []
        for index in terms.indices {
            if let group = groups.firstIndex(where: { $0.allSatisfy { terms[$0].word.qubitWiseCommutes(with:terms[index].word) } }) {
                groups[group].append(index)
            } else { groups.append([index]) }
        }
        return groups
    }
}

public enum VivoJordanWigner {
    public static func map(_ h: VivoEmbeddedHamiltonian, coefficientThreshold: Double = 1e-13) throws -> VivoQubitHamiltonian {
        try h.validate()
        guard coefficientThreshold.isFinite, (0...1e-6).contains(coefficientThreshold) else { throw VivoQMError.invalid("Pauli coefficient threshold") }
        let fermions = try VivoFermionAlgebra.terms(for:h)
        var accumulated: [VivoPauliWord:VivoComplex64] = [.init(xMask:0,zMask:0):.init(h.constantEnergyHartree)]
        for term in fermions {
            var expansion: [(VivoPauliWord,VivoComplex64)] = [(.init(xMask:0,zMask:0),.init(term.coefficient))]
            for ladder in term.operators {
                let bit = UInt64(1)<<ladder.mode, parity = bit-1
                let factors: [(VivoPauliWord,VivoComplex64)] = [
                    (.init(xMask:bit,zMask:parity),.init(0.5)),
                    (.init(xMask:bit,zMask:parity|bit),.init(0,ladder.creation ? -0.5 : 0.5))]
                var next: [(VivoPauliWord,VivoComplex64)] = []
                for (a,ca) in expansion { for (b,cb) in factors {
                    let (product,phase) = a.multiplied(by:b)
                    next.append((product,ca*cb*phase))
                } }
                expansion = next
            }
            for (word,value) in expansion { accumulated[word] = (accumulated[word] ?? .init(0))+value }
        }
        var terms: [VivoPauliTerm] = [], discarded = 0.0
        for word in accumulated.keys.sorted(by: { $0.xMask == $1.xMask ? $0.zMask < $1.zMask : $0.xMask < $1.xMask }) {
            let value = accumulated[word]!
            guard value.isFinite, abs(value.imaginary) < 1e-10*max(1,abs(value.real)) else {
                throw VivoQMError.invalid("Jordan-Wigner output is not Hermitian")
            }
            if abs(value.real) <= coefficientThreshold { discarded += abs(value.real) }
            else { terms.append(.init(word:word,coefficient:value.real)) }
        }
        let result = VivoQubitHamiltonian(qubitCount:2*h.orbitalCount,terms:terms,discardedCoefficientL1:discarded)
        try result.validate()
        return result
    }
}

/// FP64 CPU conformance engine for the portable qubit boundary. This is NOT
/// described as a Metal simulator; production GPU simulation remains separate.
public struct VivoQuantumReferenceState: Codable, Sendable {
    public let qubitCount: Int
    public private(set) var amplitudes: [VivoComplex64]
    public init(qubitCount: Int, amplitudes: [VivoComplex64]) throws {
        self.qubitCount = qubitCount; self.amplitudes = amplitudes; try validate()
    }
    public init(determinantState: VivoDeterminantState) throws {
        try determinantState.validate()
        let count = 2*determinantState.spatialOrbitalCount
        guard count <= 20 else { throw VivoQMError.capacity("reference state vector is limited to 20 qubits") }
        qubitCount = count; amplitudes = Array(repeating:.init(0),count:1<<count)
        for (index,bits) in determinantState.determinants.enumerated() { amplitudes[Int(bits)] = .init(determinantState.amplitudes[index]) }
        try validate()
    }
    public func validate() throws {
        guard (1...20).contains(qubitCount), amplitudes.count == 1<<qubitCount,
              amplitudes.allSatisfy(\.isFinite), abs(amplitudes.reduce(0) { $0+$1.squaredNorm }-1) < 1e-8 else {
            throw VivoQMError.invalid("reference state-vector dimensions, finite values, or normalization")
        }
    }
    private func action(_ word: VivoPauliWord, on index: Int) -> (Int,VivoComplex64) {
        let ys = (word.xMask&word.zMask).nonzeroBitCount%4
        let phase = [VivoComplex64(1),.init(0,1),.init(-1),.init(0,-1)][ys]
        let sign = (UInt64(index)&word.zMask).nonzeroBitCount%2 == 0 ? 1.0 : -1.0
        return (index^Int(word.xMask),phase*sign)
    }
    public func expectation(of h: VivoQubitHamiltonian, maximumTermAmplitudes: Int = 100_000_000) throws -> Double {
        try validate(); try h.validate()
        guard h.qubitCount == qubitCount, maximumTermAmplitudes > 0,
              h.terms.count <= maximumTermAmplitudes/amplitudes.count else { throw VivoQMError.capacity("Pauli expectation dimension or work budget") }
        var energy = 0.0
        for term in h.terms {
            var value = VivoComplex64(0)
            for i in amplitudes.indices {
                let (target,phase) = action(term.word,on:i)
                value = value+amplitudes[target].conjugate*phase*amplitudes[i]
            }
            guard abs(value.imaginary) < 1e-8 else { throw VivoQMError.convergence("imaginary expectation of Hermitian Pauli word") }
            energy += term.coefficient*value.real
        }
        return energy
    }
    /// Exact exp(-i theta P / 2), useful for excitation/ADAPT reference checks.
    /// Does not claim to construct UCCSD or select ADAPT operators by itself.
    public mutating func applyPauliRotation(_ word: VivoPauliWord, angle: Double) throws {
        try validate()
        guard angle.isFinite, (word.xMask|word.zMask)>>qubitCount == 0 else { throw VivoQMError.invalid("Pauli rotation") }
        let c = cos(angle/2), s = VivoComplex64(0,-sin(angle/2))
        var next = amplitudes.map { $0*c }
        for i in amplitudes.indices {
            let (target,phase) = action(word,on:i)
            next[target] = next[target]+s*phase*amplitudes[i]
        }
        amplitudes = next; try validate()
    }
}
