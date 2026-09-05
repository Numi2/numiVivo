import Foundation

public struct VivoQuantumBudget: Codable, Sendable, Equatable {
    public var maximumStateBytes: Int
    public var maximumQubits: Int
    public var maximumTerms: Int
    public var maximumShots: Int
    public init(maximumStateBytes: Int = 512 * 1024 * 1024, maximumQubits: Int = 24,
                maximumTerms: Int = 2_000_000, maximumShots: Int = 1_000_000) {
        self.maximumStateBytes = maximumStateBytes; self.maximumQubits = maximumQubits
        self.maximumTerms = maximumTerms; self.maximumShots = maximumShots
    }
    public func validate() throws {
        guard maximumStateBytes > 0, (1...62).contains(maximumQubits), maximumTerms > 0, maximumShots > 0 else {
            throw VivoChemistryError.invalid("quantum resource budget")
        }
    }
    public func stateCount(qubits: Int) throws -> Int {
        try validate()
        guard qubits >= 0, qubits <= maximumQubits, qubits < Int.bitWidth - 1 else {
            throw VivoChemistryError.resourceLimit("state-vector qubit budget")
        }
        let count = 1 << qubits
        let (bytes, overflow) = count.multipliedReportingOverflow(by: MemoryLayout<VivoComplex>.stride)
        guard !overflow, bytes <= maximumStateBytes else {
            throw VivoChemistryError.resourceLimit("state vector exceeds \(maximumStateBytes) bytes")
        }
        return count
    }
}

public struct VivoStateVector: Codable, Sendable, Equatable {
    public let qubitCount: Int
    public private(set) var amplitudes: [VivoComplex]

    public init(qubitCount: Int, amplitudes: [VivoComplex], budget: VivoQuantumBudget = .init()) throws {
        let count = try budget.stateCount(qubits: qubitCount)
        guard amplitudes.count == count,
              amplitudes.allSatisfy({ $0.real.isFinite && $0.imaginary.isFinite }) else {
            throw VivoChemistryError.invalid("state-vector dimensions or amplitudes")
        }
        let norm = amplitudes.reduce(0.0) { $0 + $1.magnitudeSquared }
        guard norm.isFinite, abs(norm - 1) <= 1e-10 else {
            throw VivoChemistryError.invalid("state vector must be normalized")
        }
        self.qubitCount = qubitCount; self.amplitudes = amplitudes
    }

    public static func basis(qubitCount: Int, index: UInt64 = 0,
                             budget: VivoQuantumBudget = .init()) throws -> Self {
        let count = try budget.stateCount(qubits: qubitCount)
        guard index < UInt64(count) else { throw VivoChemistryError.invalid("basis-state index") }
        var amplitudes = [VivoComplex](repeating: .zero, count: count)
        amplitudes[Int(index)] = .init(1)
        return try .init(qubitCount: qubitCount, amplitudes: amplitudes, budget: budget)
    }

    public static func fromCI(_ state: VivoCIState, budget: VivoQuantumBudget = .init()) throws -> Self {
        try state.validate()
        let qubits = 2 * state.orbitalCount
        let count = try budget.stateCount(qubits: qubits)
        var amplitudes = [VivoComplex](repeating: .zero, count: count)
        for (determinant, coefficient) in zip(state.determinants, state.coefficients) {
            guard determinant < UInt64(count) else { throw VivoChemistryError.invalid("CI determinant exceeds state-vector width") }
            amplitudes[Int(determinant)] = .init(coefficient)
        }
        return try .init(qubitCount: qubits, amplitudes: amplitudes, budget: budget)
    }

    public func expectation(_ key: VivoPauliKey) throws -> VivoComplex {
        let validMask = qubitCount == 64 ? UInt64.max : (UInt64(1) << qubitCount) - 1
        guard key.xMask & ~validMask == 0, key.zMask & ~validMask == 0 else {
            throw VivoChemistryError.invalid("Pauli operator exceeds state-vector width")
        }
        let yCount = (key.xMask & key.zMask).nonzeroBitCount
        let yPhase: VivoComplex
        switch yCount & 3 {
        case 0: yPhase = .init(1)
        case 1: yPhase = .init(0, 1)
        case 2: yPhase = .init(-1)
        default: yPhase = .init(0, -1)
        }
        var value = VivoComplex.zero
        for source in amplitudes.indices {
            let sourceBits = UInt64(source)
            let destination = Int(sourceBits ^ key.xMask)
            let parity = (sourceBits & key.zMask).nonzeroBitCount & 1
            let phase = parity == 0 ? yPhase : yPhase * -1
            value = value + amplitudes[destination].conjugate * phase * amplitudes[source]
        }
        guard value.real.isFinite, value.imaginary.isFinite else {
            throw VivoChemistryError.convergence("nonfinite Pauli expectation")
        }
        return value
    }

    public func expectation(_ hamiltonian: VivoPauliHamiltonian,
                            budget: VivoQuantumBudget = .init()) throws -> Double {
        try hamiltonian.validate()
        try budget.validate()
        guard hamiltonian.qubitCount == qubitCount, hamiltonian.terms.count <= budget.maximumTerms else {
            throw VivoChemistryError.resourceLimit("Hamiltonian/state-vector contract")
        }
        var total = VivoComplex.zero
        for term in hamiltonian.terms {
            total = total + term.coefficient * (try expectation(term.key))
        }
        guard total.real.isFinite, abs(total.imaginary) <= 1e-9 * max(1, abs(total.real)) else {
            throw VivoChemistryError.invalid("Hermitian Hamiltonian expectation is not real")
        }
        return total.real
    }

    public mutating func applyRY(qubit: Int, angle: Double) throws {
        guard qubit >= 0, qubit < qubitCount, angle.isFinite else { throw VivoChemistryError.invalid("RY gate") }
        let bit = 1 << qubit, c = cos(0.5 * angle), s = sin(0.5 * angle)
        for base in amplitudes.indices where base & bit == 0 {
            let pair = base | bit, a = amplitudes[base], b = amplitudes[pair]
            amplitudes[base] = a * c - b * s
            amplitudes[pair] = a * s + b * c
        }
    }

    public mutating func applyRZ(qubit: Int, angle: Double) throws {
        guard qubit >= 0, qubit < qubitCount, angle.isFinite else { throw VivoChemistryError.invalid("RZ gate") }
        let bit = 1 << qubit, half = 0.5 * angle
        let low = VivoComplex(cos(half), -sin(half)), high = VivoComplex(cos(half), sin(half))
        for i in amplitudes.indices { amplitudes[i] = amplitudes[i] * (i & bit == 0 ? low : high) }
    }

    public mutating func applyCNOT(control: Int, target: Int) throws {
        guard control >= 0, target >= 0, control < qubitCount, target < qubitCount, control != target else {
            throw VivoChemistryError.invalid("CNOT gate")
        }
        let c = 1 << control, t = 1 << target
        for i in amplitudes.indices where i & c != 0 && i & t == 0 { amplitudes.swapAt(i, i | t) }
    }

    public func sample(shots: Int, seed: UInt64 = 0x4e554d495649564f,
                       budget: VivoQuantumBudget = .init()) throws -> [UInt64:Int] {
        try budget.validate()
        guard shots >= 0, shots <= budget.maximumShots else { throw VivoChemistryError.resourceLimit("shot budget") }
        if shots == 0 { return [:] }
        var cumulative = [Double](repeating: 0, count: amplitudes.count), running = 0.0
        for i in amplitudes.indices { running += amplitudes[i].magnitudeSquared; cumulative[i] = running }
        guard abs(running - 1) <= 1e-10 else { throw VivoChemistryError.invalid("state norm before sampling") }
        cumulative[cumulative.count - 1] = 1
        var random = VivoSplitMix64(seed: seed), counts: [UInt64:Int] = [:]
        for _ in 0..<shots {
            let u = random.nextUnit()
            var low = 0, high = cumulative.count - 1
            while low < high {
                let mid = low + (high - low) / 2
                if cumulative[mid] > u { high = mid } else { low = mid + 1 }
            }
            counts[UInt64(low), default: 0] += 1
        }
        return counts
    }
}

private struct VivoSplitMix64 {
    var state: UInt64
    init(seed: UInt64) { state = seed }
    mutating func next() -> UInt64 {
        state &+= 0x9e3779b97f4a7c15
        var z = state
        z = (z ^ (z >> 30)) &* 0xbf58476d1ce4e5b9
        z = (z ^ (z >> 27)) &* 0x94d049bb133111eb
        return z ^ (z >> 31)
    }
    mutating func nextUnit() -> Double { Double(next() >> 11) * 0x1.0p-53 }
}

public struct VivoPauliMeasurementGroup: Codable, Sendable, Equatable {
    public let basis: VivoPauliKey
    public let termIndices: [Int]
}

public enum VivoQubitWiseCommutingGrouping {
    private static func compatible(_ a: VivoPauliKey, _ b: VivoPauliKey, qubits: Int) -> Bool {
        for q in 0..<qubits {
            let bit = UInt64(1) << q
            let ax = a.xMask & bit != 0, az = a.zMask & bit != 0
            let bx = b.xMask & bit != 0, bz = b.zMask & bit != 0
            let aIdentity = !ax && !az, bIdentity = !bx && !bz
            if !aIdentity && !bIdentity && (ax != bx || az != bz) { return false }
        }
        return true
    }
    private static func merged(_ a: VivoPauliKey, _ b: VivoPauliKey, qubits: Int) -> VivoPauliKey {
        var x: UInt64 = 0, z: UInt64 = 0
        for q in 0..<qubits {
            let bit = UInt64(1) << q
            let ax = a.xMask & bit != 0, az = a.zMask & bit != 0
            let bx = b.xMask & bit != 0, bz = b.zMask & bit != 0
            let useA = ax || az
            let px = useA ? ax : bx, pz = useA ? az : bz
            if px { x |= bit }; if pz { z |= bit }
        }
        return .init(xMask: x, zMask: z)
    }
    public static func group(_ hamiltonian: VivoPauliHamiltonian) throws -> [VivoPauliMeasurementGroup] {
        try hamiltonian.validate()
        let order = hamiltonian.terms.indices.sorted {
            let wa = (hamiltonian.terms[$0].key.xMask | hamiltonian.terms[$0].key.zMask).nonzeroBitCount
            let wb = (hamiltonian.terms[$1].key.xMask | hamiltonian.terms[$1].key.zMask).nonzeroBitCount
            return wa == wb ? $0 < $1 : wa > wb
        }
        var groups: [(VivoPauliKey,[Int])] = []
        for index in order {
            let key = hamiltonian.terms[index].key
            if let slot = groups.indices.first(where: { compatible(groups[$0].0, key, qubits: hamiltonian.qubitCount) }) {
                groups[slot].0 = merged(groups[slot].0, key, qubits: hamiltonian.qubitCount)
                groups[slot].1.append(index)
            } else { groups.append((key,[index])) }
        }
        return groups.map { .init(basis: $0.0, termIndices: $0.1.sorted()) }
    }
}
