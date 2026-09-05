import Foundation

public struct VivoFermionicExcitation: Codable, Sendable, Equatable, Hashable {
    /// Right-hand annihilation modes, in operator application order.
    public let occupiedModes: [Int]
    /// Created modes, in the corresponding target set. Creation is applied in reverse order.
    public let virtualModes: [Int]
    public init(occupiedModes: [Int], virtualModes: [Int]) {
        self.occupiedModes = occupiedModes; self.virtualModes = virtualModes
    }
    public func validate(qubitCount: Int) throws {
        guard (1...2).contains(occupiedModes.count), occupiedModes.count == virtualModes.count,
              Set(occupiedModes).count == occupiedModes.count, Set(virtualModes).count == virtualModes.count,
              Set(occupiedModes).isDisjoint(with: Set(virtualModes)),
              (occupiedModes + virtualModes).allSatisfy({ $0 >= 0 && $0 < qubitCount }) else {
            throw VivoChemistryError.invalid("fermionic single/double excitation")
        }
        let occupiedSpins = occupiedModes.map { $0 & 1 }.sorted()
        let virtualSpins = virtualModes.map { $0 & 1 }.sorted()
        guard occupiedSpins == virtualSpins else {
            throw VivoChemistryError.invalid("excitation must conserve spin projection")
        }
    }
}

public enum VivoUCCSDPool {
    public static func referenceDeterminant(spatialOrbitals: Int, alphaElectrons: Int, betaElectrons: Int) throws -> UInt64 {
        guard spatialOrbitals > 0, spatialOrbitals <= 31,
              (0...spatialOrbitals).contains(alphaElectrons), (0...spatialOrbitals).contains(betaElectrons) else {
            throw VivoChemistryError.invalid("UCC reference occupations")
        }
        var determinant: UInt64 = 0
        for p in 0..<alphaElectrons { determinant |= UInt64(1) << (2*p) }
        for p in 0..<betaElectrons { determinant |= UInt64(1) << (2*p+1) }
        return determinant
    }

    public static func make(spatialOrbitals n: Int, alphaElectrons na: Int, betaElectrons nb: Int) throws -> [VivoFermionicExcitation] {
        let reference = try referenceDeterminant(spatialOrbitals: n, alphaElectrons: na, betaElectrons: nb)
        let qubits = 2*n
        let occupied = (0..<qubits).filter { reference & (UInt64(1) << $0) != 0 }
        let virtual = (0..<qubits).filter { reference & (UInt64(1) << $0) == 0 }
        var pool: [VivoFermionicExcitation] = []
        for i in occupied { for a in virtual where (i & 1) == (a & 1) {
            pool.append(.init(occupiedModes:[i],virtualModes:[a]))
        } }
        if occupied.count >= 2 && virtual.count >= 2 {
            for ii in 0..<(occupied.count-1) { for jj in (ii+1)..<occupied.count {
                for aa in 0..<(virtual.count-1) { for bb in (aa+1)..<virtual.count {
                    let o = [occupied[ii],occupied[jj]], v = [virtual[aa],virtual[bb]]
                    if o.map({$0&1}).sorted() == v.map({$0&1}).sorted() {
                        pool.append(.init(occupiedModes:o,virtualModes:v))
                    }
                } }
            } }
        }
        for excitation in pool { try excitation.validate(qubitCount:qubits) }
        return pool.sorted {
            if $0.occupiedModes.count != $1.occupiedModes.count { return $0.occupiedModes.count < $1.occupiedModes.count }
            if $0.occupiedModes != $1.occupiedModes { return $0.occupiedModes.lexicographicallyPrecedes($1.occupiedModes) }
            return $0.virtualModes.lexicographicallyPrecedes($1.virtualModes)
        }
    }
}

public enum VivoFermionicUnitary {
    /// Applies exp(theta * (tau - tau†)) exactly for one single/double excitation.
    /// The full product ansatz is therefore unitary and conserves N and Sz exactly.
    public static func applying(_ excitation: VivoFermionicExcitation, angle: Double,
                                to input: VivoStateVector, budget: VivoQuantumBudget = .init()) throws -> VivoStateVector {
        try excitation.validate(qubitCount:input.qubitCount)
        guard angle.isFinite else { throw VivoChemistryError.invalid("excitation angle") }
        var amplitudes = input.amplitudes
        let annihilate = excitation.occupiedModes.map { VivoFermionAction(mode:$0,creation:false) }
        let create = excitation.virtualModes.reversed().map { VivoFermionAction(mode:$0,creation:true) }
        let actions = annihilate + create
        let c = cos(angle), s = sin(angle)
        let occupiedMask = excitation.occupiedModes.reduce(UInt64(0)) { $0 | (UInt64(1)<<$1) }
        let virtualMask = excitation.virtualModes.reduce(UInt64(0)) { $0 | (UInt64(1)<<$1) }
        for source in amplitudes.indices {
            let bits = UInt64(source)
            guard bits & occupiedMask == occupiedMask, bits & virtualMask == 0 else { continue }
            guard let (destinationBits, sign) = vivoApplyFermions(bits, actions) else {
                throw VivoChemistryError.invalid("excitation applicability/sign inconsistency")
            }
            let destination = Int(destinationBits)
            let a = amplitudes[source], b = amplitudes[destination]
            amplitudes[source] = a*c - b*(sign*s)
            amplitudes[destination] = a*(sign*s) + b*c
        }
        return try .init(qubitCount:input.qubitCount,amplitudes:amplitudes,budget:budget)
    }

    public static func product(reference: VivoStateVector, excitations: [VivoFermionicExcitation],
                               parameters: [Double], budget: VivoQuantumBudget = .init()) throws -> VivoStateVector {
        guard excitations.count == parameters.count else { throw VivoChemistryError.invalid("ansatz parameter count") }
        var state = reference
        for (excitation,parameter) in zip(excitations,parameters) {
            state = try applying(excitation,angle:parameter,to:state,budget:budget)
        }
        return state
    }
}

public struct VivoVQEConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var maximumEnergyEvaluations: Int
    public var energyToleranceHartree: Double
    public var parameterTolerance: Double
    public var initialSearchRadius: Double
    public var minimumSearchRadius: Double
    public init(maximumIterations: Int = 128, maximumEnergyEvaluations: Int = 10_000,
                energyToleranceHartree: Double = 1e-10, parameterTolerance: Double = 1e-6,
                initialSearchRadius: Double = 0.25, minimumSearchRadius: Double = 1e-6) {
        self.maximumIterations=maximumIterations; self.maximumEnergyEvaluations=maximumEnergyEvaluations
        self.energyToleranceHartree=energyToleranceHartree; self.parameterTolerance=parameterTolerance
        self.initialSearchRadius=initialSearchRadius; self.minimumSearchRadius=minimumSearchRadius
    }
    public func validate() throws {
        guard maximumIterations > 0, maximumEnergyEvaluations > 0,
              [energyToleranceHartree,parameterTolerance,initialSearchRadius,minimumSearchRadius].allSatisfy({$0.isFinite && $0>0}),
              minimumSearchRadius <= initialSearchRadius else { throw VivoChemistryError.invalid("VQE optimizer configuration") }
    }
}
public struct VivoVQEResult: Codable, Sendable, Equatable {
    public let energyHartree: Double
    public let parameters: [Double]
    public let excitations: [VivoFermionicExcitation]
    public let evaluations: Int
    public let iterations: Int
    public let converged: Bool
    public let termination: String
    public let measurementGroups: Int
}

public enum VivoVQE {
    public static func optimize(hamiltonian: VivoPauliHamiltonian, reference: VivoStateVector,
                                excitations: [VivoFermionicExcitation], initialParameters: [Double]? = nil,
                                configuration cfg: VivoVQEConfiguration = .init(),
                                budget: VivoQuantumBudget = .init()) throws -> VivoVQEResult {
        try hamiltonian.validate(); try cfg.validate(); try budget.validate()
        guard hamiltonian.qubitCount == reference.qubitCount, excitations.count <= 4096 else {
            throw VivoChemistryError.resourceLimit("VQE ansatz/state contract")
        }
        for excitation in excitations { try excitation.validate(qubitCount:reference.qubitCount) }
        var parameters = initialParameters ?? [Double](repeating:0,count:excitations.count)
        guard parameters.count==excitations.count,parameters.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("VQE initial parameters") }
        var evaluations=0
        func energy(_ p:[Double]) throws -> Double {
            guard evaluations < cfg.maximumEnergyEvaluations else { throw VivoChemistryError.resourceLimit("VQE energy-evaluation budget") }
            evaluations += 1
            return try VivoFermionicUnitary.product(reference:reference,excitations:excitations,parameters:p,budget:budget).expectation(hamiltonian,budget:budget)
        }
        var value=try energy(parameters),radius=cfg.initialSearchRadius,iterations=0,converged=excitations.isEmpty,termination=excitations.isEmpty ? "empty-ansatz" : "iteration-limit"
        for iteration in 0..<cfg.maximumIterations where !converged {
            iterations=iteration+1
            let starting=value
            var largestMove=0.0
            for i in parameters.indices {
                var plus=parameters,minus=parameters; plus[i]+=radius; minus[i]-=radius
                let ep=try energy(plus), em=try energy(minus)
                if ep < value && ep <= em { parameters=plus; value=ep; largestMove=max(largestMove,radius) }
                else if em < value { parameters=minus; value=em; largestMove=max(largestMove,radius) }
                if evaluations >= cfg.maximumEnergyEvaluations { termination="evaluation-budget"; break }
            }
            if abs(starting-value)<=cfg.energyToleranceHartree {
                radius *= 0.5
                if radius < cfg.minimumSearchRadius || largestMove < cfg.parameterTolerance && radius <= cfg.parameterTolerance {
                    converged=true;termination="energy-and-parameter-tolerance";break
                }
            }
            if evaluations >= cfg.maximumEnergyEvaluations { break }
        }
        return .init(energyHartree:value,parameters:parameters,excitations:excitations,evaluations:evaluations,iterations:iterations,
                     converged:converged,termination:termination,measurementGroups:try VivoQubitWiseCommutingGrouping.group(hamiltonian).count)
    }
}

public struct VivoADAPTVQEConfiguration: Codable, Sendable, Equatable {
    public var maximumOperators: Int
    public var selectionGradientTolerance: Double
    public var finiteDifferenceStep: Double
    public var optimizer: VivoVQEConfiguration
    public init(maximumOperators: Int = 64, selectionGradientTolerance: Double = 1e-5,
                finiteDifferenceStep: Double = 1e-4, optimizer: VivoVQEConfiguration = .init()) {
        self.maximumOperators=maximumOperators; self.selectionGradientTolerance=selectionGradientTolerance
        self.finiteDifferenceStep=finiteDifferenceStep; self.optimizer=optimizer
    }
}
public enum VivoADAPTVQE {
    public static func solve(hamiltonian: VivoPauliHamiltonian, reference: VivoStateVector,
                             pool: [VivoFermionicExcitation], configuration cfg: VivoADAPTVQEConfiguration = .init(),
                             budget: VivoQuantumBudget = .init()) throws -> VivoVQEResult {
        guard cfg.maximumOperators>0,cfg.selectionGradientTolerance.isFinite,cfg.selectionGradientTolerance>0,
              cfg.finiteDifferenceStep.isFinite,cfg.finiteDifferenceStep>0 else { throw VivoChemistryError.invalid("ADAPT configuration") }
        try cfg.optimizer.validate()
        var selected:[VivoFermionicExcitation]=[],parameters:[Double]=[],remaining=pool
        var latest=try VivoVQE.optimize(hamiltonian:hamiltonian,reference:reference,excitations:selected,configuration:cfg.optimizer,budget:budget)
        for _ in 0..<min(cfg.maximumOperators,pool.count) {
            let current=try VivoFermionicUnitary.product(reference:reference,excitations:selected,parameters:parameters,budget:budget)
            var bestIndex:Int?,bestGradient=0.0
            for index in remaining.indices {
                let plus=try VivoFermionicUnitary.applying(remaining[index],angle:cfg.finiteDifferenceStep,to:current,budget:budget)
                let minus=try VivoFermionicUnitary.applying(remaining[index],angle:-cfg.finiteDifferenceStep,to:current,budget:budget)
                let gradient=(try plus.expectation(hamiltonian,budget:budget)-minus.expectation(hamiltonian,budget:budget))/(2*cfg.finiteDifferenceStep)
                if abs(gradient)>abs(bestGradient) { bestGradient=gradient;bestIndex=index }
            }
            guard let index=bestIndex,abs(bestGradient)>=cfg.selectionGradientTolerance else { return latest }
            selected.append(remaining.remove(at:index));parameters.append(0)
            latest=try VivoVQE.optimize(hamiltonian:hamiltonian,reference:reference,excitations:selected,initialParameters:parameters,configuration:cfg.optimizer,budget:budget)
            parameters=latest.parameters
        }
        return latest
    }
}
