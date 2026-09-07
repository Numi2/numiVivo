import Foundation

/// Shared connected Slater-Condon action for full CI, selected CI and residuals.
/// The physical scalar is restored by the caller, not included in diagonal.
struct VivoDirectHamiltonian {
    let h: VivoEmbeddedHamiltonian
    let determinants: [UInt64]
    let index: [UInt64: Int]
    let occupied: [[Int]]
    let virtual: [[Int]]
    let diagonal: [Double]
    let budget: VivoChemistryBudget
    init(_ h: VivoEmbeddedHamiltonian, determinants: [UInt64], budget: VivoChemistryBudget) throws {
        try h.validate(budget: budget)
        let n = h.orbitalCount
        guard (1...31).contains(n), !determinants.isEmpty, determinants.count <= budget.maximumDeterminants else {
            throw VivoChemistryError.invalid("direct Hamiltonian sector")
        }
        _ = try budget.elements([determinants.count,2*n], simultaneousArrays: 3)
        var seed = [Double](repeating: 0, count: determinants.count); seed[0] = 1
        try VivoCIState(orbitalCount: n, alphaElectrons: h.alphaElectrons, betaElectrons: h.betaElectrons,
                        determinants: determinants, coefficients: seed).validate(budget: budget)
        self.h = h; self.determinants = determinants; self.budget = budget
        index = Dictionary(uniqueKeysWithValues: determinants.enumerated().map { ($0.element,$0.offset) })
        occupied = determinants.map { det in (0..<(2*n)).filter { det & (UInt64(1) << $0) != 0 } }
        virtual = determinants.map { det in (0..<(2*n)).filter { det & (UInt64(1) << $0) == 0 } }
        diagonal = occupied.map { Self.diagonalElement(h, occupied: $0) }
        guard diagonal.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("CI diagonal overflow") }
    }
    private static func diagonalElement(_ h: VivoEmbeddedHamiltonian, occupied: [Int]) -> Double {
        var energy = 0.0
        for p in occupied {
            energy += h.oneElectron[p/2,p/2]
            for q in occupied {
                energy += 0.5*(h.eri(p/2,p/2,q/2,q/2)-(p%2 == q%2 ? h.eri(p/2,q/2,q/2,p/2) : 0))
            }
        }
        return energy
    }
    func diagonalElement(_ determinant: UInt64) -> Double {
        Self.diagonalElement(h, occupied: (0..<(2*h.orbitalCount)).filter { determinant & (UInt64(1) << $0) != 0 })
    }
    private func integral(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        (p%2 == r%2 && q%2 == s%2 ? h.eri(p/2,r/2,q/2,s/2) : 0)
        - (p%2 == s%2 && q%2 == r%2 ? h.eri(p/2,s/2,q/2,r/2) : 0)
    }
    /// Includes all connected destinations, including those absent from `index`.
    /// No coefficient or integral threshold is used. Work is charged before emit.
    func connections(from k: Int, work: inout Int,
                     emit: (UInt64, Double) throws -> Void) throws {
        guard determinants.indices.contains(k) else { throw VivoChemistryError.invalid("CI source index") }
        func charge() throws {
            guard work < budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("direct Hamiltonian aggregate work") }
            work += 1
        }
        let det = determinants[k], occ = occupied[k], vir = virtual[k]
        try charge(); try emit(det,diagonal[k])
        for i in occ { for a in vir where i%2 == a%2 {
            try charge()
            let action = vivoApplyFermions(det,[.init(mode: i,creation: false),.init(mode: a,creation: true)])!
            var value = h.oneElectron[a/2,i/2]
            for j in occ where j != i { value += integral(a,j,i,j) }
            try emit(action.0,action.1*value)
        } }
        for ii in occ.indices { for jj in (ii+1)..<occ.count {
            let i = occ[ii], j = occ[jj]
            for aa in vir.indices { for bb in (aa+1)..<vir.count {
                let a = vir[aa], b = vir[bb]
                if i%2+j%2 != a%2+b%2 { continue }
                try charge()
                let action = vivoApplyFermions(det,[.init(mode: i,creation: false),.init(mode: j,creation: false),
                    .init(mode: b,creation: true),.init(mode: a,creation: true)])!
                try emit(action.0,action.1*integral(a,b,i,j))
            } }
        } }
    }
    func apply(_ vector: [Double], work: inout Int) throws -> [Double] {
        guard vector.count == determinants.count, vector.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("Hamiltonian action vector")
        }
        var out = [Double](repeating: 0, count: vector.count)
        for k in determinants.indices where vector[k] != 0 {
            try connections(from: k,work: &work) { destination,value in
                if let i = index[destination] { out[i] += vector[k]*value }
            }
        }
        guard out.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("Hamiltonian action overflow") }
        return out
    }
}
public extension VivoDirectCI {
    /// Complete fixed-spin-sector residual, including components absent from
    /// the sparse input state. A stored residual field is never trusted.
    static func residualNorm(hamiltonian h: VivoEmbeddedHamiltonian, state: VivoCIState,
                             energyHartree: Double, budget: VivoChemistryBudget = .init()) throws -> Double {
        try state.validate(budget: budget); try h.validate(budget: budget)
        guard state.orbitalCount == h.orbitalCount, state.alphaElectrons == h.alphaElectrons,
              state.betaElectrons == h.betaElectrons, energyHartree.isFinite else {
            throw VivoChemistryError.invalid("CI residual state/Hamiltonian binding")
        }
        let basis = try determinants(n: h.orbitalCount,na: h.alphaElectrons,nb: h.betaElectrons,budget: budget)
        let action = try VivoDirectHamiltonian(h,determinants: basis,budget: budget)
        var vector = [Double](repeating: 0,count: basis.count)
        for (det,c) in zip(state.determinants,state.coefficients) {
            guard let i = action.index[det] else { throw VivoChemistryError.invalid("state is outside its declared sector") }
            vector[i] = c
        }
        var work = 0
        let applied = try action.apply(vector,work: &work), electronic = energyHartree-h.constantEnergyHartree
        return zip(applied,vector).reduce(0.0) { hypot($0,$1.0-electronic*$1.1) }
    }
}
