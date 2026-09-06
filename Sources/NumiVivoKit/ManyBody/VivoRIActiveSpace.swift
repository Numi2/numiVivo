import Foundation

/// Indices refer to the canonical RI-RHF molecular orbitals, not AO indices.
/// Every occupied orbital must be active or explicitly frozen; excluded virtual
/// orbitals are a declared model truncation, never an automatic fallback.
public struct VivoRIActiveSpace: Codable, Sendable, Equatable {
    public var activeOrbitals: [Int]
    public var frozenDoublyOccupiedOrbitals: [Int]
    public var identifier: String
    public init(activeOrbitals: [Int], frozenDoublyOccupiedOrbitals: [Int] = [], identifier: String) {
        self.activeOrbitals = activeOrbitals
        self.frozenDoublyOccupiedOrbitals = frozenDoublyOccupiedOrbitals
        self.identifier = identifier
    }
    public func validate(orbitalCount: Int, occupied: Int) throws {
        let selected = activeOrbitals + frozenDoublyOccupiedOrbitals
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              orbitalCount > 0, occupied >= 0, occupied <= orbitalCount,
              !activeOrbitals.isEmpty, Set(selected).count == selected.count,
              selected.allSatisfy({ $0 >= 0 && $0 < orbitalCount }),
              frozenDoublyOccupiedOrbitals.allSatisfy({ $0 < occupied }),
              Set(0..<occupied).isSubset(of: Set(selected)) else {
            throw VivoChemistryError.invalid("RI active/frozen partition must be unique, disjoint and contain all occupied orbitals")
        }
    }
}

public extension VivoEmbeddedHamiltonian {
    /// Prepare the EXISTING solver boundary directly from RI factors. The parent
    /// AO and full-MO four-index tensors are never formed. Only the explicitly
    /// selected active-space m^4 tensor is materialized for the current solvers.
    ///
    /// Frozen-core J/K and its physical scalar are evaluated in the AO basis.
    /// This avoids storing core-core/core-active transformed factor blocks.
    static func fromRI(_ reference: VivoRIHartreeFockResult, space: VivoRIActiveSpace,
                       budget: VivoChemistryBudget = .init()) throws -> Self {
        _ = try VivoFactorizedHartreeFock.validateRestricted(reference,budget:budget)
        let one = reference.source.oneElectron, factors = reference.source.factors
        let hf = reference.scf, n = one.count, active = space.activeOrbitals
        let core = space.frozenDoublyOccupiedOrbitals, m = active.count
        try space.validate(orbitalCount:n,occupied:hf.alphaElectrons)
        let count = try budget.elements([m,m,m,m],simultaneousArrays:4)
        _ = try budget.elements([n,n],simultaneousArrays:16)
        _ = try budget.operatorApplications([m,m,m,m,factors.rank],label:"RI active Hamiltonian")
        var ca = VivoQMMatrix(n,m), coreDensity = VivoQMMatrix(n,n)
        for p in 0..<n {
            for (j,orbital) in active.enumerated() { ca[p,j] = hf.alphaCoefficients[p,orbital] }
            for q in 0..<n {
                for i in core { coreDensity[p,q] += hf.alphaCoefficients[p,i] * hf.alphaCoefficients[q,i] }
            }
        }
        var effective = one.coreHamiltonian, scalar = one.constantEnergyHartree
        if !core.isEmpty {
            let jk = try factors.coulombExchange(density:coreDensity,budget:budget)
            let interaction = try jk.coulomb.scaled(2).adding(jk.exchange,scale:-1)
            effective = try effective.adding(interaction)
            let coreEnergyOperator = try one.coreHamiltonian.scaled(2).adding(interaction)
            scalar += zip(coreDensity.values,coreEnergyOperator.values).reduce(0.0) { $0+$1.0*$1.1 }
        }
        let projected = try factors.transformed(by:ca,budget:budget)
        var g = [Double](repeating:0,count:count)
        // Fill one representative and its full real chemist symmetry orbit.
        // Redundant assignments on diagonal pairs are harmless and deterministic.
        for p in 0..<m { for q in 0...p { for r in 0..<m { for s in 0...r {
            if VivoCoulombFactors.pair(r,s) > VivoCoulombFactors.pair(p,q) { continue }
            let value = projected.eri(p,q,r,s)
            for (a,b,c,d) in [(p,q,r,s),(q,p,r,s),(p,q,s,r),(q,p,s,r),
                               (r,s,p,q),(s,r,p,q),(r,s,q,p),(s,r,q,p)] {
                g[((a*m+b)*m+c)*m+d] = value
            }
        } } } }
        let h = Self(orbitalIdentifiers:active.map { "ri-mo-\($0)" },
                     alphaElectrons:hf.alphaElectrons-core.count,
                     betaElectrons:hf.betaElectrons-core.count,
                     oneElectron:try effective.congruence(ca),twoElectron:g,
                     constantEnergyHartree:scalar,
                     energyReference:"RI physical electronic Hamiltonian with explicit frozen-core scalar",
                     provenance:["integralConvention":"real-spatial-chemist", "precision":"fp64",
                                 "parentBasis":one.basis.identifier,
                                 "auxiliaryBasis":reference.source.auxiliaryBasis.identifier,
                                 "factorMethod":factors.method, "factorRank":String(factors.rank),
                                 "activeSpaceIdentifier":space.identifier,
                                 "activeParentOrbitals":active.map(String.init).joined(separator:","),
                                 "frozenDoublyOccupiedOrbitals":core.map(String.init).joined(separator:","),
                                 "approximation":"density fitting plus declared frozen-core/virtual-space truncation"])
        try h.validate(budget:budget)
        return h
    }
}

public extension VivoAdvancedChemistryOperations {
    /// Artifact dependencies retain the complete RI reference, parent molecular
    /// system and explicit active-space configuration. Existing many-body and
    /// ECC/DMET operations consume this output without a new solver interface.
    static func riHamiltonian(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        func read<T: Decodable>(_ type: T.Type, _ name: String, _ inputs: [String:Data]) throws -> T {
            guard let data = inputs[name] else { throw VivoChemistryError.invalid("missing RI Hamiltonian slot \(name)") }
            return try VivoCanonicalJSON.decode(type,from:data)
        }
        return .init(identifier:"vivo.native.ri-active-hamiltonian",version:"1",implementationFingerprint:id,
                     outputs:[.init(name:"hamiltonian",kind:"vivo.embedded-hamiltonian")],
                     execute:{ cfg,inputs,budget in
            guard Set(inputs.keys) == Set(["reference"]) else { throw VivoChemistryError.invalid("RI Hamiltonian input slots") }
            let space = try VivoCanonicalJSON.decode(VivoRIActiveSpace.self,from:VivoCanonicalJSON.encode(cfg))
            let reference = try read(VivoRIHartreeFockResult.self,"reference",inputs)
            return ["hamiltonian":try VivoCanonicalJSON.encode(VivoEmbeddedHamiltonian.fromRI(reference,space:space,budget:budget))]
        }, validateOutputs:{ cfg,inputs,outputs,budget in
            guard Set(inputs.keys) == Set(["reference"]), Set(outputs.keys) == Set(["hamiltonian"]) else {
                throw VivoChemistryError.invalid("RI Hamiltonian output slots")
            }
            let space = try VivoCanonicalJSON.decode(VivoRIActiveSpace.self,from:VivoCanonicalJSON.encode(cfg))
            let reference = try read(VivoRIHartreeFockResult.self,"reference",inputs)
            let actual = try read(VivoEmbeddedHamiltonian.self,"hamiltonian",outputs)
            let rebuilt = try VivoEmbeddedHamiltonian.fromRI(reference,space:space,budget:budget)
            try actual.validate(budget:budget)
            guard actual.orbitalIdentifiers == rebuilt.orbitalIdentifiers,
                  actual.alphaElectrons == rebuilt.alphaElectrons, actual.betaElectrons == rebuilt.betaElectrons,
                  actual.provenance == rebuilt.provenance, actual.energyReference == rebuilt.energyReference,
                  abs(actual.constantEnergyHartree-rebuilt.constantEnergyHartree) < 1e-9,
                  try actual.oneElectron.adding(rebuilt.oneElectron,scale:-1).frobeniusNorm < 1e-9,
                  zip(actual.twoElectron,rebuilt.twoElectron).allSatisfy({ abs($0.0-$0.1) < 1e-9 }) else {
                throw VivoChemistryError.invalid("RI Hamiltonian does not reconstruct from its reference and active space")
            }
        })
    }
}
