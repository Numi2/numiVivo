import Foundation

/// A complete orbital frame and a CAS projector define a subspace of the same
/// global fixed-(Nalpha,Nbeta) Fock sector. Different fragments may overlap.
/// They are combined as amplitudes, NEVER by adding their density matrices.
public struct VivoFockFragment: Codable, Sendable, Equatable {
    public let identifier: String
    public let partition: VivoActiveSpace
    /// Columns of the fragment frame in the common orthonormal spatial frame.
    public let orbitalRotation: VivoQMMatrix?
    public init(identifier: String, partition: VivoActiveSpace, orbitalRotation: VivoQMMatrix? = nil) {
        self.identifier = identifier; self.partition = partition; self.orbitalRotation = orbitalRotation
    }
}
public struct VivoVariationalSpaceConfiguration: Codable, Sendable, Equatable {
    public var maximumDimension: Int
    public var linearDependenceTolerance: Double
    public var projectedResidualToleranceHartree: Double
    public init(maximumDimension: Int = 256, linearDependenceTolerance: Double = 1e-10,
                projectedResidualToleranceHartree: Double = 1e-9) {
        self.maximumDimension = maximumDimension; self.linearDependenceTolerance = linearDependenceTolerance
        self.projectedResidualToleranceHartree = projectedResidualToleranceHartree
    }
    public func validate() throws {
        guard (1...4096).contains(maximumDimension), linearDependenceTolerance.isFinite,
              (1e-14...1e-5).contains(linearDependenceTolerance), projectedResidualToleranceHartree.isFinite,
              projectedResidualToleranceHartree > 0, projectedResidualToleranceHartree <= 1e-5 else {
            throw VivoChemistryError.invalid("variational Fock-space dimension, rank or eigen-residual tolerance")
        }
    }
}
public struct VivoVariationalCIResult: Codable, Sendable, Equatable {
    public let energyHartree: Double
    public let state: VivoCIState
    public let variationalDimension: Int
    public let fullSectorDimension: Int
    public let projectedResidualHartree: Double
    /// This is physical residual outside the retained space, NOT an energy-error
    /// bound. A converged projected eigenpair need not be a converged full CI.
    public let externalResidualHartree: Double
    public let fullResidualHartree: Double
    public let method: String
}

/// Bounded reference implementation with matrix-free Slater-Condon action. The
/// vectors still inhabit the complete determinant sector. Report that cost,
/// even when the eigenproblem is much smaller. This is not scalable DMRG/DMET.
struct VivoFockSpace {
    let orbitalCount: Int
    let alphaElectrons: Int
    let betaElectrons: Int
    let determinants: [UInt64]
    let configuration: VivoVariationalSpaceConfiguration
    let budget: VivoChemistryBudget
    private(set) var columns: [[Double]] = []
    private(set) var dependentColumns = 0
    var dimension: Int { columns.count }
    init(hamiltonian: VivoEmbeddedHamiltonian, configuration: VivoVariationalSpaceConfiguration,
         budget: VivoChemistryBudget) throws {
        try configuration.validate(); try hamiltonian.validate(budget: budget)
        orbitalCount = hamiltonian.orbitalCount; alphaElectrons = hamiltonian.alphaElectrons
        betaElectrons = hamiltonian.betaElectrons; self.configuration = configuration; self.budget = budget
        determinants = try VivoDirectCI.determinants(n: orbitalCount, na: alphaElectrons, nb: betaElectrons, budget: budget)
        _ = try budget.elements([determinants.count, min(configuration.maximumDimension, determinants.count)], simultaneousArrays: 12)
    }
    static func dot(_ a: [Double], _ b: [Double]) -> Double { zip(a,b).reduce(0) { $0+$1.0*$1.1 } }
    static func norm(_ v: [Double]) -> Double { v.reduce(0) { hypot($0,$1) } }
    func removeProjection(_ vector: [Double]) -> [Double] {
        var q = vector
        // Reorthogonalization is essential when nearly duplicate orbital frames
        // or residual directions enter through overlapping fragments.
        for _ in 0..<2 { for b in columns {
            let p = Self.dot(b,q); for i in q.indices { q[i] -= p*b[i] }
        } }
        return q
    }
    @discardableResult mutating func append(_ vector: [Double]) throws -> Bool {
        guard vector.count == determinants.count, vector.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("variational basis vector shape or finite values")
        }
        let original = Self.norm(vector)
        guard original > 0 else { dependentColumns += 1; return false }
        var q = removeProjection(vector.map { $0/original })
        let length = Self.norm(q)
        if length <= configuration.linearDependenceTolerance { dependentColumns += 1; return false }
        guard columns.count < configuration.maximumDimension, columns.count < determinants.count else {
            throw VivoChemistryError.resourceLimit("variational space capacity; independent directions cannot be silently discarded")
        }
        q = q.map { $0/length }
        if let largest = q.indices.max(by: { abs(q[$0]) < abs(q[$1]) }), q[largest] < 0 { q = q.map { -$0 } }
        columns.append(q); return true
    }
    mutating func add(_ fragment: VivoFockFragment, hamiltonian h: VivoEmbeddedHamiltonian) throws {
        try fragment.partition.validate(for: h, budget: budget)
        guard !fragment.identifier.isEmpty, fragment.identifier.utf8.count <= 1024,
              fragment.partition.frozenOrbitals.isEmpty else {
            throw VivoChemistryError.invalid("global fragment identity; frozen optimization masks have no meaning in a fixed subspace")
        }
        let p = fragment.partition, n = orbitalCount
        let rotation = try fragment.orbitalRotation ?? .identity(n)
        guard rotation.rows == n, rotation.columns == n, rotation.values.allSatisfy(\.isFinite),
              try rotation.transposed.multiplied(by: rotation).adding(.identity(n),scale:-1).frobeniusNorm < 1e-9 else {
            throw VivoChemistryError.invalid("global fragment orbital frame is not orthogonal")
        }
        let active = try VivoDirectCI.determinants(n:p.active.count,na:alphaElectrons-p.doublyOccupiedCore.count,
            nb:betaElectrons-p.doublyOccupiedCore.count,budget:budget)
        let index = Dictionary(uniqueKeysWithValues: determinants.enumerated().map { ($0.element,$0.offset) })
        let isIdentity = try rotation.adding(.identity(n),scale:-1).frobeniusNorm < 1e-14
        for local in active {
            var bits: UInt64 = 0
            for i in p.doublyOccupiedCore { bits |= UInt64(3) << (2*i) }
            for i in p.active.indices { for spin in 0..<2 where local & (UInt64(1) << (2*i+spin)) != 0 {
                bits |= UInt64(1) << (2*p.active[i]+spin)
            } }
            var vector = [Double](repeating:0,count:determinants.count)
            if isIdentity { vector[index[bits]!] = 1 }
            else {
                // A basis determinant's sign is immaterial to its span. The
                // exact exterior-power rotation handles all relative phases.
                let localState = VivoCIState(orbitalCount:n,alphaElectrons:alphaElectrons,betaElectrons:betaElectrons,
                    determinants:[bits],coefficients:[1])
                let lifted = try VivoCIOrbitalFrame.rotated(localState,by:rotation.transposed,budget:budget)
                for (d,c) in zip(lifted.determinants,lifted.coefficients) { vector[index[d]!] = c }
            }
            try append(vector)
        }
    }
    func solve(_ h: VivoEmbeddedHamiltonian, work: inout Int) throws -> VivoVariationalCIResult {
        guard dimension > 0, h.orbitalCount == orbitalCount, h.alphaElectrons == alphaElectrons,
              h.betaElectrons == betaElectrons else { throw VivoChemistryError.invalid("empty or mismatched global Fock space") }
        let action = try VivoDirectHamiltonian(h,determinants:determinants,budget:budget)
        let images = try columns.map { try action.apply($0,work:&work) }
        var projected = VivoQMMatrix(dimension,dimension)
        var gramError = 0.0
        for i in 0..<dimension { for j in 0...i {
            let a = Self.dot(columns[i],images[j]), b = Self.dot(columns[j],images[i])
            guard abs(a-b) <= 1e-9*max(1,abs(a),abs(b)) else { throw VivoChemistryError.invalid("projected Hamiltonian lost Hermiticity") }
            let value = 0.5*(a+b); projected[i,j] = value; projected[j,i] = value
            gramError = max(gramError,abs(Self.dot(columns[i],columns[j])-(i == j ? 1:0)))
        } }
        guard gramError < 1e-9 else { throw VivoChemistryError.convergence("global subspace orthogonality") }
        let spectrum = try VivoQMDenseAlgebra.symmetricEigen(projected,tolerance:1e-13,maximumSweeps:256)
        var vector = [Double](repeating:0,count:determinants.count)
        for j in 0..<dimension { for i in vector.indices { vector[i] += columns[j][i]*spectrum.vectors[j,0] } }
        let image = try action.apply(vector,work:&work)
        let electronic = Self.dot(vector,image)
        let residual = zip(image,vector).map { $0-electronic*$1 }
        let inside = Self.norm(columns.map { Self.dot($0,residual) })
        let outside = Self.norm(removeProjection(residual))
        guard inside <= configuration.projectedResidualToleranceHartree else {
            throw VivoChemistryError.convergence("variational CI projected residual did not converge")
        }
        let state = VivoCIState(orbitalCount:orbitalCount,alphaElectrons:alphaElectrons,betaElectrons:betaElectrons,
            determinants:determinants,coefficients:vector)
        try state.validate(budget:budget)
        return .init(energyHartree:electronic+h.constantEnergyHartree,state:state,variationalDimension:dimension,
            fullSectorDimension:determinants.count,projectedResidualHartree:inside,externalResidualHartree:outside,
            fullResidualHartree:Self.norm(residual),method:"global variational CI in an explicit orthonormal Fock subspace; not full-CI convergence unless the full residual also converges")
    }
    /// Add physical Hamiltonian residuals at ALL geometries to one shared space.
    /// No full-CI energies/vectors enter the expansion or rank decision.
    mutating func enrich(hamiltonians: [VivoEmbeddedHamiltonian], states: [VivoVariationalCIResult],
                         work: inout Int) throws -> Int {
        guard hamiltonians.count == states.count else { throw VivoChemistryError.invalid("path residual cardinality") }
        var residuals: [[Double]] = []
        for (h,state) in zip(hamiltonians,states) {
            guard state.state.determinants == determinants else { throw VivoChemistryError.invalid("path residual determinant convention") }
            let action = try VivoDirectHamiltonian(h,determinants:determinants,budget:budget)
            let v = state.state.coefficients
            let image = try action.apply(v,work:&work)
            let energy = state.energyHartree-h.constantEnergyHartree
            let r = removeProjection(zip(image,v).map { $0-energy*$1 })
            if Self.norm(r) > configuration.projectedResidualToleranceHartree { residuals.append(r) }
        }
        // Collect from the same old space before appending: path ordering must
        // not change which physical residuals were requested.
        let previous = dimension
        for r in residuals { try append(r) }
        return dimension-previous
    }
    var matrix: VivoQMMatrix {
        var result = VivoQMMatrix(determinants.count,dimension)
        for j in columns.indices { for i in determinants.indices { result[i,j] = columns[j][i] } }
        return result
    }
}
public enum VivoCIOneParticleDensity {
    /// A normalized global CI state supplies the density, including interference
    /// between overlapping fragment amplitudes. No positive-density projection
    /// or post-hoc electron-number repair is performed.
    public static func spatial(_ state: VivoCIState,budget:VivoChemistryBudget = .init()) throws -> VivoQMMatrix {
        try state.validate(budget:budget)
        let n = state.orbitalCount
        _ = try budget.elements([n,n],simultaneousArrays:4)
        let (work,overflow) = (2*n*n).multipliedReportingOverflow(by:state.determinants.count)
        guard !overflow, work <= budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("global density contraction work") }
        let index = Dictionary(uniqueKeysWithValues:state.determinants.enumerated().map { ($0.element,$0.offset) })
        var d = VivoQMMatrix(n,n)
        for p in 0..<n { for q in 0..<n { for spin in 0..<2 {
            d[p,q] += VivoCIDensityMatrices.expectation(state,index,[.init(mode:2*q+spin,creation:false),.init(mode:2*p+spin,creation:true)])
        } } }
        let spectrum = try VivoQMDenseAlgebra.symmetricEigen(d,tolerance:1e-13)
        guard spectrum.values.allSatisfy({ $0 >= -1e-9 && $0 <= 2+1e-9 }),
              abs(spectrum.values.reduce(0,+)-Double(state.alphaElectrons+state.betaElectrons)) < 1e-8 else {
            throw VivoChemistryError.convergence("global CI density trace or occupation bounds")
        }
        return d
    }
}
