import Foundation

public enum VivoReactiveAtomSource: Codable, Sendable, Equatable {
    /// Explicit atom correspondence is inherited from the ordered nuclei.
    case explicit(atomIndices: [Int])
    /// Endpoint atom array positions MUST already share the same mapping.
    /// No graph isomorphism, bond-distance inference or protonation guessing.
    case mappedEndpoints(reactant: VivoMolecularStructure, product: VivoMolecularStructure,
                         additionalAtomIndices: [Int])
}
public struct VivoReactiveAtomEvidence: Codable, Sendable, Equatable {
    public let atomIndex: Int
    public let atomIdentifier: String
    public let reasons: [String]
}
public enum VivoReactionAtomSeeds {
    public static func derive(source: VivoReactiveAtomSource, atomIdentifiers: [String],
                              system: VivoElectronicSystem) throws -> [VivoReactiveAtomEvidence] {
        let n = system.nuclei.count
        guard n > 0, atomIdentifiers.count == n, Set(atomIdentifiers).count == n,
              atomIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }) else {
            throw VivoChemistryError.invalid("reaction seed atom mapping")
        }
        var reasons: [Int:Set<String>] = [:]
        func add(_ indices: [Int], reason: String) throws {
            guard Set(indices).count == indices.count, indices.allSatisfy({ (0..<n).contains($0) }) else {
                throw VivoChemistryError.invalid("reaction seed references duplicate or unknown atoms")
            }
            for index in indices { reasons[index,default: []].insert(reason) }
        }
        switch source {
        case .explicit(let indices): try add(indices,reason: "declared reactive atom")
        case .mappedEndpoints(let reactant, let product, let additional):
            guard reactant.schemaVersion == VivoMolecularStructure.schemaVersion,
                  product.schemaVersion == VivoMolecularStructure.schemaVersion,
                  reactant.atoms.count == n, product.atoms.count == n,
                  reactant.bonds.count <= n*n, product.bonds.count <= n*n else {
                throw VivoChemistryError.invalid("mapped endpoint topology size or schema")
            }
            for i in 0..<n {
                let a = reactant.atoms[i], b = product.atoms[i]
                guard a.index == UInt32(i), b.index == UInt32(i), a.element.atomicNumber == b.element.atomicNumber,
                      Int(a.element.atomicNumber) == system.nuclei[i].atomicNumber,
                      a.isotopeMassNumber == b.isotopeMassNumber else {
                    throw VivoChemistryError.invalid("endpoint atom/isotope correspondence is not the declared nuclear mapping")
                }
                if a.formalCharge != b.formalCharge { reasons[i,default: []].insert("mapped formal-charge change") }
            }
            func bonds(_ structure: VivoMolecularStructure) throws -> [String:VivoMolecularBond] {
                var result: [String:VivoMolecularBond] = [:]
                for bond in structure.bonds {
                    let i = Int(bond.canonicalPair.0), j = Int(bond.canonicalPair.1)
                    guard (0..<n).contains(i), (0..<n).contains(j), i != j, bond.order != .unknown else {
                        throw VivoChemistryError.invalid("endpoint bond identity or unknown order")
                    }
                    let key = "\(i):\(j)"
                    guard result.updateValue(bond,forKey: key) == nil else { throw VivoChemistryError.invalid("duplicate endpoint bond") }
                }
                return result
            }
            let left = try bonds(reactant), right = try bonds(product)
            for key in Set(left.keys).union(right.keys).sorted() {
                let a = left[key], b = right[key]
                if a?.order != b?.order || a?.stereo != b?.stereo {
                    let bond = a ?? b!
                    try add([Int(bond.atomA),Int(bond.atomB)],reason: "mapped bond/order/stereochemistry change")
                }
            }
            try add(additional,reason: "declared additional reactive atom")
        }
        guard !reasons.isEmpty else {
            throw VivoChemistryError.invalid("no mapped reactive atoms; unchanged topology requires explicit reaction-center atoms")
        }
        return reasons.keys.sorted().map { .init(atomIndex: $0,atomIdentifier: atomIdentifiers[$0],reasons: reasons[$0]!.sorted()) }
    }
}
public struct VivoMolecularSpacePreparationRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-space-preparation/v1"
    public let schema: String
    public let identifier: String
    public let atomIdentifiers: [String]
    public let coordinateUnit: String
    public let snapshots: [VivoMolecularPathSnapshot]
    public let basis: VivoGaussianBasis
    public let reference: VivoSCFConfiguration
    public let reactiveAtoms: VivoReactiveAtomSource
    /// Nil means all AOs on reactive atoms, including any core/polarization AOs.
    /// Explicit shell indices can specify a physically intended valence target.
    public let targetShellIndices: [Int]?
    public let projection: VivoAtomicSpaceProjectionConfiguration
    public let discoveryPointIdentifiers: [String]
    public let confirmationPointIdentifiers: [String]
    public let target: VivoElectronicProfileTarget
    public let solver: VivoSpaceRefinementSolver
    public let maximumActiveOrbitals: Int
    public let maximumRounds: Int
    public let maximumPointEvaluations: Int
    public let maximumPreparationCalls: Int
    public let minimumTransportSingularValue: Double
    public let minimumStateOverlapSquared: Double
    public let budget: VivoChemistryBudget
    public init(identifier: String, atomIdentifiers: [String], coordinateUnit: String,
                snapshots: [VivoMolecularPathSnapshot], basis: VivoGaussianBasis,
                reactiveAtoms: VivoReactiveAtomSource, targetShellIndices: [Int]? = nil,
                reference: VivoSCFConfiguration = .init(), projection: VivoAtomicSpaceProjectionConfiguration = .init(),
                discoveryPointIdentifiers: [String], confirmationPointIdentifiers: [String],
                target: VivoElectronicProfileTarget, solver: VivoSpaceRefinementSolver = .directCI(configuration: .init()),
                maximumActiveOrbitals: Int = 31, maximumRounds: Int = 8, maximumPointEvaluations: Int = 512,
                maximumPreparationCalls: Int = 160, minimumTransportSingularValue: Double = 0.5,
                minimumStateOverlapSquared: Double = 0.1, budget: VivoChemistryBudget = .init()) {
        schema = Self.schema; self.identifier = identifier; self.atomIdentifiers = atomIdentifiers
        self.coordinateUnit = coordinateUnit; self.snapshots = snapshots; self.basis = basis
        self.reactiveAtoms = reactiveAtoms; self.targetShellIndices = targetShellIndices; self.reference = reference
        self.projection = projection; self.discoveryPointIdentifiers = discoveryPointIdentifiers
        self.confirmationPointIdentifiers = confirmationPointIdentifiers; self.target = target; self.solver = solver
        self.maximumActiveOrbitals = maximumActiveOrbitals; self.maximumRounds = maximumRounds
        self.maximumPointEvaluations = maximumPointEvaluations; self.maximumPreparationCalls = maximumPreparationCalls
        self.minimumTransportSingularValue = minimumTransportSingularValue; self.minimumStateOverlapSquared = minimumStateOverlapSquared
        self.budget = budget
    }
    public func validate() throws {
        try budget.validate(); try reference.validate(); try projection.validate()
        guard schema == Self.schema, !identifier.isEmpty, identifier.utf8.count <= 512,
              !coordinateUnit.isEmpty, coordinateUnit.utf8.count <= 128, (4...32).contains(snapshots.count),
              Set(snapshots.map(\.identifier)).count == snapshots.count,
              snapshots.allSatisfy({ !$0.identifier.isEmpty && $0.identifier.utf8.count <= 256 }),
              reference.reference == .restricted, (1...31).contains(maximumActiveOrbitals),
              (1...32).contains(maximumRounds), (1...100_000).contains(maximumPointEvaluations),
              (1...160).contains(maximumPreparationCalls), 5*snapshots.count-1 <= maximumPreparationCalls,
              [minimumTransportSingularValue,minimumStateOverlapSquared].allSatisfy({ $0.isFinite && $0 > 0 && $0 <= 1 }) else {
            throw VivoChemistryError.invalid("molecular space-preparation identity, scope, restricted reference or capacity")
        }
        let names = snapshots.map(\.identifier), discovery = Set(discoveryPointIdentifiers), holdout = Set(confirmationPointIdentifiers)
        guard discovery.count == discoveryPointIdentifiers.count, holdout.count == confirmationPointIdentifiers.count,
              discovery.count >= 3, !holdout.isEmpty, discovery.isDisjoint(with: holdout),
              discovery.union(holdout) == Set(names), discovery.contains(names[0]), discovery.contains(names.last!),
              discovery.contains(target.barrierPointIdentifier), names.dropFirst().dropLast().contains(target.barrierPointIdentifier),
              [target.maximumBarrierShiftHartree,target.maximumRelativeProfileShiftHartree].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoChemistryError.invalid("molecular discovery/confirmation split or electronic target")
        }
        try solver.statePolicy.validate()
        let first = snapshots[0].system
        try basis.validate(nucleusCount: first.nuclei.count)
        let evidence = try VivoReactionAtomSeeds.derive(source: reactiveAtoms,atomIdentifiers: atomIdentifiers,system: first)
        let atoms = Set(evidence.map(\.atomIndex))
        if let shells = targetShellIndices {
            guard !shells.isEmpty, Set(shells).count == shells.count,
                  shells.allSatisfy({ basis.shells.indices.contains($0) && atoms.contains(basis.shells[$0].nucleusIndex) }),
                  Set(shells.map { basis.shells[$0].nucleusIndex }) == atoms else {
                throw VivoChemistryError.invalid("explicit target shells must cover exactly the declared reactive atoms")
            }
        }
        var previous = -Double.infinity
        for snapshot in snapshots {
            try snapshot.system.validate()
            guard snapshot.coordinate.isFinite, snapshot.coordinate > previous,
                  snapshot.system.nuclei.map(\.atomicNumber) == first.nuclei.map(\.atomicNumber),
                  snapshot.system.nuclei.map(\.structureAtomIndex) == first.nuclei.map(\.structureAtomIndex),
                  snapshot.system.alphaElectrons == first.alphaElectrons,
                  snapshot.system.betaElectrons == first.betaElectrons,
                  first.alphaElectrons == first.betaElectrons,
                  snapshot.system.pointCharges == first.pointCharges else {
                throw VivoChemistryError.invalid("molecular refinement changes atom mapping, spin/charge or frozen environment")
            }
            previous = snapshot.coordinate
        }
        let n = try VivoGaussianIntegralEngine.expanded(system: first,basis: basis,budget: budget).count
        guard n <= 31 else { throw VivoChemistryError.resourceLimit("molecular refinement currently uses the finite 31-orbital determinant representation") }
        _ = try budget.elements([snapshots.count,n,n,n,n],simultaneousArrays: 24)
    }
}
public struct VivoMolecularSpacePreparationResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-space-preparation-result/v1"
    public let schema: String
    public let source: VivoMolecularSpacePreparationRequest
    public let reactiveAtomEvidence: [VivoReactiveAtomEvidence]
    public let projections: [VivoAtomicSpaceProjectionResult]
    public let refinementRequest: VivoPropertyDirectedSpaceRequest
    public let preparationCalls: Int
    public let interpretation: String
}
public enum VivoMolecularSpacePreparation {
    public static func prepare(_ request: VivoMolecularSpacePreparationRequest) throws -> VivoMolecularSpacePreparationResult {
        try request.validate()
        let integrals = try request.snapshots.map {
            try VivoGaussianIntegralEngine.compute(system: $0.system,basis: request.basis,budget: request.budget)
        }
        let references = try request.snapshots.indices.map {
            try VivoHartreeFock.solve(system: request.snapshots[$0].system,integrals: integrals[$0],
                configuration: request.reference,budget: request.budget)
        }
        return try prepare(request,integrals: integrals,references: references)
    }
    /// Workflow callers reuse the authoritative, content-addressed AO/HF stages.
    /// Supplied primitives must be bound to the exact system, basis and settings.
    public static func prepare(_ request: VivoMolecularSpacePreparationRequest,
                               integrals: [VivoAOIntegrals], references: [VivoHartreeFockResult]) throws -> VivoMolecularSpacePreparationResult {
        try request.validate()
        guard integrals.count == request.snapshots.count, references.count == request.snapshots.count else {
            throw VivoChemistryError.invalid("molecular preparation primitive count")
        }
        let evidence = try VivoReactionAtomSeeds.derive(source: request.reactiveAtoms,
            atomIdentifiers: request.atomIdentifiers,system: request.snapshots[0].system)
        let atoms = Set(evidence.map(\.atomIndex)), shells = request.targetShellIndices.map(Set.init)
        var projections: [VivoAtomicSpaceProjectionResult] = []
        for i in request.snapshots.indices {
            let snapshot = request.snapshots[i], ao = integrals[i], hf = references[i]
            guard ao.sourceSystem == snapshot.system, ao.sourceBasis == request.basis else {
                throw VivoChemistryError.invalid("molecular preparation AO source binding")
            }
            try VivoHartreeFock.validate(result: hf,system: snapshot.system,integrals: ao,
                configuration: request.reference,budget: request.budget)
            let targets = ao.orbitals.indices.filter { i in
                atoms.contains(ao.orbitals[i].nucleusIndex) && (shells?.contains(ao.orbitals[i].shellIndex) ?? true)
            }
            projections.append(try VivoAtomicSpaceProjection.select(overlap: ao.overlap,coefficients: hf.alphaCoefficients,
                occupied: snapshot.system.alphaElectrons,targetAOIndices: targets,
                configuration: request.projection,budget: request.budget))
        }
        let points = try VivoMolecularOrbitalFrames.prepare(snapshots: request.snapshots,basis: request.basis,
            reference: request.reference,integrals: integrals,references: references,
            rotations: projections.map(\.rotation),orbitalPrefix: "atomic-space",budget: request.budget)
        let initial = projections[0]
        // The first (discovery) geometry fixes thresholds, groups and candidates.
        // Holdout data cannot change them. All other orbitals remain candidates.
        let refinement = VivoPropertyDirectedSpaceRequest(identifier: request.identifier,points: points,
            transportGroups: initial.groups,initialSpace: initial.initialSpace,
            mandatoryActiveOrbitals: initial.initialSpace.active,candidateBlocks: initial.candidateBlocks,
            discoveryPointIdentifiers: request.discoveryPointIdentifiers,confirmationPointIdentifiers: request.confirmationPointIdentifiers,
            target: request.target,solver: request.solver,maximumRounds: request.maximumRounds,
            maximumActiveOrbitals: request.maximumActiveOrbitals,maximumPointEvaluations: request.maximumPointEvaluations,
            minimumTransportSingularValue: request.minimumTransportSingularValue,
            minimumStateOverlapSquared: request.minimumStateOverlapSquared,budget: request.budget)
        try refinement.validate()
        return .init(schema: VivoMolecularSpacePreparationResult.schema,source: request,reactiveAtomEvidence: evidence,
            projections: projections,refinementRequest: refinement,preparationCalls: 5*points.count-1,
            interpretation: "source-bound molecular preparation and conservative AO-projection seeding; preparation has a bounded count of per-call-budget numerical stages, while downstream refinement has its own aggregate work limit; not a mechanism search or a chemical-accuracy claim")
    }
    public static func validate(_ result: VivoMolecularSpacePreparationResult, request: VivoMolecularSpacePreparationRequest) throws {
        guard result == (try prepare(request)) else { throw VivoChemistryError.invalid("molecular space-preparation source reconstruction") }
    }
}
