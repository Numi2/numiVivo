import Foundation

public struct VivoECCSolventClosureConfiguration: Codable, Sendable, Equatable {
    public var maximumOuterIterations: Int
    public var projectorTolerance: Double
    public var projectedStationarityToleranceHartree: Double
    public var maximumAggregateECCFrameEvaluations: Int
    public init(maximumOuterIterations: Int = 32, projectorTolerance: Double = 1e-7,
                projectedStationarityToleranceHartree: Double = 1e-8,
                maximumAggregateECCFrameEvaluations: Int = 100000) {
        self.maximumOuterIterations = maximumOuterIterations
        self.projectorTolerance = projectorTolerance
        self.projectedStationarityToleranceHartree = projectedStationarityToleranceHartree
        self.maximumAggregateECCFrameEvaluations = maximumAggregateECCFrameEvaluations
    }
    public func validate() throws {
        guard (2...256).contains(maximumOuterIterations), projectorTolerance.isFinite, projectorTolerance > 0,
              projectedStationarityToleranceHartree.isFinite, projectedStationarityToleranceHartree > 0,
              projectedStationarityToleranceHartree <= 1e-5,
              (2...10_000_000).contains(maximumAggregateECCFrameEvaluations) else {
            throw VivoChemistryError.invalid("ECC/solvent outer iteration, projector, stationarity or aggregate work contract")
        }
    }
}

public struct VivoECCSolventClosureRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ecc-solvent-closure/v1"
    public let schema: String
    public let molecule: VivoCorrelatedSolventRequest
    public let embedding: VivoECCDMETConfiguration
    /// One list per ECC cluster. Columns are in that cluster's complete frame
    /// and are explicitly declared doubly occupied inactive orbitals.
    public let inactiveDoublyOccupiedColumns: [[Int]]
    public let globalSpace: VivoVariationalSpaceConfiguration
    public let configuration: VivoECCSolventClosureConfiguration
    public init(molecule: VivoCorrelatedSolventRequest, embedding: VivoECCDMETConfiguration,
                inactiveDoublyOccupiedColumns: [[Int]],
                globalSpace: VivoVariationalSpaceConfiguration = .init(),
                configuration: VivoECCSolventClosureConfiguration = .init()) {
        schema = Self.schema; self.molecule = molecule; self.embedding = embedding
        self.inactiveDoublyOccupiedColumns = inactiveDoublyOccupiedColumns
        self.globalSpace = globalSpace; self.configuration = configuration
    }
    public func validate() throws {
        try molecule.validate(); try globalSpace.validate(); try configuration.validate()
        guard schema == Self.schema, molecule.configuration.partition == nil,
              inactiveDoublyOccupiedColumns.count == embedding.fragments.count,
              configuration.projectedStationarityToleranceHartree <= max(1e-8, molecule.configuration.energyToleranceHartree * 100),
              configuration.maximumAggregateECCFrameEvaluations >= configuration.maximumOuterIterations else {
            throw VivoChemistryError.invalid("ECC/solvent molecular partition, cluster assumptions or tolerance binding")
        }
    }
}

public struct VivoECCSolventClosureIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let gasEnergyHartree: Double
    public let eccDemocraticEnergyHartree: Double
    public let densityResidual: Double
    public let potentialResidualHartree: Double
    public let projectorResidual: Double
    public let energyChangeHartree: Double?
    public let eccMomentResidual: Double
    public let eccElectronDefect: Double
    public let projectedResidualHartree: Double
    public let externalResidualHartree: Double
    public let variationalDimension: Int
    public let fullSectorDimension: Int
    public let eccFrameEvaluations: Int
}

public struct VivoECCSolventClosureResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ecc-solvent-closure-result/v1"
    public let schema: String
    public let request: VivoECCSolventClosureRequest
    public let coefficients: VivoQMMatrix
    public let ecc: VivoECCDMETResult
    public let globalSubspace: VivoQMMatrix
    public let state: VivoCIState
    public let densityAO: VivoQMMatrix
    public let occupations: [Double]
    public let equilibriumField: VivoSmoothCPCMResult
    public let gasEnergyHartree: Double
    public let energyHartree: Double
    public let densityResidual: Double
    public let potentialResidualHartree: Double
    public let projectorResidual: Double
    public let projectedResidualHartree: Double
    public let externalResidualHartree: Double
    public let aggregateECCFrameEvaluations: Int
    public let hamiltonianOperatorApplications: Int
    public let history: [VivoECCSolventClosureIteration]
    public let method: String
}

/// Alternates the existing ECC-DMET orbital/bath/moment cycle with a single
/// globally normalized variational CI state and smooth C-PCM. The solvent sees
/// one N-representable density including interference between overlapping ECC
/// projectors. The reported physical energy is <H_gas> + G_pol[D].
///
/// This deliberately does NOT relabel the democratic ECC-DMET energy as a
/// variational solvent functional. ECC chooses a correlated projector family;
/// the global CI union supplies the density/energy authority for solvent closure.
public enum VivoECCSolventClosure {
    public static let method = "ECC-DMET-guided coherent global Fock union + variational CI + equilibrium smooth C-PCM; one N-representable solvent density; democratic ECC energy retained only as a diagnostic; not a proof of complete ECC-DMET/PCM or a kinetic-rate model"

    private struct Workspace {
        let coefficients: VivoQMMatrix
        let gas: VivoEmbeddedHamiltonian
        let pcm: VivoSmoothCPCMOperator
        init(_ molecule: VivoCorrelatedSolventRequest) throws {
            let ao = try VivoGaussianIntegralEngine.compute(system: molecule.system, basis: molecule.basis, budget: molecule.budget)
            if let specified = molecule.coefficients { coefficients = specified }
            else {
                let eig = try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
                guard eig.values.first! >= 1e-8 else { throw VivoChemistryError.invalid("ECC/solvent AO linear dependence") }
                var c = eig.vectors
                for i in 0..<ao.count { for j in 0..<ao.count { c[i,j] /= sqrt(eig.values[j]) } }
                coefficients = c
            }
            guard coefficients.rows == ao.count, coefficients.columns == ao.count,
                  coefficients.values.allSatisfy(\.isFinite),
                  try ao.overlap.congruence(coefficients).adding(.identity(ao.count), scale: -1).frobeniusNorm < 1e-8 else {
                throw VivoChemistryError.invalid("ECC/solvent requires one complete orthonormal molecular frame")
            }
            gas = try VivoEmbeddedHamiltonian.fromAO(ao, coefficients: coefficients,
                alphaElectrons: molecule.system.alphaElectrons, betaElectrons: molecule.system.betaElectrons,
                orbitalIdentifiers: (0..<ao.count).map { "ecc-solvent-orbital-\($0)" },
                energyReference: "physical gas electronic Hamiltonian; AO scalar once; equilibrium solvent polarization accounted separately",
                budget: molecule.budget)
            pcm = try .init(system: molecule.system, basis: molecule.basis,
                configuration: molecule.solvent, budget: molecule.budget)
        }
    }

    private static func effective(_ gas: VivoEmbeddedHamiltonian, potential: VivoQMMatrix) throws -> VivoEmbeddedHamiltonian {
        .init(orbitalIdentifiers: gas.orbitalIdentifiers, alphaElectrons: gas.alphaElectrons,
              betaElectrons: gas.betaElectrons, oneElectron: try gas.oneElectron.adding(potential),
              twoElectron: gas.twoElectron, constantEnergyHartree: gas.constantEnergyHartree,
              energyReference: gas.energyReference + "; equilibrium solvent trial field")
    }

    private static func fragments(_ result: VivoECCDMETResult, request: VivoECCSolventClosureRequest,
                                  h: VivoEmbeddedHamiltonian) throws -> [VivoFockFragment] {
        guard result.frame.clusters.count == request.inactiveDoublyOccupiedColumns.count else {
            throw VivoChemistryError.invalid("ECC/solvent cluster/core assumption cardinality")
        }
        var output: [VivoFockFragment] = []
        for (i, cluster) in result.frame.clusters.enumerated() {
            let k = cluster.coefficients.columns, core = request.inactiveDoublyOccupiedColumns[i]
            guard core.allSatisfy({ $0 >= k && $0 < h.orbitalCount }), Set(core).count == core.count,
                  h.alphaElectrons - core.count == cluster.fragment.clusterAlphaElectrons,
                  h.betaElectrons - core.count == cluster.fragment.clusterBetaElectrons else {
                throw VivoChemistryError.invalid("ECC/solvent inactive occupations are absent, fractional or inconsistent with cluster electron counts")
            }
            let rotation = try result.orbitalRotation.multiplied(by: cluster.completeFrame)
            output.append(.init(identifier: cluster.fragment.identifier,
                partition: .init(doublyOccupiedCore: core, active: Array(0..<k)), orbitalRotation: rotation))
        }
        return output
    }

    private static func buildSpace(ecc: VivoECCDMETResult, request: VivoECCSolventClosureRequest,
                                   h: VivoEmbeddedHamiltonian) throws -> VivoFockSpace {
        var space = try VivoFockSpace(hamiltonian: h, configuration: request.globalSpace, budget: request.molecule.budget)
        for fragment in try fragments(ecc, request: request, h: h) { try space.add(fragment, hamiltonian: h) }
        guard space.dimension > 0 else { throw VivoChemistryError.convergence("ECC/solvent produced an empty global projector union") }
        return space
    }

    private static func projector(_ space: VivoFockSpace) throws -> VivoQMMatrix {
        let b = space.matrix
        return try b.multiplied(by: b.transposed)
    }

    private static func aoDensity(_ p: VivoQMMatrix, coefficients c: VivoQMMatrix) throws -> VivoQMMatrix {
        try c.multiplied(by: p).multiplied(by: c.transposed)
    }

    public static func run(_ request: VivoECCSolventClosureRequest) throws -> VivoECCSolventClosureResult {
        try request.validate()
        let w = try Workspace(request.molecule), cfg = request.molecule.configuration, budget = request.molecule.budget
        try request.embedding.validate(hamiltonian: w.gas, budget: budget)

        // Gas-phase ECC supplies the first physically declared projector family;
        // no FCI/solvent answer is used to initialize the outer closure.
        var ecc = try VivoECCDMET.solve(w.gas, configuration: request.embedding, budget: budget)
        guard ecc.converged else { throw VivoChemistryError.convergence("initial gas ECC-DMET: \(ecc.termination)") }
        var aggregateFrames = ecc.frameEvaluations
        guard aggregateFrames <= request.configuration.maximumAggregateECCFrameEvaluations else {
            throw VivoChemistryError.resourceLimit("ECC/solvent initial frame-evaluation budget")
        }
        var space = try buildSpace(ecc: ecc, request: request, h: w.gas), work = 0
        var solved = try space.solve(w.gas, work: &work)
        var densityMO = try VivoCIOneParticleDensity.spatial(solved.state, budget: budget)
        var densityAO = try aoDensity(densityMO, coefficients: w.coefficients)
        var oldProjector = try projector(space)
        var previousEnergy: Double?, history: [VivoECCSolventClosureIteration] = []
        var previousRotation = ecc.orbitalRotation
        var previousPotential = ecc.correlationPotentialHartree

        for iteration in 1...request.configuration.maximumOuterIterations {
            let inputDensity = densityAO
            let oldField = try w.pcm.evaluate(totalDensity: inputDensity)
            let oldPotential = try oldField.reactionPotentialMatrix.congruence(w.coefficients)
            let hEff = try effective(w.gas, potential: oldPotential)
            ecc = try VivoECCDMET.solve(hEff, configuration: request.embedding,
                initialRotation: previousRotation,
                initialPotentialHartree: previousPotential.isEmpty ? nil : previousPotential,
                budget: budget)
            guard ecc.converged else { throw VivoChemistryError.convergence("ECC/solvent outer ECC cycle \(iteration): \(ecc.termination)") }
            aggregateFrames += ecc.frameEvaluations
            guard aggregateFrames <= request.configuration.maximumAggregateECCFrameEvaluations else {
                throw VivoChemistryError.resourceLimit("ECC/solvent aggregate frame-evaluation budget")
            }
            space = try buildSpace(ecc: ecc, request: request, h: hEff)
            let newProjector = try projector(space)
            let projectorResidual = try newProjector.adding(oldProjector, scale: -1).frobeniusNorm
            solved = try space.solve(hEff, work: &work)
            let nextMO = try VivoCIOneParticleDensity.spatial(solved.state, budget: budget)
            let nextAO = try aoDensity(nextMO, coefficients: w.coefficients)
            let field = try w.pcm.evaluate(totalDensity: nextAO)
            let potential = try field.reactionPotentialMatrix.congruence(w.coefficients)
            let densityResidual = try nextMO.adding(densityMO, scale: -1).frobeniusNorm
            let potentialResidual = try potential.adding(oldPotential, scale: -1).frobeniusNorm
            let gasRDM = try VivoCIDensityMatrices.compute(solved.state, budget: budget)
            let gasEnergy = try gasRDM.energy(of: w.gas, budget: budget)
            let energy = gasEnergy + field.polarizationEnergyHartree
            let energyChange = previousEnergy.map { abs(energy - $0) }
            let momentResidual = ecc.matching?.momentResiduals.reduce(0.0) { hypot($0,$1) } ?? 0
            history.append(.init(iteration: iteration, energyHartree: energy, gasEnergyHartree: gasEnergy,
                eccDemocraticEnergyHartree: ecc.energyHartree, densityResidual: densityResidual,
                potentialResidualHartree: potentialResidual, projectorResidual: projectorResidual,
                energyChangeHartree: energyChange, eccMomentResidual: momentResidual,
                eccElectronDefect: ecc.frame.electronCountDefect,
                projectedResidualHartree: solved.projectedResidualHartree,
                externalResidualHartree: solved.externalResidualHartree,
                variationalDimension: solved.variationalDimension, fullSectorDimension: solved.fullSectorDimension,
                eccFrameEvaluations: ecc.frameEvaluations))

            if let energyChange,
               densityResidual <= cfg.densityTolerance,
               potentialResidual <= cfg.potentialToleranceHartree,
               projectorResidual <= request.configuration.projectorTolerance,
               energyChange <= cfg.energyToleranceHartree,
               solved.projectedResidualHartree <= request.configuration.projectedStationarityToleranceHartree,
               abs(ecc.frame.electronCountDefect) <= request.embedding.electronTolerance,
               momentResidual <= request.embedding.matching.momentTolerance {
                // Re-solve in the returned field/projector and verify closure,
                // rather than accepting a state generated by the previous field.
                let finalH = try effective(w.gas, potential: potential)
                let finalECC = try VivoECCDMET.solve(finalH, configuration: request.embedding,
                    initialRotation: ecc.orbitalRotation,
                    initialPotentialHartree: ecc.correlationPotentialHartree.isEmpty ? nil : ecc.correlationPotentialHartree,
                    budget: budget)
                guard finalECC.converged else { throw VivoChemistryError.convergence("ECC/solvent returned-field ECC closure") }
                aggregateFrames += finalECC.frameEvaluations
                guard aggregateFrames <= request.configuration.maximumAggregateECCFrameEvaluations else {
                    throw VivoChemistryError.resourceLimit("ECC/solvent returned-field frame budget")
                }
                let finalSpace = try buildSpace(ecc: finalECC, request: request, h: finalH)
                let finalProjector = try projector(finalSpace)
                guard try finalProjector.adding(newProjector, scale: -1).frobeniusNorm <= 2 * request.configuration.projectorTolerance else {
                    throw VivoChemistryError.convergence("ECC/solvent projector changes in returned reaction field")
                }
                let closure = try finalSpace.solve(finalH, work: &work)
                let closureMO = try VivoCIOneParticleDensity.spatial(closure.state, budget: budget)
                guard try closureMO.adding(nextMO, scale: -1).frobeniusNorm <= 2 * cfg.densityTolerance else {
                    throw VivoChemistryError.convergence("ECC/solvent returned-field density closure")
                }
                let closureAO = try aoDensity(closureMO, coefficients: w.coefficients)
                let closureField = try w.pcm.evaluate(totalDensity: closureAO)
                let closurePotential = try closureField.reactionPotentialMatrix.congruence(w.coefficients)
                guard try closurePotential.adding(potential, scale: -1).frobeniusNorm <= 2 * cfg.potentialToleranceHartree else {
                    throw VivoChemistryError.convergence("ECC/solvent returned-field potential closure")
                }
                let closureRDM = try VivoCIDensityMatrices.compute(closure.state, budget: budget)
                let closureGas = try closureRDM.energy(of: w.gas, budget: budget)
                let closureEnergy = closureGas + closureField.polarizationEnergyHartree
                guard abs(closureEnergy - energy) <= max(2 * cfg.energyToleranceHartree, 1e-10) else {
                    throw VivoChemistryError.convergence("ECC/solvent returned-field energy closure")
                }
                let occupations = try VivoQMDenseAlgebra.symmetricEigen(closureMO, tolerance: 1e-13).values
                return .init(schema: VivoECCSolventClosureResult.schema, request: request,
                    coefficients: w.coefficients, ecc: finalECC, globalSubspace: finalSpace.matrix,
                    state: closure.state, densityAO: closureAO, occupations: occupations,
                    equilibriumField: closureField, gasEnergyHartree: closureGas, energyHartree: closureEnergy,
                    densityResidual: densityResidual, potentialResidualHartree: potentialResidual,
                    projectorResidual: projectorResidual,
                    projectedResidualHartree: closure.projectedResidualHartree,
                    externalResidualHartree: closure.externalResidualHartree,
                    aggregateECCFrameEvaluations: aggregateFrames, hamiltonianOperatorApplications: work,
                    history: history, method: method)
            }

            previousEnergy = energy
            densityMO = try nextMO.scaled(1 - cfg.damping).adding(densityMO, scale: cfg.damping)
            densityAO = try aoDensity(densityMO, coefficients: w.coefficients)
            oldProjector = newProjector
            previousRotation = ecc.orbitalRotation
            previousPotential = ecc.correlationPotentialHartree
        }
        throw VivoChemistryError.convergence("ECC-DMET/global-density/C-PCM closure did not reach joint stationarity")
    }

    public static func validate(_ result: VivoECCSolventClosureResult, request: VivoECCSolventClosureRequest) throws {
        guard result.schema == VivoECCSolventClosureResult.schema, result.request == request,
              result.method == method, result.history.count >= 2,
              result.history.count <= request.configuration.maximumOuterIterations,
              result.aggregateECCFrameEvaluations <= request.configuration.maximumAggregateECCFrameEvaluations,
              result.densityResidual <= request.molecule.configuration.densityTolerance,
              result.potentialResidualHartree <= request.molecule.configuration.potentialToleranceHartree,
              result.projectorResidual <= request.configuration.projectorTolerance,
              result.projectedResidualHartree <= request.configuration.projectedStationarityToleranceHartree,
              abs(result.ecc.frame.electronCountDefect) <= request.embedding.electronTolerance else {
            throw VivoChemistryError.invalid("ECC/solvent stored result convergence, work or method binding")
        }
        let rebuilt = try run(request)
        guard rebuilt == result else {
            throw VivoChemistryError.invalid("ECC/solvent wavefunction, projector, density, field, energy or ECC closure differs on reconstruction")
        }
    }

    /// Full-bath H2 is deliberately an integration/control problem: it should
    /// reproduce the full-CI correlated solvent result while exercising the
    /// actual ECC -> global projector -> solvent feedback loop.
    public static func hydrogenControl() -> VivoECCSolventClosureRequest {
        let system = VivoElectronicSystem(nuclei: [
            .init(atomicNumber: 1, positionBohr: .init(0,0,-0.7)),
            .init(atomicNumber: 1, positionBohr: .init(0,0,0.7))], alphaElectrons: 1, betaElectrons: 1)
        let molecule = VivoCorrelatedSolventRequest(system: system,
            basis: .hydrogenSTO3G(nucleusIndices: [0,1]), solvent: .init(dielectricConstant: 4, angularPoints: 50),
            configuration: .init(maximumIterations: 64, densityTolerance: 1e-9,
                                 potentialToleranceHartree: 1e-9, energyToleranceHartree: 1e-11,
                                 ciResidualTolerance: 1e-12))
        let ecc = VivoECCDMETConfiguration(mode: .singleFragment,
            fragments: [.init(identifier: "H-fragment", orbitals: [0], maximumBathOrbitals: 1,
                              clusterAlphaElectrons: 1, clusterBetaElectrons: 1)],
            matching: .init(referenceMethod: .fci), bathSelection: .init(minimumBathOrbitals: 1))
        return .init(molecule: molecule, embedding: ecc, inactiveDoublyOccupiedColumns: [[]],
                     globalSpace: .init(maximumDimension: 16),
                     configuration: .init(maximumOuterIterations: 32, projectorTolerance: 1e-8,
                                          projectedStationarityToleranceHartree: 1e-9,
                                          maximumAggregateECCFrameEvaluations: 4096))
    }
}
