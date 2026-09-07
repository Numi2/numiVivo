import Foundation

public struct VivoNCMCLambdaSchedule: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ncmc-lambda-schedule/v2"
    public var schema: String
    /// Strictly increasing values including 0 and 1.
    public var lambdas: [Double]
    /// Steps at each interval midpoint, between two thermodynamic perturbations.
    public var propagationStepsPerLambda: UInt64

    public init(lambdas: [Double], propagationStepsPerLambda: UInt64) {
        schema = Self.schema
        self.lambdas = lambdas
        self.propagationStepsPerLambda = propagationStepsPerLambda
    }

    public static func linear(intervals: Int, propagationStepsPerLambda: UInt64) throws -> Self {
        guard intervals >= 1, intervals <= 1_000_000 else {
            throw VivoChemistryError.invalid("NCMC lambda interval count")
        }
        return .init(lambdas: (0...intervals).map { Double($0) / Double(intervals) },
                     propagationStepsPerLambda: propagationStepsPerLambda)
    }

    /// The reverse traversal must mirror a nonuniform schedule, not reuse it.
    public func reversed() throws -> Self {
        try validate()
        return .init(lambdas: lambdas.reversed().map { 1 - $0 },
                     propagationStepsPerLambda: propagationStepsPerLambda)
    }

    public func validate() throws {
        guard schema == Self.schema,
              lambdas.count >= 2, lambdas.count <= 1_000_001,
              lambdas.first == 0, lambdas.last == 1,
              lambdas.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }),
              zip(lambdas, lambdas.dropFirst()).allSatisfy({ $0.0 < $0.1 }),
              propagationStepsPerLambda > 0,
              propagationStepsPerLambda <= 10_000_000 else {
            throw VivoChemistryError.invalid("NCMC lambda schedule")
        }
        let work = (lambdas.count - 1).multipliedReportingOverflow(by: Int(propagationStepsPerLambda))
        guard !work.overflow, work.partialValue <= 100_000_000 else {
            throw VivoChemistryError.resourceLimit("NCMC schedule propagation work")
        }
    }
}

public struct VivoNCMCPropagationResult: Codable, Sendable, Equatable {
    public let finalPhysicalState: VivoConstantPHPhysicalState
    /// For deterministic reversible volume-preserving propagation this is the
    /// FULL change in potential plus kinetic energy at fixed lambda. Symplectic
    /// or reversible does not imply energy-conserving: finite-step work is not zero.
    public let shadowWorkKJPerMol: Double
    /// v2 requires zero: stochastic heat/path-action accounting is not inferred
    /// from this scalar and must not be added again to heat-subtracted work.
    public let logReverseOverForwardPathProbability: Double
    public let steps: UInt64
    public let method: String

    public init(finalPhysicalState: VivoConstantPHPhysicalState,
                shadowWorkKJPerMol: Double,
                logReverseOverForwardPathProbability: Double = 0,
                steps: UInt64,
                method: String) {
        self.finalPhysicalState = finalPhysicalState
        self.shadowWorkKJPerMol = shadowWorkKJPerMol
        self.logReverseOverForwardPathProbability = logReverseOverForwardPathProbability
        self.steps = steps
        self.method = method
    }

    public func validate(initial: VivoConstantPHPhysicalState, expectedSteps: UInt64,
                         maximumAbsoluteWorkKJPerMol: Double) throws {
        try initial.validate(); try finalPhysicalState.validate()
        let (expectedStep, overflow) = initial.stepIndex.addingReportingOverflow(expectedSteps)
        guard !overflow, steps == expectedSteps, expectedSteps > 0,
              finalPhysicalState.positionsNM.count == initial.positionsNM.count,
              finalPhysicalState.velocitiesNMPerPS.count == initial.velocitiesNMPerPS.count,
              finalPhysicalState.stepIndex == expectedStep,
              finalPhysicalState.timePS > initial.timePS,
              finalPhysicalState.periodicCell == initial.periodicCell,
              shadowWorkKJPerMol.isFinite,
              abs(shadowWorkKJPerMol) <= maximumAbsoluteWorkKJPerMol,
              logReverseOverForwardPathProbability.isFinite,
              !method.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoChemistryError.invalid("NCMC propagation result")
        }
        guard logReverseOverForwardPathProbability == 0 else {
            throw VivoChemistryError.unsupported("NCMC v2 requires deterministic fixed-cell propagation; raw stochastic path ratios require explicit heat/action accounting")
        }
    }
}

/// Lambda-dependent complete physical Hamiltonian on one fixed particle/mass
/// manifold. The lambda interpolation itself is part of the fingerprinted model.
public struct VivoNCMCLambdaHamiltonian: Sendable {
    public let fingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let fromStateIdentifier: String
    public let toStateIdentifier: String
    public let completePotentialEnergyKJPerMol: @Sendable (VivoConstantPHPhysicalState, Double) async throws -> Double
    public let propagate: @Sendable (VivoConstantPHPhysicalState, Double, UInt64) async throws -> VivoNCMCPropagationResult

    public init(fingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                fromStateIdentifier: String,
                toStateIdentifier: String,
                completePotentialEnergyKJPerMol: @escaping @Sendable (VivoConstantPHPhysicalState, Double) async throws -> Double,
                propagate: @escaping @Sendable (VivoConstantPHPhysicalState, Double, UInt64) async throws -> VivoNCMCPropagationResult) {
        self.fingerprint = fingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.fromStateIdentifier = fromStateIdentifier
        self.toStateIdentifier = toStateIdentifier
        self.completePotentialEnergyKJPerMol = completePotentialEnergyKJPerMol
        self.propagate = propagate
    }
}

public struct VivoNCMCProtocolStep: Codable, Sendable, Equatable {
    public let index: Int
    public let lambdaBefore: Double
    public let lambdaAfter: Double
    public let prePerturbationStateFingerprint: VivoFingerprint
    public let propagationLambda: Double
    public let postPropagationStateFingerprint: VivoFingerprint
    public let prePropagationPerturbationWorkKJPerMol: Double
    public let postPropagationPerturbationWorkKJPerMol: Double
    public let postPropagationEnergyBeforeKJPerMol: Double
    public let postPropagationEnergyAfterKJPerMol: Double
    public let energyBeforeKJPerMol: Double
    public let energyAfterKJPerMol: Double
    public let perturbationWorkKJPerMol: Double
    public let propagationShadowWorkKJPerMol: Double
    public let logReverseOverForwardPathProbability: Double
    public let propagationSteps: UInt64
}

public struct VivoNCMCProtocolResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ncmc-protocol-result/v2"
    public let schema: String
    public let hamiltonianFingerprint: VivoFingerprint
    public let scheduleFingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let initialPhysicalStateFingerprint: VivoFingerprint
    public let finalPhysicalState: VivoConstantPHPhysicalState
    public let finalPhysicalStateFingerprint: VivoFingerprint
    public let perturbationWorkKJPerMol: Double
    public let shadowWorkKJPerMol: Double
    public let logReverseOverForwardPathProbability: Double
    public let steps: [VivoNCMCProtocolStep]
    public let interpretation: String

    public var totalNonequilibriumWorkKJPerMol: Double {
        perturbationWorkKJPerMol + shadowWorkKJPerMol
    }
}

/// Symmetric perturb-propagate-perturb NCMC. Each interval evaluates the
/// complete potential at its endpoints and midpoint, at the appropriate fixed
/// coordinates. A reversed mirrored schedule with reversed momenta follows the
/// reverse proposal. Deterministic segment work includes the kinetic change.
public enum VivoNCMCProtocolKernel {
    public static let interpretation = "Symmetric midpoint perturb-propagate-perturb NCMC on a fixed mass/cell manifold; every perturbation uses identical coordinates on its two sides. Deterministic fixed-lambda shadow work includes potential AND kinetic energy changes. Raw stochastic path ratios and unaccounted heat are rejected. Numerical reversibility remains precision- and convergence-qualified."

    public static func run(initial: VivoConstantPHPhysicalState,
                           hamiltonian: VivoNCMCLambdaHamiltonian,
                           schedule: VivoNCMCLambdaSchedule,
                           maximumAbsoluteWorkKJPerMol: Double = 1.0e8) async throws -> VivoNCMCProtocolResult {
        try initial.validate(); try schedule.validate()
        guard maximumAbsoluteWorkKJPerMol.isFinite, maximumAbsoluteWorkKJPerMol > 0,
              maximumAbsoluteWorkKJPerMol <= 1.0e12 else {
            throw VivoChemistryError.invalid("NCMC protocol work bound")
        }
        guard !hamiltonian.fromStateIdentifier.isEmpty, !hamiltonian.toStateIdentifier.isEmpty,
              hamiltonian.fromStateIdentifier != hamiltonian.toStateIdentifier else {
            throw VivoChemistryError.invalid("NCMC endpoint identifiers")
        }
        let scheduleID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(schedule))
        var physical = initial
        var perturbationWork = 0.0
        var shadowWork = 0.0
        var logPath = 0.0
        var trace: [VivoNCMCProtocolStep] = []
        trace.reserveCapacity(schedule.lambdas.count - 1)

        for index in 1..<schedule.lambdas.count {
            try Task.checkCancellation()
            let beforeLambda = schedule.lambdas[index - 1]
            let afterLambda = schedule.lambdas[index]
            let midpoint = beforeLambda + 0.5 * (afterLambda - beforeLambda)
            let perturbationState = physical
            let before = try await hamiltonian.completePotentialEnergyKJPerMol(perturbationState, beforeLambda)
            let atMidpoint = try await hamiltonian.completePotentialEnergyKJPerMol(perturbationState, midpoint)
            guard before.isFinite, atMidpoint.isFinite else {
                throw VivoChemistryError.convergence("nonfinite NCMC pre-propagation energy")
            }
            try Task.checkCancellation()
            let propagated = try await hamiltonian.propagate(perturbationState, midpoint, schedule.propagationStepsPerLambda)
            try propagated.validate(initial: perturbationState,
                                    expectedSteps: schedule.propagationStepsPerLambda,
                                    maximumAbsoluteWorkKJPerMol: maximumAbsoluteWorkKJPerMol)
            let post = propagated.finalPhysicalState
            let postMidpoint = try await hamiltonian.completePotentialEnergyKJPerMol(post, midpoint)
            let after = try await hamiltonian.completePotentialEnergyKJPerMol(post, afterLambda)
            guard [before, atMidpoint, postMidpoint, after].allSatisfy(\.isFinite) else {
                throw VivoChemistryError.convergence("nonfinite NCMC perturbation energy")
            }
            let firstWork = atMidpoint - before, lastWork = after - postMidpoint
            let delta = firstWork + lastWork
            perturbationWork += delta
            shadowWork += propagated.shadowWorkKJPerMol
            logPath += propagated.logReverseOverForwardPathProbability
            guard [firstWork, lastWork, perturbationWork, shadowWork, perturbationWork + shadowWork].allSatisfy({
                $0.isFinite && abs($0) <= maximumAbsoluteWorkKJPerMol
            }), logPath.isFinite else {
                throw VivoChemistryError.convergence("NCMC work accumulation exceeds numerical bounds")
            }
            trace.append(.init(index: index - 1,
                lambdaBefore: beforeLambda, lambdaAfter: afterLambda,
                prePerturbationStateFingerprint: try perturbationState.fingerprint(),
                propagationLambda: midpoint,
                postPropagationStateFingerprint: try post.fingerprint(),
                prePropagationPerturbationWorkKJPerMol: firstWork,
                postPropagationPerturbationWorkKJPerMol: lastWork,
                postPropagationEnergyBeforeKJPerMol: postMidpoint,
                postPropagationEnergyAfterKJPerMol: after,
                energyBeforeKJPerMol: before, energyAfterKJPerMol: atMidpoint,
                perturbationWorkKJPerMol: delta,
                propagationShadowWorkKJPerMol: propagated.shadowWorkKJPerMol,
                logReverseOverForwardPathProbability: propagated.logReverseOverForwardPathProbability,
                propagationSteps: propagated.steps))
            physical = post
        }
        try Task.checkCancellation()

        return .init(schema: VivoNCMCProtocolResult.schema,
            hamiltonianFingerprint: hamiltonian.fingerprint,
            scheduleFingerprint: scheduleID,
            physicalManifoldFingerprint: hamiltonian.physicalManifoldFingerprint,
            initialPhysicalStateFingerprint: try initial.fingerprint(),
            finalPhysicalState: physical,
            finalPhysicalStateFingerprint: try physical.fingerprint(),
            perturbationWorkKJPerMol: perturbationWork,
            shadowWorkKJPerMol: shadowWork,
            logReverseOverForwardPathProbability: logPath,
            steps: trace,
            interpretation: interpretation)
    }

    public static func switchEngine(fingerprint: VivoFingerprint,
                                    physicalManifoldFingerprint: VivoFingerprint,
                                    schedule: VivoNCMCLambdaSchedule,
                                    buildHamiltonian: @escaping @Sendable (
                                        VivoConstantPHExecutableState,
                                        VivoConstantPHExecutableState
                                    ) async throws -> VivoNCMCLambdaHamiltonian,
                                    maximumAbsoluteWorkKJPerMol: Double = 1.0e8) throws -> VivoConstantPHNCMCSwitchEngine {
        try schedule.validate()
        struct Identity: Codable {
            let schema: String; let implementation: VivoFingerprint; let manifold: VivoFingerprint
            let schedule: VivoNCMCLambdaSchedule; let workBound: Double
        }
        let engineID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            schema: "numivivo.org/ncmc-switch-engine/v2", implementation: fingerprint,
            manifold: physicalManifoldFingerprint, schedule: schedule, workBound: maximumAbsoluteWorkKJPerMol)))
        return VivoConstantPHNCMCSwitchEngine(fingerprint: engineID,
            physicalManifoldFingerprint: physicalManifoldFingerprint) { initial, from, to in
            let lambda = try await buildHamiltonian(from, to)
            guard lambda.physicalManifoldFingerprint == physicalManifoldFingerprint,
                  lambda.fromStateIdentifier == from.identifier,
                  lambda.toStateIdentifier == to.identifier else {
                throw VivoChemistryError.invalid("NCMC lambda Hamiltonian identity")
            }
            guard from.physicalManifoldFingerprint == physicalManifoldFingerprint,
                  to.physicalManifoldFingerprint == physicalManifoldFingerprint else {
                throw VivoChemistryError.invalid("NCMC endpoint physical manifold")
            }
            let directedSchedule = try from.identifier < to.identifier ? schedule : schedule.reversed()
            let result = try await run(initial: initial, hamiltonian: lambda,
                                       schedule: directedSchedule,
                                       maximumAbsoluteWorkKJPerMol: maximumAbsoluteWorkKJPerMol)
            return try .init(switchEngineFingerprint: engineID,
                physicalManifoldFingerprint: physicalManifoldFingerprint,
                fromStateIdentifier: from.identifier,
                toStateIdentifier: to.identifier,
                initialPhysicalState: initial,
                finalPhysicalState: result.finalPhysicalState,
                protocolWorkKJPerMol: result.perturbationWorkKJPerMol,
                shadowWorkKJPerMol: result.shadowWorkKJPerMol,
                logReverseOverForwardPathProbability: result.logReverseOverForwardPathProbability,
                switchingSteps: result.steps.reduce(UInt64(0)) { $0 + $1.propagationSteps },
                interpretation: result.interpretation)
        }
    }
}
