import Foundation

public struct VivoNCMCLambdaSchedule: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ncmc-lambda-schedule/v1"
    public var schema: String
    /// Strictly increasing values including 0 and 1.
    public var lambdas: [Double]
    /// Propagation steps after each perturbation to lambdas[i], i>0.
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
    /// Integrator/discretization work needed for exact path acceptance. A kernel
    /// that is analytically reversible and measure-preserving may report zero.
    public let shadowWorkKJPerMol: Double
    /// log(P_reverse/P_forward) for stochastic propagation at this fixed lambda.
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
        guard steps == expectedSteps,
              finalPhysicalState.positionsNM.count == initial.positionsNM.count,
              finalPhysicalState.velocitiesNMPerPS.count == initial.velocitiesNMPerPS.count,
              finalPhysicalState.stepIndex >= initial.stepIndex,
              finalPhysicalState.timePS >= initial.timePS,
              shadowWorkKJPerMol.isFinite,
              abs(shadowWorkKJPerMol) <= maximumAbsoluteWorkKJPerMol,
              logReverseOverForwardPathProbability.isFinite,
              !method.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoChemistryError.invalid("NCMC propagation result")
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
    public let energyBeforeKJPerMol: Double
    public let energyAfterKJPerMol: Double
    public let perturbationWorkKJPerMol: Double
    public let propagationShadowWorkKJPerMol: Double
    public let logReverseOverForwardPathProbability: Double
    public let propagationSteps: UInt64
}

public struct VivoNCMCProtocolResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ncmc-protocol-result/v1"
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

/// Exact bookkeeping for a perturb-propagate NCMC schedule. At each lambda
/// update, protocol work is U(lambda_new,x)-U(lambda_old,x) at the identical
/// pre-propagation state. Fixed-lambda propagation contributes independently
/// reported shadow work and stochastic path probability. No endpoint-only work
/// approximation is used.
public enum VivoNCMCProtocolKernel {
    public static let interpretation = "Perturb-propagate NCMC work accounting over a fingerprinted lambda Hamiltonian. Protocol work is accumulated at every instantaneous lambda change on the same coordinates; each fixed-lambda propagation segment supplies its own shadow work and reverse/forward path probability."

    public static func run(initial: VivoConstantPHPhysicalState,
                           hamiltonian: VivoNCMCLambdaHamiltonian,
                           schedule: VivoNCMCLambdaSchedule,
                           maximumAbsoluteWorkKJPerMol: Double = 1.0e8) async throws -> VivoNCMCProtocolResult {
        try initial.validate(); try schedule.validate()
        guard maximumAbsoluteWorkKJPerMol.isFinite, maximumAbsoluteWorkKJPerMol > 0,
              maximumAbsoluteWorkKJPerMol <= 1.0e12 else {
            throw VivoChemistryError.invalid("NCMC protocol work bound")
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
            // Both energies address one immutable pre-propagation snapshot.
            // Sequential evaluation also avoids overlapping access to a shared
            // resource-owning Metal evaluator and Swift 6 mutable captures.
            let perturbationState = physical
            let before = try await hamiltonian.completePotentialEnergyKJPerMol(perturbationState, beforeLambda)
            let after = try await hamiltonian.completePotentialEnergyKJPerMol(perturbationState, afterLambda)
            guard before.isFinite, after.isFinite else {
                throw VivoChemistryError.convergence("nonfinite NCMC perturbation energy")
            }
            let delta = after - before
            perturbationWork += delta
            guard perturbationWork.isFinite, abs(perturbationWork) <= maximumAbsoluteWorkKJPerMol else {
                throw VivoChemistryError.convergence("NCMC perturbation work overflow")
            }
            let preFingerprint = try perturbationState.fingerprint()
            let propagated = try await hamiltonian.propagate(perturbationState, afterLambda, schedule.propagationStepsPerLambda)
            try propagated.validate(initial: perturbationState,
                                    expectedSteps: schedule.propagationStepsPerLambda,
                                    maximumAbsoluteWorkKJPerMol: maximumAbsoluteWorkKJPerMol)
            shadowWork += propagated.shadowWorkKJPerMol
            logPath += propagated.logReverseOverForwardPathProbability
            guard shadowWork.isFinite, abs(shadowWork) <= maximumAbsoluteWorkKJPerMol,
                  logPath.isFinite else {
                throw VivoChemistryError.convergence("NCMC propagation bookkeeping overflow")
            }
            trace.append(.init(index: index - 1,
                lambdaBefore: beforeLambda, lambdaAfter: afterLambda,
                prePerturbationStateFingerprint: preFingerprint,
                energyBeforeKJPerMol: before, energyAfterKJPerMol: after,
                perturbationWorkKJPerMol: delta,
                propagationShadowWorkKJPerMol: propagated.shadowWorkKJPerMol,
                logReverseOverForwardPathProbability: propagated.logReverseOverForwardPathProbability,
                propagationSteps: propagated.steps))
            physical = propagated.finalPhysicalState
        }

        guard physical.positionsNM.count == initial.positionsNM.count,
              physical.velocitiesNMPerPS.count == initial.velocitiesNMPerPS.count else {
            throw VivoChemistryError.invalid("NCMC protocol changed particle manifold")
        }
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
                                    maximumAbsoluteWorkKJPerMol: Double = 1.0e8) -> VivoConstantPHNCMCSwitchEngine {
        VivoConstantPHNCMCSwitchEngine(fingerprint: fingerprint,
            physicalManifoldFingerprint: physicalManifoldFingerprint) { initial, from, to in
            let lambda = try await buildHamiltonian(from, to)
            guard lambda.physicalManifoldFingerprint == physicalManifoldFingerprint,
                  lambda.fromStateIdentifier == from.identifier,
                  lambda.toStateIdentifier == to.identifier else {
                throw VivoChemistryError.invalid("NCMC lambda Hamiltonian identity")
            }
            let result = try await run(initial: initial, hamiltonian: lambda,
                                       schedule: schedule,
                                       maximumAbsoluteWorkKJPerMol: maximumAbsoluteWorkKJPerMol)
            return try .init(switchEngineFingerprint: fingerprint,
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
