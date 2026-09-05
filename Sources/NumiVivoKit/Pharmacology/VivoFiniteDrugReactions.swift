import Foundation

/// A local, fixed-volume reaction operator, not a second circulation/transport
/// engine. Concentrations and material ledgers are in M; time is in seconds.
/// Bound drug removed with its target becomes the explicitly inactive metabolite.
public struct VivoFiniteDrugModel: Codable, Equatable, Sendable {
    public let context: VivoKineticContext
    public let association: VivoKineticParameter
    public let dissociation: VivoKineticParameter
    public let inactivation: VivoKineticParameter
    public let targetTurnover: VivoKineticParameter
    public let baselineTarget: VivoKineticParameter
    public let freeDrugClearance: VivoKineticParameter
    public let metabolism: VivoKineticParameter
    public let metaboliteClearance: VivoKineticParameter
    public let maximumFreeDrugM: Double
    public init(context: VivoKineticContext, association: VivoKineticParameter,
                dissociation: VivoKineticParameter, inactivation: VivoKineticParameter,
                targetTurnover: VivoKineticParameter, baselineTarget: VivoKineticParameter,
                freeDrugClearance: VivoKineticParameter, metabolism: VivoKineticParameter,
                metaboliteClearance: VivoKineticParameter, maximumFreeDrugM: Double) {
        self.context = context; self.association = association; self.dissociation = dissociation
        self.inactivation = inactivation; self.targetTurnover = targetTurnover; self.baselineTarget = baselineTarget
        self.freeDrugClearance = freeDrugClearance; self.metabolism = metabolism
        self.metaboliteClearance = metaboliteClearance; self.maximumFreeDrugM = maximumFreeDrugM
    }
    public func validate() throws {
        try context.validate()
        try association.validate(unit: .perMolarSecond, label: "finite-pool association")
        try baselineTarget.validate(unit: .molar, label: "finite-pool target baseline", positive: true)
        for p in [dissociation, inactivation, targetTurnover, freeDrugClearance, metabolism, metaboliteClearance] {
            try p.validate(unit: .perSecond, label: "finite-pool rate")
        }
        guard maximumFreeDrugM.isFinite, maximumFreeDrugM > 0,
              (dissociation.value + inactivation.value).isFinite,
              (freeDrugClearance.value + metabolism.value).isFinite,
              (targetTurnover.value * baselineTarget.value).isFinite else {
            throw VivoKineticsError.invalid("finite-pool domain or derived rate")
        }
    }
}

public struct VivoFiniteDrugState: Codable, Equatable, Sendable {
    public var freeDrugM: Double
    public var freeTargetM: Double
    public var reversibleComplexM: Double
    public var covalentComplexM: Double
    public var inactiveMetaboliteM: Double
    public var eliminatedDrugEquivalentM: Double
    public var synthesizedTargetM: Double
    public var removedTargetM: Double
    public init(freeDrugM: Double, freeTargetM: Double, reversibleComplexM: Double = 0,
                covalentComplexM: Double = 0, inactiveMetaboliteM: Double = 0,
                eliminatedDrugEquivalentM: Double = 0, synthesizedTargetM: Double = 0,
                removedTargetM: Double = 0) {
        self.freeDrugM = freeDrugM; self.freeTargetM = freeTargetM; self.reversibleComplexM = reversibleComplexM
        self.covalentComplexM = covalentComplexM; self.inactiveMetaboliteM = inactiveMetaboliteM
        self.eliminatedDrugEquivalentM = eliminatedDrugEquivalentM
        self.synthesizedTargetM = synthesizedTargetM; self.removedTargetM = removedTargetM
    }
    public var values: [Double] { [freeDrugM, freeTargetM, reversibleComplexM, covalentComplexM,
        inactiveMetaboliteM, eliminatedDrugEquivalentM, synthesizedTargetM, removedTargetM] }
    public var totalTargetM: Double { freeTargetM + reversibleComplexM + covalentComplexM }
    public var drugBalanceM: Double { freeDrugM + reversibleComplexM + covalentComplexM + inactiveMetaboliteM + eliminatedDrugEquivalentM }
    public var targetBalanceM: Double { totalTargetM + removedTargetM - synthesizedTargetM }
    public var occupancy: Double? { totalTargetM > 0 ? (reversibleComplexM + covalentComplexM) / totalTargetM : nil }
    public func validate() throws {
        guard values.allSatisfy({ $0.isFinite && $0 >= 0 }), totalTargetM.isFinite,
              drugBalanceM.isFinite, targetBalanceM.isFinite else {
            throw VivoKineticsError.invalid("finite-pool state or material ledger")
        }
    }
}

public struct VivoFiniteDrugPolicy: Codable, Equatable, Sendable {
    public let relativeTolerance: Double
    public let absoluteToleranceM: Double
    public let balanceRelativeTolerance: Double
    public let maximumStepSeconds: Double
    public let minimumStepSeconds: Double
    public let maximumAttempts: Int
    public init(relativeTolerance: Double = 1e-6, absoluteToleranceM: Double = 1e-13,
                balanceRelativeTolerance: Double = 1e-10, maximumStepSeconds: Double = 1,
                minimumStepSeconds: Double = 1e-12, maximumAttempts: Int = 100_000) {
        self.relativeTolerance = relativeTolerance; self.absoluteToleranceM = absoluteToleranceM
        self.balanceRelativeTolerance = balanceRelativeTolerance; self.maximumStepSeconds = maximumStepSeconds
        self.minimumStepSeconds = minimumStepSeconds; self.maximumAttempts = maximumAttempts
    }
    public func validate() throws {
        guard relativeTolerance.isFinite, (1e-12...1e-2).contains(relativeTolerance),
              absoluteToleranceM.isFinite, absoluteToleranceM > 0,
              balanceRelativeTolerance.isFinite, (1e-14...1e-4).contains(balanceRelativeTolerance),
              maximumStepSeconds.isFinite, minimumStepSeconds.isFinite, minimumStepSeconds > 0,
              maximumStepSeconds >= minimumStepSeconds, (1...1_000_000).contains(maximumAttempts) else {
            throw VivoKineticsError.invalid("finite-pool numerical policy")
        }
    }
}

public struct VivoFiniteDrugStep: Codable, Equatable, Sendable {
    public let initial: VivoFiniteDrugState
    public let candidate: VivoFiniteDrugState
    public let durationSeconds: Double
    public let acceptedSubsteps: Int
    public let rejectedAttempts: Int
    public let maximumAcceptedLocalError: Double
    public let maximumRelativeMaterialError: Double
}

/// Positive symmetric splitting of exact local reactions with deterministic
/// step doubling. The caller owns commit/rollback and any dosing/transport.
/// No stochastic redraw, concentration clipping, or conservation repair occurs.
public enum VivoFiniteDrugReactions {
    public static func advance(_ initial: VivoFiniteDrugState, model: VivoFiniteDrugModel,
                               durationSeconds: Double, policy: VivoFiniteDrugPolicy = .init()) throws -> VivoFiniteDrugStep {
        try model.validate(); try initial.validate(); try policy.validate()
        guard durationSeconds.isFinite, durationSeconds >= 0, initial.freeDrugM <= model.maximumFreeDrugM else {
            throw VivoKineticsError.invalid("finite-pool duration or free-drug domain")
        }
        var state = initial, time = 0.0, step = min(policy.maximumStepSeconds, durationSeconds)
        var accepted = 0, rejected = 0, maximumError = 0.0, maximumBalance = 0.0
        let drugScale = max(initial.drugBalanceM, policy.absoluteToleranceM)
        let targetScale = max(initial.totalTargetM, model.baselineTarget.value, policy.absoluteToleranceM)
        while time < durationSeconds {
            try Task.checkCancellation()
            guard accepted + rejected < policy.maximumAttempts else { throw VivoKineticsError.capacity("finite-pool attempt budget") }
            step = min(step, durationSeconds - time)
            guard step > 0, time + step > time else { throw VivoKineticsError.numerical("finite-pool clock cannot advance") }
            let coarse = try split(state, model, step)
            let fine = try split(split(state, model, step * 0.5), model, step * 0.5)
            let error = zip(coarse.values, fine.values).map {
                abs($0 - $1) / (3 * (policy.absoluteToleranceM + policy.relativeTolerance * max(abs($0), abs($1))))
            }.max() ?? 0
            guard error.isFinite else { throw VivoKineticsError.numerical("finite-pool local error overflow") }
            if error > 1 {
                rejected += 1
                let reduced = step * max(0.1, 0.9 * pow(error, -1.0 / 3))
                guard reduced >= policy.minimumStepSeconds, reduced < step else {
                    throw VivoKineticsError.numerical("finite-pool tolerance requires a subminimum step")
                }
                step = reduced; continue
            }
            let balance = max(abs(fine.drugBalanceM - initial.drugBalanceM) / drugScale,
                              abs(fine.targetBalanceM - initial.targetBalanceM) / targetScale)
            guard balance.isFinite, balance <= policy.balanceRelativeTolerance,
                  fine.freeDrugM <= model.maximumFreeDrugM else {
                throw VivoKineticsError.numerical("finite-pool material balance or applicability domain failed")
            }
            maximumError = max(maximumError, error); maximumBalance = max(maximumBalance, balance)
            state = fine; time += step; accepted += 1
            step = min(policy.maximumStepSeconds, step * (error > 0 ? min(2, max(0.5, 0.9 * pow(error, -1.0 / 3))) : 2))
        }
        return .init(initial: initial, candidate: state, durationSeconds: durationSeconds,
                     acceptedSubsteps: accepted, rejectedAttempts: rejected,
                     maximumAcceptedLocalError: maximumError, maximumRelativeMaterialError: maximumBalance)
    }

    private static func split(_ input: VivoFiniteDrugState, _ model: VivoFiniteDrugModel, _ h: Double) throws -> VivoFiniteDrugState {
        var s = input
        turnover(&s, model, h * 0.5)
        drugDisposition(&s, model, h * 0.5)
        metaboliteClearance(&s, model, h * 0.5)
        complexConversion(&s, model, h * 0.5)
        try association(&s, model, h)
        complexConversion(&s, model, h * 0.5)
        metaboliteClearance(&s, model, h * 0.5)
        drugDisposition(&s, model, h * 0.5)
        turnover(&s, model, h * 0.5)
        try s.validate()
        return s
    }
    private static func turnover(_ s: inout VivoFiniteDrugState, _ m: VivoFiniteDrugModel, _ h: Double) {
        let x = m.targetTurnover.value * h, survival = exp(-x), loss = -expm1(-x)
        let initialTarget = s.totalTargetM, initialBound = s.reversibleComplexM + s.covalentComplexM
        // Stable x - (1-exp(-x)) for synthesis lost within the same substep.
        let sourceLoss = x < 1e-3 ? x * x * (0.5 + x * (-1.0/6 + x * (1.0/24 + x * (-1.0/120 + x/720)))) : x - loss
        s.freeTargetM = s.freeTargetM * survival + m.baselineTarget.value * loss
        s.reversibleComplexM *= survival; s.covalentComplexM *= survival
        s.inactiveMetaboliteM += initialBound * loss
        s.synthesizedTargetM += m.baselineTarget.value * x
        s.removedTargetM += initialTarget * loss + m.baselineTarget.value * sourceLoss
    }
    private static func drugDisposition(_ s: inout VivoFiniteDrugState, _ m: VivoFiniteDrugModel, _ h: Double) {
        let rate = m.freeDrugClearance.value + m.metabolism.value
        if rate == 0 { return }
        let lost = s.freeDrugM * (-expm1(-rate * h))
        s.freeDrugM *= exp(-rate * h)
        s.inactiveMetaboliteM += lost * (m.metabolism.value / rate)
        s.eliminatedDrugEquivalentM += lost * (m.freeDrugClearance.value / rate)
    }
    private static func metaboliteClearance(_ s: inout VivoFiniteDrugState, _ m: VivoFiniteDrugModel, _ h: Double) {
        let lost = s.inactiveMetaboliteM * (-expm1(-m.metaboliteClearance.value * h))
        s.inactiveMetaboliteM *= exp(-m.metaboliteClearance.value * h)
        s.eliminatedDrugEquivalentM += lost
    }
    private static func complexConversion(_ s: inout VivoFiniteDrugState, _ m: VivoFiniteDrugModel, _ h: Double) {
        let rate = m.dissociation.value + m.inactivation.value
        if rate == 0 { return }
        let lost = s.reversibleComplexM * (-expm1(-rate * h))
        let released = lost * (m.dissociation.value / rate)
        s.reversibleComplexM *= exp(-rate * h)
        s.freeDrugM += released; s.freeTargetM += released
        s.covalentComplexM += lost * (m.inactivation.value / rate)
    }
    private static func association(_ s: inout VivoFiniteDrugState, _ m: VivoFiniteDrugModel, _ h: Double) throws {
        let smaller = min(s.freeDrugM, s.freeTargetM), delta = abs(s.freeDrugM - s.freeTargetM)
        if smaller == 0 || m.association.value == 0 || h == 0 { return }
        let exposure = m.association.value * h
        guard exposure.isFinite, (exposure * max(smaller, delta)).isFinite else {
            throw VivoKineticsError.numerical("finite bimolecular rate-time product overflow")
        }
        let remaining: Double
        if delta == 0 || exposure * delta == 0 { remaining = smaller / (1 + exposure * smaller) }
        else {
            let x = exposure * delta
            remaining = smaller * (delta / (delta + smaller * (-expm1(-x)))) * exp(-x)
        }
        let extent = smaller - remaining
        guard remaining >= 0, extent >= 0 else { throw VivoKineticsError.numerical("bimolecular exact extent invalid") }
        if s.freeDrugM >= s.freeTargetM { s.freeDrugM = remaining + delta; s.freeTargetM = remaining }
        else { s.freeDrugM = remaining; s.freeTargetM = remaining + delta }
        s.reversibleComplexM += extent
    }
}
