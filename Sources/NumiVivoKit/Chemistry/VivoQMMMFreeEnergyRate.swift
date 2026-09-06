import Foundation

public enum VivoQMMMRateEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
}

public struct VivoQMMMFreeEnergyRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate/v1"
    public var schema:String
    public var context:VivoKineticContext
    public var environment:VivoQMMMRateEnvironment
    public var freeEnergy:VivoQMMMActivationFreeEnergyResult
    public var transmissionProbability:Double
    public var transmissionOrigin:VivoKineticOrigin
    public var transmissionEvidence:VivoKineticEvidence
    public var samplingDescription:String
    public init(context:VivoKineticContext,environment:VivoQMMMRateEnvironment,freeEnergy:VivoQMMMActivationFreeEnergyResult,
                transmissionProbability:Double=1,transmissionOrigin:VivoKineticOrigin = .assumed,
                transmissionEvidence:VivoKineticEvidence,samplingDescription:String) {
        schema=Self.schema;self.context=context;self.environment=environment;self.freeEnergy=freeEnergy
        self.transmissionProbability=transmissionProbability;self.transmissionOrigin=transmissionOrigin
        self.transmissionEvidence=transmissionEvidence;self.samplingDescription=samplingDescription
    }
}

public struct VivoQMMMFreeEnergyRateResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let barrier:VivoActivationBarrier
    public let rateRequest:VivoTransitionStateRateRequest
    public let estimate:VivoTransitionStateRateEstimate
    public let parameter:VivoKineticParameter
    public let evidenceFingerprint:VivoFingerprint
    public let evidenceData:Data
    public let limitations:[String]
}

/// This is the protein/explicit-solution bridge that the local-RRHO prepared-rate
/// adapter intentionally could not provide. The activation Gibbs free energy is
/// taken from a converged environment-consistent PMF, not from electrostatic
/// solvent convergence or an isolated stationary-point difference.
public enum VivoQMMMFreeEnergyRate {
    public static func calculate(_ request:VivoQMMMFreeEnergyRateRequest)throws->VivoQMMMFreeEnergyRateResult {
        try request.context.validate();try request.transmissionEvidence.validate(origin:request.transmissionOrigin)
        try VivoQMMMFreeEnergy.validate(request.freeEnergy)
        guard request.schema==VivoQMMMFreeEnergyRateRequest.schema,request.freeEnergy.converged,
              request.context.temperatureK==request.freeEnergy.temperatureK,
              request.transmissionProbability.isFinite,request.transmissionProbability>0,request.transmissionProbability<=1,
              !request.samplingDescription.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,
              request.samplingDescription.utf8.count<=16384 else {
            throw VivoKineticsError.invalid("QM/MM free-energy rate context, convergence or transmission")
        }
        let data=try VivoCanonicalJSON.encode(request.freeEnergy),source=try VivoCanonicalJSON.fingerprint(data)
        let environment=request.environment == .proteinEnvironment ? "protein" : "explicit solution"
        let evidence=VivoKineticEvidence(source:"NumiVivo QM/MM umbrella/MBAR activation free energy",
            locator:"\(environment) PMF; retained window traces and reconstructible MBAR analysis",
            sourceFingerprint:source.hex)
        let barrier=VivoActivationBarrier(context:request.context,quantity:.activationGibbsFreeEnergy,
            referenceState:.preReactiveBoundComplex,value:request.freeEnergy.activationFreeEnergyKJPerMol,
            unit:.kilojoulesPerMol,conditionalStandardDeviation:request.freeEnergy.conditionalStandardDeviationKJPerMol,
            method:"periodic QM/MM conservative umbrella sampling + unbinned MBAR",
            samplingDescription:request.samplingDescription,origin:.calculated,evidence:evidence)
        let rateRequest=VivoTransitionStateRateRequest(barrier:barrier,transmissionProbability:request.transmissionProbability,
            transmissionOrigin:request.transmissionOrigin,transmissionEvidence:request.transmissionEvidence)
        let derived=try VivoTransitionStateDerivation.calculate(rateRequest)
        let limitations=[
            "The reported uncertainty is conditional on the sampled Hamiltonian and declared reaction coordinate.",
            "Alternative protonation states, reactive conformers and mechanisms require explicit population/pathway integration.",
            "Classical umbrella time is not interpreted as physical reaction time.",
            request.transmissionOrigin == .assumed ? "Dynamical recrossing and tunnelling remain unresolved because transmission is assumed." : "Transmission evidence is external to the PMF reconstruction.",
            "Force-field, QM level, adaptive-region calibration and finite-size errors are not included in the conditional sampling SD."
        ]
        return .init(schema:VivoQMMMFreeEnergyRateResult.schema,
            requestFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),barrier:barrier,
            rateRequest:rateRequest,estimate:derived.estimate,parameter:derived.parameter,
            evidenceFingerprint:derived.evidenceFingerprint,evidenceData:derived.evidenceData,limitations:limitations)
    }
    public static func validate(_ result:VivoQMMMFreeEnergyRateResult,request:VivoQMMMFreeEnergyRateRequest)throws {
        guard result == (try calculate(request)) else { throw VivoKineticsError.invalid("QM/MM free-energy rate does not reconstruct") }
    }
    public static func applying(_ result:VivoQMMMFreeEnergyRateResult,request:VivoQMMMFreeEnergyRateRequest,
                                to model:VivoCovalentKineticPack)throws->VivoCovalentKineticPack {
        try validate(result,request:request)
        return try VivoTransitionStateDerivation.calculate(result.rateRequest).applyingInactivation(to:model)
    }
}
