import Foundation

public enum VivoQMMMRateEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
    fileprivate var freeEnergyEnvironment: VivoQMMMFreeEnergyEnvironment {
        switch self { case .explicitSolution:return .explicitSolution;case .proteinEnvironment:return .proteinEnvironment }
    }
}

public struct VivoQMMMFreeEnergyRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate/v2"
    public var schema:String
    public var context:VivoKineticContext
    public var environment:VivoQMMMRateEnvironment
    public var freeEnergy:VivoQMMMQualifiedActivationFreeEnergy
    public var transmissionProbability:Double
    public var transmissionOrigin:VivoKineticOrigin
    public var transmissionEvidence:VivoKineticEvidence
    public var samplingDescription:String
    public init(context:VivoKineticContext,environment:VivoQMMMRateEnvironment,freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,
                transmissionProbability:Double=1,transmissionOrigin:VivoKineticOrigin = .assumed,
                transmissionEvidence:VivoKineticEvidence,samplingDescription:String) {
        schema=Self.schema;self.context=context;self.environment=environment;self.freeEnergy=freeEnergy
        self.transmissionProbability=transmissionProbability;self.transmissionOrigin=transmissionOrigin
        self.transmissionEvidence=transmissionEvidence;self.samplingDescription=samplingDescription
    }
}

public struct VivoQMMMFreeEnergyRateResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate-result/v2"
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

/// Protein/explicit-solution bridge. The barrier must be a converged PMF bound
/// to the exact sampled structure, classical system, BO provider and MD profile;
/// a bare PMF or local stationary-point result is not sufficient.
public enum VivoQMMMFreeEnergyRate {
    public static func calculate(_ request:VivoQMMMFreeEnergyRateRequest)throws->VivoQMMMFreeEnergyRateResult {
        try request.context.validate();try request.transmissionEvidence.validate(origin:request.transmissionOrigin)
        try request.freeEnergy.validate()
        let pmf=request.freeEnergy.analysis,provenance=request.freeEnergy.provenance
        guard request.schema==VivoQMMMFreeEnergyRateRequest.schema,pmf.converged,
              request.context.temperatureK==pmf.temperatureK,
              request.context.chemicalState==provenance.chemicalState,
              request.context.hostContext==provenance.environmentIdentifier,
              request.environment.freeEnergyEnvironment==provenance.environment,
              request.transmissionProbability.isFinite,request.transmissionProbability>0,request.transmissionProbability<=1,
              !request.samplingDescription.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,
              request.samplingDescription.utf8.count<=16384 else {
            throw VivoKineticsError.invalid("QM/MM rate context/environment differs from the qualified PMF or transmission is invalid")
        }
        let environment=request.environment == .proteinEnvironment ? "protein" : "explicit solution"
        let evidence=VivoKineticEvidence(source:"NumiVivo qualified QM/MM activation free energy",
            locator:"\(environment) PMF; exact system/provider/dynamics provenance + retained MBAR traces",
            sourceFingerprint:request.freeEnergy.evidenceFingerprint.hex)
        let barrier=VivoActivationBarrier(context:request.context,quantity:.activationGibbsFreeEnergy,
            referenceState:.preReactiveBoundComplex,value:pmf.activationFreeEnergyKJPerMol,
            unit:.kilojoulesPerMol,conditionalStandardDeviation:pmf.conditionalStandardDeviationKJPerMol,
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
