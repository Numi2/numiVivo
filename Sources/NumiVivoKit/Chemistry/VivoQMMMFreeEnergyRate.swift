import Foundation

public enum VivoQMMMRateEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
    fileprivate var freeEnergyEnvironment: VivoQMMMFreeEnergyEnvironment {
        switch self { case .explicitSolution:return .explicitSolution;case .proteinEnvironment:return .proteinEnvironment }
    }
}

public struct VivoQMMMFreeEnergyRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate/v3"
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
    public static let schema="numivivo.org/qmmm-free-energy-rate-result/v3"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    /// Positive-velocity conventional TST flux before the explicit transmission probability.
    public let untransmittedFluxTSTRatePerSecond:Double
    public let positiveCoordinateVelocityNMPerPS:Double
    public let surfaceToReactantDensityPerNM:Double
    /// Eyring-equivalent barrier that reproduces the direct PMF flux exactly.
    /// This is not the simple PMF peak-minus-basin-minimum diagnostic.
    public let barrier:VivoActivationBarrier
    public let rateRequest:VivoTransitionStateRateRequest
    public let estimate:VivoTransitionStateRateEstimate
    public let parameter:VivoKineticParameter
    public let evidenceFingerprint:VivoFingerprint
    public let evidenceData:Data
    public let limitations:[String]
}

/// Protein/explicit-solution bridge. The configurational PMF contributes the
/// dividing-surface density / reactant-basin population. The coordinate mass
/// metric supplies the canonical positive velocity. Their product is a true
/// conventional classical TST flux with units of inverse time. Only then is an
/// Eyring-equivalent activation free energy constructed for the existing kinetic
/// contract. Biased trajectory elapsed time is never interpreted as a rate.
public enum VivoQMMMFreeEnergyRate {
    private static let gasConstantKJ=0.00831446261815324
    public static func calculate(_ request:VivoQMMMFreeEnergyRateRequest)throws->VivoQMMMFreeEnergyRateResult {
        try request.context.validate();try request.transmissionEvidence.validate(origin:request.transmissionOrigin)
        try request.freeEnergy.validate()
        let pmf=request.freeEnergy.analysis,provenance=request.freeEnergy.provenance,flux=request.freeEnergy.fluxNormalization
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
        let rt=gasConstantKJ*request.context.temperatureK
        let positiveVelocity=sqrt(rt*flux.inverseMassMetricPerDa/(2*Double.pi))
        let fluxRate=flux.surfaceToReactantDensityPerNM*positiveVelocity*1e12
        let thermal=VivoAtomicUnits.boltzmannJPerK*request.context.temperatureK/VivoTransitionStateTheory.planckJouleSecond
        guard positiveVelocity.isFinite,positiveVelocity>0,fluxRate.isFinite,fluxRate>0,
              thermal.isFinite,thermal>0,fluxRate<=thermal else {
            throw VivoKineticsError.unsupported("PMF flux is nonfinite or implies a barrierless/negative Eyring-equivalent activation free energy")
        }
        let equivalentBarrier = -rt*log(fluxRate/thermal)
        guard equivalentBarrier.isFinite,equivalentBarrier>=0 else { throw VivoKineticsError.numerical("PMF flux-to-barrier conversion") }
        let environment=request.environment == .proteinEnvironment ? "protein" : "explicit solution"
        let evidence=VivoKineticEvidence(source:"NumiVivo qualified QM/MM activation PMF",
            locator:"\(environment) PMF + mapped connectivity + exact Hamiltonian provenance + coordinate mass metric",
            sourceFingerprint:request.freeEnergy.evidenceFingerprint.hex)
        let barrier=VivoActivationBarrier(context:request.context,quantity:.activationGibbsFreeEnergy,
            referenceState:.preReactiveBoundComplex,value:equivalentBarrier,
            unit:.kilojoulesPerMol,conditionalStandardDeviation:pmf.conditionalStandardDeviationKJPerMol,
            method:"flux-normalized periodic QM/MM umbrella/MBAR conventional TST; Eyring-equivalent barrier",
            samplingDescription:request.samplingDescription,origin:.calculated,evidence:evidence)
        let rateRequest=VivoTransitionStateRateRequest(barrier:barrier,transmissionProbability:request.transmissionProbability,
            transmissionOrigin:request.transmissionOrigin,transmissionEvidence:request.transmissionEvidence)
        let derived=try VivoTransitionStateDerivation.calculate(rateRequest)
        let directLog=log(request.transmissionProbability)+log(fluxRate)
        guard abs(derived.estimate.naturalLogRatePerSecond-directLog)<1e-10 else {
            throw VivoKineticsError.numerical("PMF flux and Eyring-equivalent kinetic derivation disagree")
        }
        let limitations=[
            "The PMF profile height is diagnostic; rate normalization uses dividing-surface density divided by integrated reactant-basin population.",
            "The reported uncertainty is conditional on the sampled Hamiltonian, declared coordinate and approximate finite-sample PMF uncertainty.",
            "Alternative protonation states, reactive conformers and mechanisms require explicit population/pathway integration.",
            "Classical umbrella time is not interpreted as physical reaction time.",
            request.transmissionOrigin == .assumed ? "Dynamical recrossing and tunnelling remain unresolved because transmission is assumed." : "Transmission evidence is external to the PMF reconstruction.",
            "Force-field, QM level, adaptive-region calibration and finite-size errors are not included in the conditional sampling SD."
        ]
        return .init(schema:VivoQMMMFreeEnergyRateResult.schema,
            requestFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),
            untransmittedFluxTSTRatePerSecond:fluxRate,positiveCoordinateVelocityNMPerPS:positiveVelocity,
            surfaceToReactantDensityPerNM:flux.surfaceToReactantDensityPerNM,barrier:barrier,
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
