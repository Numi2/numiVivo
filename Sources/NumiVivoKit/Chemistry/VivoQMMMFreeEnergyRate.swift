import Foundation

public enum VivoQMMMRateEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
    fileprivate var freeEnergyEnvironment: VivoQMMMFreeEnergyEnvironment {
        switch self { case .explicitSolution:return .explicitSolution;case .proteinEnvironment:return .proteinEnvironment }
    }
}

public struct VivoQMMMSampledFluxNormalization: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-sampled-flux-normalization/v1"
    public let schema:String
    public let freeEnergyEvidenceFingerprint:VivoFingerprint
    public let surfaceEvidenceFingerprint:VivoFingerprint
    public let positiveCoordinateVelocityNMPerPS:Double
    public let positiveCoordinateVelocityStandardErrorNMPerPS:Double
    public let effectiveSurfaceSamples:Double
    public init(freeEnergyEvidenceFingerprint:VivoFingerprint,surfaceEvidenceFingerprint:VivoFingerprint,
                positiveCoordinateVelocityNMPerPS:Double,positiveCoordinateVelocityStandardErrorNMPerPS:Double,
                effectiveSurfaceSamples:Double) {
        schema=Self.schema;self.freeEnergyEvidenceFingerprint=freeEnergyEvidenceFingerprint
        self.surfaceEvidenceFingerprint=surfaceEvidenceFingerprint
        self.positiveCoordinateVelocityNMPerPS=positiveCoordinateVelocityNMPerPS
        self.positiveCoordinateVelocityStandardErrorNMPerPS=positiveCoordinateVelocityStandardErrorNMPerPS
        self.effectiveSurfaceSamples=effectiveSurfaceSamples
    }
    public func validate(freeEnergy:VivoQMMMQualifiedActivationFreeEnergy)throws {
        guard schema==Self.schema,freeEnergyEvidenceFingerprint==freeEnergy.evidenceFingerprint,
              positiveCoordinateVelocityNMPerPS.isFinite,positiveCoordinateVelocityNMPerPS>0,
              positiveCoordinateVelocityStandardErrorNMPerPS.isFinite,positiveCoordinateVelocityStandardErrorNMPerPS>=0,
              effectiveSurfaceSamples.isFinite,effectiveSurfaceSamples>=2 else {
            throw VivoKineticsError.invalid("sampled QM/MM surface flux normalization")
        }
    }
}

public enum VivoQMMMSurfaceFluxNormalization {
    public static func make(surface:VivoQMMMSurfaceEnsembleResult,
                            request:VivoQMMMSurfaceEnsembleRequest)throws->VivoQMMMSampledFluxNormalization {
        try VivoQMMMSurfaceEnsemble.validate(surface,request:request)
        guard surface.converged,surface.freeEnergyEvidenceFingerprint==request.freeEnergy.evidenceFingerprint else {
            throw VivoKineticsError.invalid("nonconverged or mismatched dividing-surface ensemble")
        }
        return .init(freeEnergyEvidenceFingerprint:request.freeEnergy.evidenceFingerprint,
                     surfaceEvidenceFingerprint:surface.evidenceFingerprint,
                     positiveCoordinateVelocityNMPerPS:surface.positiveCoordinateVelocityNMPerPS,
                     positiveCoordinateVelocityStandardErrorNMPerPS:surface.positiveCoordinateVelocityStandardErrorNMPerPS,
                     effectiveSurfaceSamples:surface.effectiveSurfaceSamples)
    }
}

public struct VivoQMMMFreeEnergyRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate/v4"
    public var schema:String
    public var context:VivoKineticContext
    public var environment:VivoQMMMRateEnvironment
    public var freeEnergy:VivoQMMMQualifiedActivationFreeEnergy
    /// Required when the coordinate mass metric is geometry dependent or when
    /// constrained dynamics should define the actual canonical velocity law.
    public var sampledFluxNormalization:VivoQMMMSampledFluxNormalization?
    public var transmissionProbability:Double
    public var transmissionOrigin:VivoKineticOrigin
    public var transmissionEvidence:VivoKineticEvidence
    public var samplingDescription:String
    public init(context:VivoKineticContext,environment:VivoQMMMRateEnvironment,freeEnergy:VivoQMMMQualifiedActivationFreeEnergy,
                sampledFluxNormalization:VivoQMMMSampledFluxNormalization?=nil,
                transmissionProbability:Double=1,transmissionOrigin:VivoKineticOrigin = .assumed,
                transmissionEvidence:VivoKineticEvidence,samplingDescription:String) {
        schema=Self.schema;self.context=context;self.environment=environment;self.freeEnergy=freeEnergy
        self.sampledFluxNormalization=sampledFluxNormalization;self.transmissionProbability=transmissionProbability
        self.transmissionOrigin=transmissionOrigin;self.transmissionEvidence=transmissionEvidence
        self.samplingDescription=samplingDescription
    }
}

public struct VivoQMMMFreeEnergyRateResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-rate-result/v4"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let untransmittedFluxTSTRatePerSecond:Double
    public let positiveCoordinateVelocityNMPerPS:Double
    public let positiveCoordinateVelocityStandardErrorNMPerPS:Double?
    public let velocityNormalization:String
    public let surfaceToReactantDensityPerNM:Double
    public let barrier:VivoActivationBarrier
    public let rateRequest:VivoTransitionStateRateRequest
    public let estimate:VivoTransitionStateRateEstimate
    public let parameter:VivoKineticParameter
    public let evidenceFingerprint:VivoFingerprint
    public let evidenceData:Data
    public let limitations:[String]
}

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
        let positiveVelocity:Double,velocitySE:Double?,velocityOrigin:String,fluxEvidence:String
        if let sampled=request.sampledFluxNormalization {
            try sampled.validate(freeEnergy:request.freeEnergy)
            positiveVelocity=sampled.positiveCoordinateVelocityNMPerPS
            velocitySE=sampled.positiveCoordinateVelocityStandardErrorNMPerPS
            velocityOrigin="sampled canonical dividing-surface velocity"
            fluxEvidence=sampled.surfaceEvidenceFingerprint.hex
        } else if let metric=flux.inverseMassMetricPerDa {
            positiveVelocity=sqrt(rt*metric/(2*Double.pi));velocitySE=nil
            velocityOrigin="analytic unconstrained Cartesian mass metric"
            fluxEvidence=request.freeEnergy.evidenceFingerprint.hex
        } else {
            throw VivoKineticsError.unsupported("geometry-dependent/shared-atom reaction coordinate requires sampled dividing-surface velocity normalization")
        }
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
            locator:"\(environment) PMF + mapped connectivity + exact Hamiltonian/sampling provenance + \(velocityOrigin); fluxEvidence=\(fluxEvidence)",
            sourceFingerprint:request.freeEnergy.evidenceFingerprint.hex)
        let barrier=VivoActivationBarrier(context:request.context,quantity:.activationGibbsFreeEnergy,
            referenceState:.preReactiveBoundComplex,value:equivalentBarrier,
            unit:.kilojoulesPerMol,conditionalStandardDeviation:nil,
            method:"flux-normalized periodic QM/MM umbrella/MBAR conventional TST; Eyring-equivalent barrier",
            samplingDescription:request.samplingDescription,origin:.calculated,evidence:evidence)
        let rateRequest=VivoTransitionStateRateRequest(barrier:barrier,transmissionProbability:request.transmissionProbability,
            transmissionOrigin:request.transmissionOrigin,transmissionEvidence:request.transmissionEvidence)
        let derived=try VivoTransitionStateDerivation.calculate(rateRequest)
        let directLog=log(request.transmissionProbability)+log(fluxRate)
        guard abs(derived.estimate.naturalLogRatePerSecond-directLog)<1e-10 else {
            throw VivoKineticsError.numerical("PMF flux and Eyring-equivalent kinetic derivation disagree")
        }
        var limitations=[
            "The PMF profile height is diagnostic; rate normalization uses dividing-surface density divided by integrated reactant-basin population.",
            "Single-PMF profile-height uncertainty is not propagated into the flux rate because it is not the same statistical observable; use independent replicated PMFs for sampling dispersion.",
            "Alternative protonation states, reactive conformers and mechanisms require explicit population/pathway integration.",
            "Classical umbrella time is not interpreted as physical reaction time.",
            request.transmissionOrigin == .assumed ? "Dynamical recrossing and tunnelling remain unresolved because transmission is assumed." : "Transmission is separately qualified against the same PMF Hamiltonian.",
            "Force-field, QM level, adaptive-region calibration and finite-size errors are outside the replicated sampling dispersion."
        ]
        if velocitySE != nil { limitations.append("Surface velocity uncertainty is reported separately and is not yet combined with PMF replica dispersion.") }
        return .init(schema:VivoQMMMFreeEnergyRateResult.schema,
            requestFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),
            untransmittedFluxTSTRatePerSecond:fluxRate,positiveCoordinateVelocityNMPerPS:positiveVelocity,
            positiveCoordinateVelocityStandardErrorNMPerPS:velocitySE,velocityNormalization:velocityOrigin,
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
