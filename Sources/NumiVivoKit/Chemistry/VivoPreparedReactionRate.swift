import Foundation

public enum VivoPreparedReactionEnvironment: String, Codable, Sendable {
    case isolatedReactantComplex, continuumReactantComplex, proteinEnvironment
}

public struct VivoPreparedReactionRateRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/prepared-reaction-rate/v1"
    public var schema:String
    public var context:VivoKineticContext
    public var environment:VivoPreparedReactionEnvironment
    public var preparation:VivoMolecularPreparationRequest
    public var sampling:VivoMolecularSamplingResult
    public var transitionState:VivoTransitionStateTheoryRequest
    /// An explicit declaration of how the classical sampling model relates to
    /// the electronic model. It does not stand in for QM reweighting or a PMF.
    public var modelTransferDescription:String
    /// Computed/experimental transmission requires actual retained evidence
    /// bytes, not a digest computed from the coefficient declaration itself.
    public var transmissionEvidence:VivoKineticEvidence?
    public var transmissionEvidenceData:Data?
    public init(context:VivoKineticContext,environment:VivoPreparedReactionEnvironment,
                preparation:VivoMolecularPreparationRequest,sampling:VivoMolecularSamplingResult,
                transitionState:VivoTransitionStateTheoryRequest,modelTransferDescription:String,
                transmissionEvidence:VivoKineticEvidence? = nil,transmissionEvidenceData:Data? = nil) {
        schema=Self.schema;self.context=context;self.environment=environment;self.preparation=preparation
        self.sampling=sampling;self.transitionState=transitionState;self.modelTransferDescription=modelTransferDescription
        self.transmissionEvidence=transmissionEvidence;self.transmissionEvidenceData=transmissionEvidenceData
    }
}

public struct VivoPreparedReactionRateResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/prepared-reaction-rate-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let preparedStructureFingerprint:VivoFingerprint
    public let classicalSystemFingerprint:VivoFingerprint
    public let molecularRate:VivoTransitionStateTheoryResult
    public let kineticRequest:VivoTransitionStateRateRequest
    public let kineticEstimate:VivoTransitionStateRateEstimate
    public let kineticParameter:VivoKineticParameter
    public let kineticEvidenceFingerprint:VivoFingerprint
    public let kineticEvidenceData:Data
    public let unresolvedContributions:[String]
    public let interpretation:String
}

/// Integration of the existing molecular TST and kinetic derivation paths.
/// A locally harmonic model rate remains an ASSUMED kinetic parameter even when
/// all numerical checks pass. Sampling a classical observable does not turn an
/// electronic polarization calculation into a complete activation free energy.
public enum VivoPreparedReactionRate {
    public static let interpretation="Context-bound conditional local harmonic TST model, reconstructed from explicit molecular preparation, replica diagnostics and mapped reaction connectivity. The kinetic parameter remains assumed with unknown total uncertainty; not a protein inactivation prediction or a converged solution reaction free energy."

    public static func calculate(_ request:VivoPreparedReactionRateRequest) throws -> VivoPreparedReactionRateResult {
        try request.context.validate()
        guard request.schema==VivoPreparedReactionRateRequest.schema,
              !request.modelTransferDescription.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,
              request.modelTransferDescription.utf8.count<=16384 else {
            throw VivoKineticsError.invalid("prepared rate schema or classical-to-electronic model declaration")
        }
        guard request.environment != .proteinEnvironment else {
            throw VivoKineticsError.unsupported("protein rates require environment-consistent reaction sampling, QM/MM free-energy corrections and constrained thermochemistry; free-molecule RRHO or electrostatic solvent convergence cannot supply them")
        }
        let transmission=request.transitionState.transmission
        try transmission.validate()
        let origin:VivoKineticOrigin
        switch transmission.kind {case .assumed:origin = .assumed;case .computedDynamics:origin = .calculated;case .experimental:origin = .measured}
        let transmissionEvidence:VivoKineticEvidence
        if let evidence=request.transmissionEvidence {
            try evidence.validate(origin:origin)
            guard let bytes=request.transmissionEvidenceData,!bytes.isEmpty,bytes.count<=1024*1024,
                  evidence.sourceFingerprint == (try VivoCanonicalJSON.fingerprint(bytes)).hex else {
                throw VivoKineticsError.invalid("transmission evidence bytes do not match their immutable source fingerprint")
            }
            transmissionEvidence=evidence
        } else {
            guard origin == .assumed,request.transmissionEvidenceData==nil else {
                throw VivoKineticsError.invalid("computed or experimental transmission requires retained evidence bytes and a source locator")
            }
            transmissionEvidence = .init(source:transmission.sourceIdentifier,
                locator:"declared transmission assumption; no independent dynamics evidence",
                sourceFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(transmission)).hex)
        }
        let prepared=try VivoMolecularPreparation.prepare(request.preparation)
        guard let compiled=prepared.compiledForceField else {
            throw VivoKineticsError.invalid("prepared rate requires an explicitly parameterized native sampling system")
        }
        let structureID=try VivoStructureCodec.fingerprint(prepared.structure),systemID=try compiled.system.fingerprint()
        try VivoMolecularSampling.validate(request.sampling)
        let sampling=request.sampling.request,context=request.context
        guard request.sampling.converged,sampling.structureFingerprint==structureID,sampling.systemFingerprint==systemID,
              sampling.contextIdentifier==context.hostContext,request.preparation.pH==context.pH,
              request.preparation.microstateIdentifier==context.chemicalState,
              sampling.replicas.allSatisfy({$0.configuration.targetTemperatureK==context.temperatureK}) else {
            throw VivoKineticsError.invalid("rate context, prepared microstate or sampling Hamiltonian/temperature mismatch")
        }
        // Temperature/energy-only convergence cannot be used as molecular
        // conformational evidence for this bridge.
        guard sampling.observables.contains(where:{ observable in
            switch observable.kind {
            case .distance,.distanceDifference,.torsionCosine,.torsionSine:return true
            default:return false
            }
        }) else { throw VivoKineticsError.invalid("prepared rate requires declared molecular geometry observables") }
        let saddle=request.transitionState.connectivity.request.saddle
        let model=saddle.request.model,tc=saddle.thermochemistry.configuration
        guard tc.temperatureK==context.temperatureK,model.system.pointCharges.isEmpty,
              (request.environment == .continuumReactantComplex) == (model.solvent != nil),
              request.transitionState.transmission.coefficient<=1 else {
            throw VivoKineticsError.unsupported("rate environment, temperature or classical transmission convention differs")
        }
        guard let endpoint=request.transitionState.connectivity.request.endpoints.first(where: {
            $0.identifier==request.transitionState.reactantEndpointIdentifier
        }),endpoint.components.count==1 else {
            throw VivoKineticsError.unsupported("bound-complex kinetic integration is unimolecular; separated reactants retain the existing molecular TST standard-state rate")
        }
        let indices=model.system.nuclei.compactMap(\.structureAtomIndex)
        guard indices.count==model.system.nuclei.count,Set(indices).count==indices.count,
              Set(indices)==Set(prepared.structure.atoms.map(\.index)) else {
            throw VivoKineticsError.invalid("isolated/continuum rate must retain every prepared physical atom with an explicit map; no protein-to-cluster relabeling")
        }
        for (nucleus,index) in zip(model.system.nuclei,indices) {
            guard nucleus.atomicNumber==Int(prepared.structure.atoms[Int(index)].element.atomicNumber) else {
                throw VivoKineticsError.invalid("prepared-to-electronic element mapping differs")
            }
        }
        let electrons=model.system.alphaElectrons+model.system.betaElectrons
        let nuclearCharge=model.system.nuclei.reduce(0){$0+$1.atomicNumber}
        guard nuclearCharge-electrons==request.preparation.expectedFormalCharge else {
            throw VivoKineticsError.invalid("electronic charge sector differs from the prepared protonation microstate")
        }
        let molecular=try VivoTransitionStateTheory.estimate(request.transitionState)
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        var missing=[
            "Global conformer and protonation-state populations are not supplied by local RRHO.",
            "Classical trajectory diagnostics do not establish QM configurational equilibrium or a reaction PMF.",
            "Sampling-to-electronic model transfer has not been free-energy reweighted.",
            "Anharmonic/soft-mode thermochemistry and total model uncertainty remain unresolved."]
        if model.solvent != nil {
            missing += ["Cavitation, dispersion/repulsion and solvent configurational/reorganization contributions are incomplete.",
                        "A concentration standard-state conversion does not supply missing solution entropy."]
        }
        if request.transitionState.transmission.kind == .assumed {
            missing.append("Dynamical recrossing and tunneling are not determined by an assumed transmission coefficient.")
        }
        let evidence=VivoKineticEvidence(source:"NumiVivo prepared molecular reaction workflow",
            locator:"prepared request, native parameters, replica measurements, mapped connectivity and harmonic TST; conditional model only",
            sourceFingerprint:requestID.hex)
        let barrier=VivoActivationBarrier(context:context,quantity:.activationGibbsFreeEnergy,
            referenceState:.preReactiveBoundComplex,value:molecular.barrier.activationGibbsHartree*VivoAtomicUnits.hartreeInKJPerMol,
            unit:.kilojoulesPerMol,method:"native local harmonic activation Gibbs MODEL; missing contributions remain explicit",
            samplingDescription:request.modelTransferDescription,origin:.assumed,evidence:evidence)
        let kinetic=VivoTransitionStateRateRequest(barrier:barrier,transmissionProbability:transmission.coefficient,
            transmissionOrigin:origin,transmissionEvidence:transmissionEvidence)
        let derived=try VivoTransitionStateDerivation.calculate(kinetic)
        guard derived.parameter.origin == .assumed,
              abs(derived.estimate.naturalLogRatePerSecond-molecular.logRateConstant)<1e-8 else {
            throw VivoKineticsError.numerical("molecular TST and canonical kinetic derivation disagree")
        }
        return .init(schema:VivoPreparedReactionRateResult.schema,requestFingerprint:requestID,
            preparedStructureFingerprint:structureID,classicalSystemFingerprint:systemID,molecularRate:molecular,
            kineticRequest:kinetic,kineticEstimate:derived.estimate,kineticParameter:derived.parameter,
            kineticEvidenceFingerprint:derived.evidenceFingerprint,kineticEvidenceData:derived.evidenceData,
            unresolvedContributions:missing,interpretation:interpretation)
    }
    public static func validate(_ result:VivoPreparedReactionRateResult,request:VivoPreparedReactionRateRequest) throws {
        guard result == (try calculate(request)) else { throw VivoKineticsError.invalid("prepared rate does not reconstruct from its complete source request") }
    }
    /// Reuses the existing exact-context inactivation update; association,
    /// dissociation and turnover are not derived from a bound-complex barrier.
    public static func applying(_ result:VivoPreparedReactionRateResult,request:VivoPreparedReactionRateRequest,
                                to model:VivoCovalentKineticPack) throws -> VivoCovalentKineticPack {
        try validate(result,request:request)
        return try VivoTransitionStateDerivation.calculate(result.kineticRequest).applyingInactivation(to:model)
    }
}
