import Foundation

public enum VivoQMMMFreeEnergyEnvironment: String, Codable, Sendable {
    case explicitSolution
    case proteinEnvironment
}

public struct VivoQMMMReactionConnectivityBinding: Codable, Sendable, Equatable {
    public var connectivityFingerprint: VivoFingerprint
    public var reactantEndpointIdentifier: String
    public var productEndpointIdentifier: String
    public var mappedReactionAtomIndices: [UInt32]
    public init(connectivityFingerprint: VivoFingerprint,reactantEndpointIdentifier: String,
                productEndpointIdentifier: String,mappedReactionAtomIndices: [UInt32]) {
        self.connectivityFingerprint=connectivityFingerprint;self.reactantEndpointIdentifier=reactantEndpointIdentifier
        self.productEndpointIdentifier=productEndpointIdentifier;self.mappedReactionAtomIndices=mappedReactionAtomIndices
    }
    public func validate(coordinate: VivoQMMMReactionCoordinate) throws {
        guard !reactantEndpointIdentifier.isEmpty,!productEndpointIdentifier.isEmpty,
              reactantEndpointIdentifier != productEndpointIdentifier,
              mappedReactionAtomIndices==coordinate.atomIndices else {
            throw VivoChemistryError.invalid("QM/MM PMF reaction-connectivity binding")
        }
    }
}

public struct VivoQMMMFreeEnergyProvenance: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-provenance/v3"
    public var schema:String
    public var structureFingerprint:VivoFingerprint
    public var systemFingerprint:VivoFingerprint
    public var baseProviderFingerprint:VivoFingerprint
    public var dynamicsFingerprint:VivoFingerprint
    public var chemicalState:String
    public var environment:VivoQMMMFreeEnergyEnvironment
    public var environmentIdentifier:String
    public var methodDescription:String
    public var reactionConnectivity:VivoQMMMReactionConnectivityBinding
    public init(structureFingerprint:VivoFingerprint,systemFingerprint:VivoFingerprint,baseProviderFingerprint:VivoFingerprint,
                dynamicsFingerprint:VivoFingerprint,chemicalState:String,environment:VivoQMMMFreeEnergyEnvironment,
                environmentIdentifier:String,methodDescription:String,reactionConnectivity:VivoQMMMReactionConnectivityBinding) {
        schema=Self.schema;self.structureFingerprint=structureFingerprint;self.systemFingerprint=systemFingerprint
        self.baseProviderFingerprint=baseProviderFingerprint;self.dynamicsFingerprint=dynamicsFingerprint
        self.chemicalState=chemicalState;self.environment=environment;self.environmentIdentifier=environmentIdentifier
        self.methodDescription=methodDescription;self.reactionConnectivity=reactionConnectivity
    }
    public func validate(coordinate:VivoQMMMReactionCoordinate) throws {
        try reactionConnectivity.validate(coordinate:coordinate)
        guard schema==Self.schema,
              !chemicalState.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,chemicalState.utf8.count<=1024,
              !environmentIdentifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,environmentIdentifier.utf8.count<=4096,
              !methodDescription.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,methodDescription.utf8.count<=16384 else {
            throw VivoChemistryError.invalid("QM/MM free-energy provenance identity")
        }
    }
}

public struct VivoQMMMFluxNormalization: Codable, Sendable, Equatable {
    /// exp[-beta A(xi*)] / integral_R exp[-beta A(xi)] dxi, in nm^-1.
    /// The arbitrary PMF additive constant cancels exactly.
    public let surfaceToReactantDensityPerNM: Double
    /// g_xi = sum_i |d xi/d r_i|^2 / m_i. For the supported distance and
    /// distance-difference coordinates every nonzero Cartesian gradient has
    /// unit norm, so this is a geometry-independent sum of inverse masses.
    public let inverseMassMetricPerDa: Double
    public init(surfaceToReactantDensityPerNM:Double,inverseMassMetricPerDa:Double) {
        self.surfaceToReactantDensityPerNM=surfaceToReactantDensityPerNM;self.inverseMassMetricPerDa=inverseMassMetricPerDa
    }
    public func validate() throws {
        guard surfaceToReactantDensityPerNM.isFinite,surfaceToReactantDensityPerNM>0,
              inverseMassMetricPerDa.isFinite,inverseMassMetricPerDa>0 else {
            throw VivoChemistryError.invalid("QM/MM PMF flux normalization or reaction-coordinate mass metric")
        }
    }
}

public struct VivoQMMMQualifiedActivationFreeEnergy: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-qualified-activation-free-energy/v3"
    public let schema:String
    public let analysis:VivoQMMMActivationFreeEnergyResult
    public let provenance:VivoQMMMFreeEnergyProvenance
    public let fluxNormalization:VivoQMMMFluxNormalization
    public let evidenceFingerprint:VivoFingerprint

    init(analysis:VivoQMMMActivationFreeEnergyResult,provenance:VivoQMMMFreeEnergyProvenance,
         fluxNormalization:VivoQMMMFluxNormalization) throws {
        try provenance.validate(coordinate:analysis.coordinate);try fluxNormalization.validate();try VivoQMMMFreeEnergy.validate(analysis)
        guard analysis.converged else { throw VivoChemistryError.convergence("only a converged PMF can be qualified") }
        self.schema=Self.schema;self.analysis=analysis;self.provenance=provenance;self.fluxNormalization=fluxNormalization
        struct Evidence:Codable {
            let analysis:VivoQMMMActivationFreeEnergyResult
            let provenance:VivoQMMMFreeEnergyProvenance
            let fluxNormalization:VivoQMMMFluxNormalization
        }
        self.evidenceFingerprint=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
            Evidence(analysis:analysis,provenance:provenance,fluxNormalization:fluxNormalization)))
    }
    public func validate() throws {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("qualified activation-free-energy schema") }
        try provenance.validate(coordinate:analysis.coordinate);try fluxNormalization.validate();try VivoQMMMFreeEnergy.validate(analysis)
        guard analysis.converged else { throw VivoChemistryError.convergence("qualified PMF is not converged") }
        let rebuilt=try VivoQMMMFreeEnergyQualification.surfaceToReactantDensityPerNM(analysis)
        guard abs(rebuilt-fluxNormalization.surfaceToReactantDensityPerNM)<=max(1e-12,abs(rebuilt)*1e-10) else {
            throw VivoChemistryError.invalid("qualified PMF surface/reactant normalization differs from its profile")
        }
        struct Evidence:Codable {
            let analysis:VivoQMMMActivationFreeEnergyResult
            let provenance:VivoQMMMFreeEnergyProvenance
            let fluxNormalization:VivoQMMMFluxNormalization
        }
        let expected=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
            Evidence(analysis:analysis,provenance:provenance,fluxNormalization:fluxNormalization)))
        guard evidenceFingerprint==expected else { throw VivoChemistryError.invalid("qualified PMF evidence fingerprint mismatch") }
    }
}

public enum VivoQMMMFreeEnergyQualification {
    private static let gasConstantKJ=0.00831446261815324

    static func surfaceToReactantDensityPerNM(_ result:VivoQMMMActivationFreeEnergyResult) throws -> Double {
        let profile=result.profile,cfg=result.configuration
        guard profile.count>=2,profile.map(\.coordinateNM)==profile.map(\.coordinateNM).sorted(),
              cfg.reactantRangeNM.lowerBound>=profile[0].coordinateNM,
              cfg.reactantRangeNM.upperBound<=profile.last!.coordinateNM,
              cfg.dividingSurfaceNM>=profile[0].coordinateNM,cfg.dividingSurfaceNM<=profile.last!.coordinateNM else {
            throw VivoChemistryError.invalid("qualified PMF profile does not fully cover reactant basin and dividing surface")
        }
        let beta=1/(gasConstantKJ*result.temperatureK)
        func energy(_ x:Double)throws->Double {
            if x==profile[0].coordinateNM { return profile[0].relativeFreeEnergyKJPerMol }
            for i in 1..<profile.count where x<=profile[i].coordinateNM {
                let a=profile[i-1],b=profile[i],width=b.coordinateNM-a.coordinateNM
                guard width>0 else { throw VivoChemistryError.invalid("qualified PMF grid is not strictly increasing") }
                let t=(x-a.coordinateNM)/width
                return a.relativeFreeEnergyKJPerMol+t*(b.relativeFreeEnergyKJPerMol-a.relativeFreeEnergyKJPerMol)
            }
            return profile.last!.relativeFreeEnergyKJPerMol
        }
        let lo=cfg.reactantRangeNM.lowerBound,hi=cfg.reactantRangeNM.upperBound
        var xs=[lo]
        xs.append(contentsOf:profile.map(\.coordinateNM).filter{$0>lo && $0<hi});xs.append(hi)
        var integral=0.0
        for i in 1..<xs.count {
            let a=xs[i-1],b=xs[i]
            let ya=exp(-beta*(try energy(a))),yb=exp(-beta*(try energy(b)))
            integral+=0.5*(ya+yb)*(b-a)
        }
        let surface=exp(-beta*(try energy(cfg.dividingSurfaceNM)))
        let ratio=surface/integral
        guard integral.isFinite,integral>0,surface.isFinite,surface>0,ratio.isFinite,ratio>0 else {
            throw VivoChemistryError.convergence("qualified PMF reactant integral or dividing-surface density")
        }
        return ratio
    }

    private static func inverseMassMetricPerDa(coordinate:VivoQMMMReactionCoordinate,system:VivoClassicalSystem)throws->Double {
        let resolved=try VivoQMMMResolvedCoordinate(source:coordinate,system:system)
        var value=0.0
        for particle in resolved.particleIndices {
            let mass=system.particles[Int(particle)].massDa
            guard mass.isFinite,mass>0 else { throw VivoChemistryError.invalid("reaction-coordinate atom has no positive physical mass") }
            value+=1/mass
        }
        guard value.isFinite,value>0 else { throw VivoChemistryError.invalid("reaction-coordinate inverse mass metric") }
        return value
    }

    private static func binding(connectivity:VivoReactionConnectivityResult,coordinate:VivoQMMMReactionCoordinate,
                                reactantEndpointIdentifier:String,productEndpointIdentifier:String)throws->VivoQMMMReactionConnectivityBinding {
        try VivoReactionConnectivity.validate(connectivity,request:connectivity.request)
        let endpoints=Set(connectivity.request.endpoints.map(\.identifier))
        guard reactantEndpointIdentifier != productEndpointIdentifier,endpoints==Set([reactantEndpointIdentifier,productEndpointIdentifier]) else {
            throw VivoChemistryError.invalid("QM/MM PMF endpoints differ from the qualified mapped reaction")
        }
        let mapped=Set(connectivity.request.saddle.request.model.system.nuclei.compactMap(\.structureAtomIndex))
        guard Set(coordinate.atomIndices).isSubset(of:mapped) else {
            throw VivoChemistryError.invalid("QM/MM reaction coordinate contains atoms absent from the qualified saddle mapping")
        }
        let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(connectivity))
        return .init(connectivityFingerprint:id,reactantEndpointIdentifier:reactantEndpointIdentifier,
                     productEndpointIdentifier:productEndpointIdentifier,mappedReactionAtomIndices:coordinate.atomIndices)
    }

    /// Builds provenance from the same system/provider/configuration used by the
    /// production runner and from a deterministically revalidated mapped reaction.
    /// A PMF cannot be relabelled as another environment, Hamiltonian or reaction.
    public static func qualify(_ result:VivoQMMMActivationFreeEnergyResult,
                               sampling:VivoQMMMFreeEnergyRunRequest,
                               system:VivoClassicalSystem,
                               baseProvider:VivoMDCandidateForceProvider,
                               reactionConnectivity:VivoReactionConnectivityResult,
                               reactantEndpointIdentifier:String,
                               productEndpointIdentifier:String,
                               chemicalState:String,
                               environment:VivoQMMMFreeEnergyEnvironment,
                               environmentIdentifier:String,
                               methodDescription:String) throws -> VivoQMMMQualifiedActivationFreeEnergy {
        let systemID=try system.fingerprint()
        guard result.coordinate==sampling.coordinate,result.temperatureK==sampling.dynamics.targetTemperatureK,
              result.traces.map(\.window)==sampling.windows,result.traces.map(\.randomSeed)==sampling.randomSeeds,
              baseProvider.retainedSystemFingerprint==systemID else {
            throw VivoChemistryError.invalid("PMF result differs from sampling/provider identity")
        }
        let dynamicsID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(sampling.dynamics))
        let reaction=try binding(connectivity:reactionConnectivity,coordinate:result.coordinate,
                                 reactantEndpointIdentifier:reactantEndpointIdentifier,productEndpointIdentifier:productEndpointIdentifier)
        let provenance=VivoQMMMFreeEnergyProvenance(structureFingerprint:system.structureFingerprint,systemFingerprint:systemID,
            baseProviderFingerprint:baseProvider.fingerprint,dynamicsFingerprint:dynamicsID,chemicalState:chemicalState,
            environment:environment,environmentIdentifier:environmentIdentifier,methodDescription:methodDescription,
            reactionConnectivity:reaction)
        let flux=VivoQMMMFluxNormalization(surfaceToReactantDensityPerNM:try surfaceToReactantDensityPerNM(result),
                                           inverseMassMetricPerDa:try inverseMassMetricPerDa(coordinate:result.coordinate,system:system))
        return try .init(analysis:result,provenance:provenance,fluxNormalization:flux)
    }
}