import Foundation

public enum VivoQMMMReactionCoordinateKind: String, Codable, Sendable {
    case distance
    case distanceDifference
}

/// Atom indices always address the canonical prepared molecular structure. The
/// coordinate is resolved to physical classical particles once and is therefore
/// independent of virtual-site particle ordering.
public struct VivoQMMMReactionCoordinate: Codable, Sendable, Equatable {
    public var identifier: String
    public var kind: VivoQMMMReactionCoordinateKind
    public var atomIndices: [UInt32]
    public init(identifier: String,kind: VivoQMMMReactionCoordinateKind,atomIndices: [UInt32]) {
        self.identifier=identifier;self.kind=kind;self.atomIndices=atomIndices
    }
    public func validate() throws {
        let expected = kind == .distance ? 2 : 4
        guard !identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,atomIndices.count==expected else {
            throw VivoChemistryError.invalid("reaction-coordinate identity or atom arity")
        }
        switch kind {
        case .distance:
            guard atomIndices[0] != atomIndices[1] else {
                throw VivoChemistryError.invalid("distance coordinate requires two distinct atoms")
            }
        case .distanceDifference:
            let first=Set([atomIndices[0],atomIndices[1]]),second=Set([atomIndices[2],atomIndices[3]])
            guard first.count==2,second.count==2,first != second else {
                throw VivoChemistryError.invalid("distance-difference coordinate requires two nonidentical atom pairs")
            }
        }
    }
}

public struct VivoQMMMResolvedCoordinate: Sendable, Equatable {
    public let source: VivoQMMMReactionCoordinate
    public let particleIndices: [UInt32]
    public init(source: VivoQMMMReactionCoordinate,system: VivoClassicalSystem) throws {
        try source.validate();try VivoClassicalSystemValidator.validate(system)
        var byAtom:[UInt32:UInt32]=[:]
        for particle in system.particles where particle.role == .atom {
            guard let atom=particle.atomIndex,byAtom[atom]==nil else {
                throw VivoChemistryError.invalid("reaction coordinate requires a unique physical particle per atom")
            }
            byAtom[atom]=particle.index
        }
        var particles:[UInt32]=[]
        for atom in source.atomIndices {
            guard let particle=byAtom[atom] else { throw VivoChemistryError.invalid("reaction-coordinate atom \(atom) is not a physical particle") }
            particles.append(particle)
        }
        self.source=source;particleIndices=particles
    }
    /// Returns xi in nm and d(xi)/d(r_particle) in 1. Coordinates use the same
    /// exact minimum-image implementation as the retained Hamiltonian.
    public func evaluate(_ geometry: VivoMDCandidateGeometry) throws -> (valueNM:Double,gradients:[UInt32:VivoVector3D]) {
        func displacement(_ a:UInt32,_ b:UInt32)throws->VivoVector3D {
            guard Int(a)<geometry.particlePositionsNM.count,Int(b)<geometry.particlePositionsNM.count else {
                throw VivoChemistryError.invalid("reaction-coordinate particle shape")
            }
            let d=geometry.particlePositionsNM[Int(a)]-geometry.particlePositionsNM[Int(b)]
            return try geometry.periodicCell?.minimumImage(d) ?? d
        }
        func bond(_ a:UInt32,_ b:UInt32)throws->(Double,VivoVector3D) {
            let d=try displacement(a,b),r=d.norm
            guard r.isFinite,r>1e-10 else { throw VivoChemistryError.invalid("reaction-coordinate zero distance") }
            return (r,d/r)
        }
        switch source.kind {
        case .distance:
            let (r,u)=try bond(particleIndices[0],particleIndices[1])
            return (r,[particleIndices[0]:u,particleIndices[1]:u*(-1)])
        case .distanceDifference:
            let (a,u)=try bond(particleIndices[0],particleIndices[1])
            let (b,v)=try bond(particleIndices[2],particleIndices[3])
            var gradients:[UInt32:VivoVector3D]=[:]
            func accumulate(_ particle:UInt32,_ contribution:VivoVector3D) {
                gradients[particle]=(gradients[particle] ?? .zero)+contribution
            }
            accumulate(particleIndices[0],u);accumulate(particleIndices[1],u*(-1))
            accumulate(particleIndices[2],v*(-1));accumulate(particleIndices[3],v)
            return (a-b,gradients)
        }
    }
}

public struct VivoQMMMUmbrellaWindow: Codable, Sendable, Equatable {
    public var identifier:String
    public var centerNM:Double
    public var forceConstantKJPerMolNM2:Double
    public init(identifier:String,centerNM:Double,forceConstantKJPerMolNM2:Double) {
        self.identifier=identifier;self.centerNM=centerNM;self.forceConstantKJPerMolNM2=forceConstantKJPerMolNM2
    }
    public func validate() throws {
        guard !identifier.isEmpty,centerNM.isFinite,forceConstantKJPerMolNM2.isFinite,forceConstantKJPerMolNM2>0 else {
            throw VivoChemistryError.invalid("umbrella-window identity, center or force constant")
        }
    }
    public func biasKJPerMol(_ coordinateNM:Double)->Double {
        0.5*forceConstantKJPerMolNM2*pow(coordinateNM-centerNM,2)
    }
}

public enum VivoQMMMUmbrellaBias {
    /// Composes the restraint with the existing complete BO provider. No QM/MM,
    /// polarization, PME or adaptive-region implementation is duplicated here.
    /// The wrapper is NVT/NVE only because cell derivatives of a minimum-image
    /// reaction coordinate are deliberately not approximated.
    public static func provider(base: VivoMDCandidateForceProvider,system: VivoClassicalSystem,
                                coordinate: VivoQMMMReactionCoordinate,window: VivoQMMMUmbrellaWindow) throws -> VivoMDCandidateForceProvider {
        try window.validate();let resolved=try VivoQMMMResolvedCoordinate(source:coordinate,system:system)
        struct Identity:Codable { let schema:String;let base:VivoFingerprint;let system:VivoFingerprint;let coordinate:VivoQMMMReactionCoordinate;let window:VivoQMMMUmbrellaWindow }
        let systemID=try system.fingerprint()
        guard base.retainedSystemFingerprint==systemID else { throw VivoChemistryError.invalid("umbrella/base retained-system mismatch") }
        let identity=Identity(schema:"numivivo.org/qmmm-umbrella-provider/v1",base:base.fingerprint,system:systemID,coordinate:coordinate,window:window)
        let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity))
        return try VivoMDCandidateForceProvider(fingerprint:id,retainedSystemFingerprint:systemID,boundary:base.boundary,
            supportsCellMoves:false,maximumAcceptedResidual:base.maximumAcceptedResidual,
            molecularConnectivitySystem:base.molecularConnectivitySystem,polarizationModelFingerprint:base.polarizationModelFingerprint) { geometry in
            let raw=try await base.evaluate(geometry);try raw.validate(geometry:geometry,provider:base,system:system)
            let reaction=try resolved.evaluate(geometry),delta=reaction.valueNM-window.centerNM
            let slope=window.forceConstantKJPerMolNM2*delta
            var force=raw.physicalParticleForcesKJPerMolNM
            for (particle,gradient) in reaction.gradients { force[Int(particle)]=force[Int(particle)]-gradient*slope }
            return try VivoMDCandidateForceEvaluation(providerFingerprint:id,geometry:geometry,
                additionalEnergyKJPerMol:raw.additionalEnergyKJPerMol+window.biasKJPerMol(reaction.valueNM),
                physicalParticleForcesKJPerMolNM:force,derivativeMethod:raw.derivativeMethod+"; analytic harmonic umbrella",
                convergenceResidual:raw.convergenceResidual,requiredResidual:raw.requiredResidual,
                additionalAffineStrainDerivativeKJPerMol:nil)
        }
    }
}

public struct VivoQMMMUmbrellaTrace: Codable, Sendable, Equatable {
    public let window:VivoQMMMUmbrellaWindow
    public let randomSeed:UInt64
    public let coordinateNM:[Double]
    public let potentialEnergyKJPerMol:[Double]
    public init(window:VivoQMMMUmbrellaWindow,randomSeed:UInt64,coordinateNM:[Double],potentialEnergyKJPerMol:[Double]) {
        self.window=window;self.randomSeed=randomSeed;self.coordinateNM=coordinateNM;self.potentialEnergyKJPerMol=potentialEnergyKJPerMol
    }
}

public struct VivoQMMMFreeEnergyAnalysisConfiguration: Codable, Sendable, Equatable {
    public var bins:Int
    public var kernelBandwidthNM:Double
    public var reactantRangeNM:ClosedRange<Double>
    public var dividingSurfaceNM:Double
    public var maximumMBARIterations:Int
    public var mbarTolerance:Double
    public var minimumDecorrelatedSamplesPerWindow:Int
    public var minimumAdjacentOverlap:Double
    public var maximumAutocorrelationLag:Int
    public init(bins:Int=160,kernelBandwidthNM:Double=0.015,reactantRangeNM:ClosedRange<Double>,dividingSurfaceNM:Double,
                maximumMBARIterations:Int=10000,mbarTolerance:Double=1e-10,minimumDecorrelatedSamplesPerWindow:Int=100,
                minimumAdjacentOverlap:Double=0.01,maximumAutocorrelationLag:Int=4096) {
        self.bins=bins;self.kernelBandwidthNM=kernelBandwidthNM;self.reactantRangeNM=reactantRangeNM
        self.dividingSurfaceNM=dividingSurfaceNM;self.maximumMBARIterations=maximumMBARIterations;self.mbarTolerance=mbarTolerance
        self.minimumDecorrelatedSamplesPerWindow=minimumDecorrelatedSamplesPerWindow;self.minimumAdjacentOverlap=minimumAdjacentOverlap
        self.maximumAutocorrelationLag=maximumAutocorrelationLag
    }
    public func validate() throws {
        guard (32...4096).contains(bins),kernelBandwidthNM.isFinite,kernelBandwidthNM>0,
              reactantRangeNM.lowerBound.isFinite,reactantRangeNM.upperBound.isFinite,
              reactantRangeNM.lowerBound<reactantRangeNM.upperBound,dividingSurfaceNM.isFinite,
              (10...100000).contains(maximumMBARIterations),mbarTolerance.isFinite,mbarTolerance>0,mbarTolerance<1e-3,
              minimumDecorrelatedSamplesPerWindow>=20,(0..<0.5).contains(minimumAdjacentOverlap),maximumAutocorrelationLag>=3 else {
            throw VivoChemistryError.invalid("QM/MM free-energy analysis configuration")
        }
    }
}

public struct VivoQMMMFreeEnergyProfilePoint: Codable, Sendable, Equatable {
    public let coordinateNM:Double
    public let relativeFreeEnergyKJPerMol:Double
    public let localEffectiveSamples:Double
}
public struct VivoQMMMWindowDiagnostics: Codable, Sendable, Equatable {
    public let identifier:String
    public let rawSamples:Int
    public let decorrelationStride:Int
    public let retainedSamples:Int
    public let adjacentOverlapLeft:Double?
    public let adjacentOverlapRight:Double?
}
public struct VivoQMMMActivationFreeEnergyResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-activation-free-energy/v1"
    public let schema:String
    public let coordinate:VivoQMMMReactionCoordinate
    public let temperatureK:Double
    public let traces:[VivoQMMMUmbrellaTrace]
    public let configuration:VivoQMMMFreeEnergyAnalysisConfiguration
    public let profile:[VivoQMMMFreeEnergyProfilePoint]
    public let reactantCoordinateNM:Double
    public let dividingSurfaceNM:Double
    /// PMF profile height at the declared dividing surface relative to the
    /// lowest PMF point in the declared reactant range. This is a diagnostic,
    /// not by itself the flux-normalized kinetic activation free energy.
    public let activationFreeEnergyKJPerMol:Double
    /// Conditional sampling uncertainty from local reweighted effective counts.
    /// It excludes Hamiltonian/model/microstate/transmission uncertainty.
    public let conditionalStandardDeviationKJPerMol:Double
    public let mbarIterations:Int
    public let mbarResidual:Double
    public let diagnostics:[VivoQMMMWindowDiagnostics]
    public let converged:Bool
    public let issues:[String]
    public let interpretation:String
}

public enum VivoQMMMFreeEnergy {
    public static let interpretation="Umbrella-sampled PMF under one fingerprinted QM/MM Born-Oppenheimer Hamiltonian, reconstructed by unbinned MBAR after conservative autocorrelation thinning. PMF profile height is diagnostic; kinetic conversion requires reactant-basin normalization and a coordinate mass metric. Conditional uncertainty covers finite reweighted samples only; protonation/conformer populations, Hamiltonian accuracy, tunnelling and dynamical recrossing remain separate."
    private static let gasConstantKJ=0.00831446261815324

    public static func analyze(coordinate:VivoQMMMReactionCoordinate,temperatureK:Double,traces:[VivoQMMMUmbrellaTrace],
                               configuration cfg:VivoQMMMFreeEnergyAnalysisConfiguration) throws -> VivoQMMMActivationFreeEnergyResult {
        try coordinate.validate();try cfg.validate()
        guard temperatureK.isFinite,temperatureK>0,traces.count>=2,traces.count<=256,
              Set(traces.map{$0.window.identifier}).count==traces.count,
              Set(traces.map{$0.randomSeed}).count==traces.count else {
            throw VivoChemistryError.invalid("free-energy temperature or independent window identities")
        }
        let ordered=traces.sorted{$0.window.centerNM<$1.window.centerNM}
        guard ordered==traces else { throw VivoChemistryError.invalid("umbrella windows must be supplied in increasing center order") }
        var samples:[[Double]]=[],strides:[Int]=[],issues:[String]=[]
        for trace in traces {
            try trace.window.validate()
            guard trace.coordinateNM.count==trace.potentialEnergyKJPerMol.count,trace.coordinateNM.count>=cfg.minimumDecorrelatedSamplesPerWindow,
                  trace.coordinateNM.allSatisfy(\.isFinite),trace.potentialEnergyKJPerMol.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("umbrella trace shape or finite values")
            }
            let stride=autocorrelationStride(trace.coordinateNM,maximumLag:cfg.maximumAutocorrelationLag)
            let kept=Swift.stride(from:0,to:trace.coordinateNM.count,by:stride).map{trace.coordinateNM[$0]}
            if kept.count<cfg.minimumDecorrelatedSamplesPerWindow { issues.append("\(trace.window.identifier): insufficient decorrelated samples") }
            samples.append(kept);strides.append(stride)
        }
        let counts=samples.map(\.count),all=samples.flatMap{$0},source=samples.enumerated().flatMap{Array(repeating:$0.offset,count:$0.element.count)}
        guard let minimum=all.min(),let maximum=all.max(),minimum<maximum,
              cfg.dividingSurfaceNM>=minimum,cfg.dividingSurfaceNM<=maximum,
              cfg.reactantRangeNM.overlaps(minimum...maximum) else {
            throw VivoChemistryError.invalid("reaction PMF does not span reactant and dividing-surface coordinates")
        }
        let beta=1/(gasConstantKJ*temperatureK),k=traces.count
        func logSumExp(_ values:[Double])->Double {
            let m=values.max()!;return m+log(values.reduce(0){$0+exp($1-m)})
        }
        func reducedBias(_ window:Int,_ x:Double)->Double { beta*traces[window].window.biasKJPerMol(x) }
        var f=[Double](repeating:0,count:k),residual=Double.infinity,iterations=0
        for iteration in 1...cfg.maximumMBARIterations {
            var next=[Double](repeating:0,count:k)
            for i in 0..<k {
                var terms:[Double]=[];terms.reserveCapacity(all.count)
                for x in all {
                    let denominator=logSumExp((0..<k).map{log(Double(counts[$0]))+f[$0]-reducedBias($0,x)})
                    terms.append(-reducedBias(i,x)-denominator)
                }
                next[i] = -logSumExp(terms)
            }
            let shift=next[0];for i in 0..<k { next[i]-=shift }
            residual=zip(next,f).map{abs($0-$1)}.max() ?? 0;f=next;iterations=iteration
            if residual<=cfg.mbarTolerance { break }
        }
        if residual>cfg.mbarTolerance { issues.append("MBAR fixed-point residual exceeds tolerance") }
        var logDen=[Double]();logDen.reserveCapacity(all.count)
        for x in all { logDen.append(logSumExp((0..<k).map{log(Double(counts[$0]))+f[$0]-reducedBias($0,x)})) }
        let rawLogWeights=logDen.map{-1*$0},weightNorm=logSumExp(rawLogWeights)
        let weights=rawLogWeights.map{exp($0-weightNorm)}
        var overlap=Array(repeating:Array(repeating:0.0,count:k),count:k)
        for n in all.indices {
            let responsibilities=(0..<k).map{exp(log(Double(counts[$0]))+f[$0]-reducedBias($0,all[n])-logDen[n])}
            for j in 0..<k { overlap[source[n]][j]+=responsibilities[j]/Double(counts[source[n]]) }
        }
        for i in 0..<(k-1) {
            if min(overlap[i][i+1],overlap[i+1][i])<cfg.minimumAdjacentOverlap {
                issues.append("\(traces[i].window.identifier)↔\(traces[i+1].window.identifier): insufficient adjacent phase-space overlap")
            }
        }
        let width=(maximum-minimum)/Double(cfg.bins-1),band=cfg.kernelBandwidthNM
        func density(_ x:Double)->(Double,Double) {
            var terms=[Double]();terms.reserveCapacity(all.count)
            for n in all.indices { let z=(all[n]-x)/band;terms.append(weights[n]*exp(-0.5*z*z)) }
            let sum=terms.reduce(0,+)/(band*sqrt(2*Double.pi))
            let square=terms.reduce(0){$0+$1*$1}
            let ess=square>0 ? pow(terms.reduce(0,+),2)/square:0
            return (sum,ess)
        }
        var rawProfile:[(Double,Double,Double)]=[]
        for bin in 0..<cfg.bins {
            let x=minimum+Double(bin)*width,(p,ess)=density(x)
            guard p.isFinite,p>0,ess.isFinite else { throw VivoChemistryError.convergence("nonfinite PMF density") }
            rawProfile.append((x,-log(p)/beta,ess))
        }
        let offset=rawProfile.map{$0.1}.min()!,profile=rawProfile.map{VivoQMMMFreeEnergyProfilePoint(coordinateNM:$0.0,relativeFreeEnergyKJPerMol:$0.1-offset,localEffectiveSamples:$0.2)}
        let reactantCandidates=rawProfile.filter{cfg.reactantRangeNM.contains($0.0)}
        guard let reactant=reactantCandidates.min(by:{$0.1<$1.1}) else { throw VivoChemistryError.invalid("reactant PMF basin has no grid points") }
        let ts=rawProfile.min(by:{abs($0.0-cfg.dividingSurfaceNM)<abs($1.0-cfg.dividingSurfaceNM)})!
        let barrier=ts.1-reactant.1
        guard barrier.isFinite,barrier>=0 else { issues.append("declared dividing surface is not above the reactant PMF basin")
            return .init(schema:VivoQMMMActivationFreeEnergyResult.schema,coordinate:coordinate,temperatureK:temperatureK,traces:traces,
                configuration:cfg,profile:profile,reactantCoordinateNM:reactant.0,dividingSurfaceNM:cfg.dividingSurfaceNM,
                activationFreeEnergyKJPerMol:max(0,barrier),conditionalStandardDeviationKJPerMol:Double.greatestFiniteMagnitude,
                mbarIterations:iterations,mbarResidual:residual,diagnostics:[],converged:false,issues:issues,interpretation:interpretation) }
        if reactant.2<Double(cfg.minimumDecorrelatedSamplesPerWindow) { issues.append("reactant basin has low local effective sample count") }
        if ts.2<Double(cfg.minimumDecorrelatedSamplesPerWindow) { issues.append("dividing surface has low local effective sample count") }
        let sd=gasConstantKJ*temperatureK*sqrt(1/max(1,reactant.2)+1/max(1,ts.2))
        var diagnostics:[VivoQMMMWindowDiagnostics]=[]
        for i in 0..<k {
            diagnostics.append(.init(identifier:traces[i].window.identifier,rawSamples:traces[i].coordinateNM.count,
                decorrelationStride:strides[i],retainedSamples:counts[i],adjacentOverlapLeft:i>0 ? overlap[i][i-1]:nil,
                adjacentOverlapRight:i+1<k ? overlap[i][i+1]:nil))
        }
        return .init(schema:VivoQMMMActivationFreeEnergyResult.schema,coordinate:coordinate,temperatureK:temperatureK,traces:traces,
            configuration:cfg,profile:profile,reactantCoordinateNM:reactant.0,dividingSurfaceNM:cfg.dividingSurfaceNM,
            activationFreeEnergyKJPerMol:barrier,conditionalStandardDeviationKJPerMol:sd,mbarIterations:iterations,mbarResidual:residual,
            diagnostics:diagnostics,converged:issues.isEmpty,issues:issues,interpretation:interpretation)
    }

    public static func validate(_ result:VivoQMMMActivationFreeEnergyResult) throws {
        let rebuilt=try analyze(coordinate:result.coordinate,temperatureK:result.temperatureK,traces:result.traces,configuration:result.configuration)
        guard rebuilt==result else { throw VivoChemistryError.invalid("QM/MM free-energy result does not reconstruct from retained traces") }
    }

    private static func autocorrelationStride(_ x:[Double],maximumLag:Int)->Int {
        guard x.count>3 else { return 1 }
        let mean=x.reduce(0,+)/Double(x.count),variance=x.reduce(0){$0+pow($1-mean,2)}/Double(x.count)
        guard variance>0,variance.isFinite else { return x.count }
        var tau=1.0
        for lag in 1...min(maximumLag,x.count-2) {
            var covariance=0.0
            for i in 0..<(x.count-lag) { covariance+=(x[i]-mean)*(x[i+lag]-mean) }
            let rho=covariance/Double(x.count-lag)/variance
            if !rho.isFinite || rho<=0 { break };tau+=2*rho
        }
        return max(1,min(x.count,Int(ceil(tau))))
    }
}

public struct VivoQMMMFreeEnergyRunRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-run/v1"
    public var schema:String
    public var coordinate:VivoQMMMReactionCoordinate
    public var windows:[VivoQMMMUmbrellaWindow]
    public var initialStates:[VivoClassicalInitialState]
    public var randomSeeds:[UInt64]
    public var dynamics:VivoMDConfiguration
    public var equilibrationSteps:UInt64
    public var productionSteps:UInt64
    public var sampleEvery:UInt64
    public var analysis:VivoQMMMFreeEnergyAnalysisConfiguration
    public init(coordinate:VivoQMMMReactionCoordinate,windows:[VivoQMMMUmbrellaWindow],initialStates:[VivoClassicalInitialState],
                randomSeeds:[UInt64],dynamics:VivoMDConfiguration,equilibrationSteps:UInt64,productionSteps:UInt64,sampleEvery:UInt64,
                analysis:VivoQMMMFreeEnergyAnalysisConfiguration) {
        schema=Self.schema;self.coordinate=coordinate;self.windows=windows;self.initialStates=initialStates;self.randomSeeds=randomSeeds
        self.dynamics=dynamics;self.equilibrationSteps=equilibrationSteps;self.productionSteps=productionSteps;self.sampleEvery=sampleEvery;self.analysis=analysis
    }
}

public enum VivoQMMMFreeEnergyRunner {
    /// Runs each window through the same production Metal integrator and composed
    /// BO provider. Windows are sequential to bound Apple unified-memory pressure.
    public static func run(_ request:VivoQMMMFreeEnergyRunRequest,system:VivoClassicalSystem,
                           baseProvider:VivoMDCandidateForceProvider) async throws -> VivoQMMMActivationFreeEnergyResult {
        try request.coordinate.validate();try request.analysis.validate();try request.dynamics.validate()
        guard request.schema==VivoQMMMFreeEnergyRunRequest.schema,request.windows.count>=2,
              request.windows.count==request.initialStates.count,request.windows.count==request.randomSeeds.count,
              Set(request.randomSeeds).count==request.randomSeeds.count,request.productionSteps>0,request.sampleEvery>0,
              request.productionSteps%request.sampleEvery==0,request.dynamics.ensemble == .nvt,
              request.dynamics.thermostat == .langevinMiddle,request.dynamics.targetTemperatureK != nil,
              request.dynamics.frictionPerPS.map({$0>0}) == true else {
            throw VivoChemistryError.invalid("umbrella execution requires independent windows and Langevin NVT")
        }
        let systemID=try system.fingerprint();guard baseProvider.retainedSystemFingerprint==systemID else { throw VivoChemistryError.invalid("free-energy provider/system mismatch") }
        let coordinate=try VivoQMMMResolvedCoordinate(source:request.coordinate,system:system)
        var traces:[VivoQMMMUmbrellaTrace]=[]
        for i in request.windows.indices {
            try request.windows[i].validate();try request.initialStates[i].validate(particleCount:system.particles.count)
            guard request.initialStates[i].systemFingerprint==systemID else { throw VivoChemistryError.invalid("umbrella initial state/system mismatch") }
            var dynamics=request.dynamics;dynamics.randomSeed=request.randomSeeds[i]
            let biased=try VivoQMMMUmbrellaBias.provider(base:baseProvider,system:system,coordinate:request.coordinate,window:request.windows[i])
            let runtime=try await VivoMDMetalRuntime.make(system:system,initialState:request.initialStates[i],configuration:dynamics,forceProvider:biased)
            _ = try await runtime.thermalize(temperatureK:dynamics.targetTemperatureK!,seed:request.randomSeeds[i]^0x554d4252454c4c41)
            for _ in 0..<request.equilibrationSteps {
                try Task.checkCancellation();guard try await runtime.step().committed else { throw VivoChemistryError.convergence("umbrella equilibration candidate rejected") }
            }
            var values:[Double]=[],energies:[Double]=[]
            values.reserveCapacity(Int(request.productionSteps/request.sampleEvery));energies.reserveCapacity(values.capacity)
            for step in 1...request.productionSteps {
                try Task.checkCancellation();guard try await runtime.step().committed else { throw VivoChemistryError.convergence("umbrella production candidate rejected") }
                if step%request.sampleEvery==0 {
                    let sample=try await runtime.sample(includeObservables:true),geometry=try VivoMDCandidateGeometry(particlePositionsNM:sample.state.positionsNM,periodicCell:sample.state.periodicCell)
                    values.append(try coordinate.evaluate(geometry).valueNM)
                    guard let energy=sample.observables?.potentialEnergyKJPerMol else { throw VivoChemistryError.invalid("umbrella sample omitted potential energy") }
                    energies.append(energy)
                }
            }
            traces.append(.init(window:request.windows[i],randomSeed:request.randomSeeds[i],coordinateNM:values,potentialEnergyKJPerMol:energies))
        }
        return try VivoQMMMFreeEnergy.analyze(coordinate:request.coordinate,temperatureK:request.dynamics.targetTemperatureK!,traces:traces,configuration:request.analysis)
    }
}
