import Foundation

public enum VivoMolecularObservableKind: Codable, Sendable, Equatable {
    case distance(atomA: UInt32, atomB: UInt32)
    case distanceDifference(a: UInt32, b: UInt32, c: UInt32, d: UInt32)
    case torsionCosine(a: UInt32, b: UInt32, c: UInt32, d: UInt32)
    case torsionSine(a: UInt32, b: UInt32, c: UInt32, d: UInt32)
    case potentialEnergyKJPerMol
    case temperatureK
    case volumeNM3
}
public struct VivoMolecularObservable: Codable, Sendable, Equatable {
    public var identifier: String
    public var kind: VivoMolecularObservableKind
    /// Absolute precision in the observable's declared units, not a percentage.
    public var maximumMeanStandardError: Double
    public init(identifier: String, kind: VivoMolecularObservableKind, maximumMeanStandardError: Double) {
        self.identifier=identifier; self.kind=kind; self.maximumMeanStandardError=maximumMeanStandardError
    }
    public func validate() throws {
        guard !identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,
              maximumMeanStandardError.isFinite, maximumMeanStandardError>0 else {
            throw VivoArtifactValidationError.invalid("molecular observable identity or precision")
        }
    }
    /// Geometry references are STRUCTURE atoms. The explicit atom->particle map
    /// avoids confusing virtual-site particles with physical atom identities.
    public func measure(positionsNM: [VivoVector3D], periodicCell: VivoPeriodicCell?,
                        particleByAtom: [Int], potentialEnergyKJPerMol: Double? = nil,
                        temperatureK: Double? = nil) throws -> Double {
        try validate()
        func position(_ atom: UInt32) throws -> VivoVector3D {
            guard Int(atom)<particleByAtom.count, particleByAtom[Int(atom)]>=0,
                  particleByAtom[Int(atom)]<positionsNM.count else {
                throw VivoArtifactValidationError.invalid("observable atom has no mapped physical particle")
            }
            let p=positionsNM[particleByAtom[Int(atom)]]
            guard p.isFinite else { throw VivoArtifactValidationError.invalid("observable position is nonfinite") }
            return p
        }
        func delta(_ a: UInt32,_ b: UInt32) throws -> VivoVector3D {
            guard a != b else { throw VivoArtifactValidationError.invalid("observable displacement has identical atoms") }
            let d=try position(a)-position(b)
            return try periodicCell?.minimumImage(d) ?? d
        }
        let value: Double
        switch kind {
        case .distance(let a,let b): value=try delta(a,b).norm
        case .distanceDifference(let a,let b,let c,let d): value=try delta(a,b).norm-delta(c,d).norm
        case .torsionCosine(let a,let b,let c,let d),.torsionSine(let a,let b,let c,let d):
            guard Set([a,b,c,d]).count==4 else { throw VivoArtifactValidationError.invalid("torsion observable needs four distinct atoms") }
            let x=try delta(b,a),y=try delta(c,b),z=try delta(d,c)
            let n=x.cross(y),m=y.cross(z)
            guard n.norm>1e-14,m.norm>1e-14,y.norm>1e-14 else {
                throw VivoArtifactValidationError.invalid("torsion observable is undefined for collinear atoms")
            }
            if case .torsionCosine = kind { value=n.dot(m)/(n.norm*m.norm) }
            else { value=(y/y.norm).dot(n.cross(m))/(n.norm*m.norm) }
        case .potentialEnergyKJPerMol:
            guard let energy=potentialEnergyKJPerMol else { throw VivoArtifactValidationError.unresolved("sampling frame has no potential energy") }
            value=energy
        case .temperatureK:
            guard let temperature=temperatureK else { throw VivoArtifactValidationError.unresolved("sampling frame has no temperature") }
            value=temperature
        case .volumeNM3:
            guard let cell=periodicCell,cell.isValid else { throw VivoArtifactValidationError.unresolved("volume observable requires a valid periodic cell") }
            value=cell.volumeNM3
        }
        guard value.isFinite else { throw VivoArtifactValidationError.invalid("nonfinite molecular observable") }
        return value
    }
}

public struct VivoMolecularSamplingConfiguration: Codable, Sendable, Equatable {
    public var discardFramesPerReplica: Int
    public var minimumRetainedFramesPerReplica: Int
    public var minimumReplicas: Int
    public var maximumRHat: Double
    public var minimumEffectiveSamplesPerReplica: Double
    public var maximumAutocorrelationLag: Int
    public var maximumAutocorrelationWork: Int
    public var maximumTotalSamples: Int
    public init(discardFramesPerReplica: Int = 0, minimumRetainedFramesPerReplica: Int = 256,
                minimumReplicas: Int = 4, maximumRHat: Double = 1.01,
                minimumEffectiveSamplesPerReplica: Double = 100,
                maximumAutocorrelationLag: Int = 4096, maximumAutocorrelationWork: Int = 100_000_000,
                maximumTotalSamples: Int = 2_000_000) {
        self.discardFramesPerReplica=discardFramesPerReplica
        self.minimumRetainedFramesPerReplica=minimumRetainedFramesPerReplica; self.minimumReplicas=minimumReplicas
        self.maximumRHat=maximumRHat; self.minimumEffectiveSamplesPerReplica=minimumEffectiveSamplesPerReplica
        self.maximumAutocorrelationLag=maximumAutocorrelationLag; self.maximumAutocorrelationWork=maximumAutocorrelationWork
        self.maximumTotalSamples=maximumTotalSamples
    }
    public func validate() throws {
        guard discardFramesPerReplica>=0,minimumRetainedFramesPerReplica>=64,(2...64).contains(minimumReplicas),
              maximumRHat.isFinite,maximumRHat>=1,maximumRHat<=1.1,
              minimumEffectiveSamplesPerReplica.isFinite,minimumEffectiveSamplesPerReplica>=20,
              maximumAutocorrelationLag>=3,maximumAutocorrelationWork>0,maximumTotalSamples>0 else {
            throw VivoArtifactValidationError.invalid("molecular sampling convergence configuration")
        }
    }
}

/// Scalar observations retain the verified trajectory/checkpoint identity without
/// retaining every protein coordinate in RAM. Values are observable-major.
public struct VivoMolecularReplicaSeries: Codable, Sendable, Equatable {
    public var identifier: String
    public var sourceFingerprint: VivoFingerprint
    public var configuration: VivoMDConfiguration
    public var steps: [UInt64]
    public var timesPS: [Double]
    public var valuesByObservable: [[Double]]
    public init(identifier: String, sourceFingerprint: VivoFingerprint, configuration: VivoMDConfiguration,
                steps: [UInt64], timesPS: [Double], valuesByObservable: [[Double]]) {
        self.identifier=identifier;self.sourceFingerprint=sourceFingerprint;self.configuration=configuration
        self.steps=steps;self.timesPS=timesPS;self.valuesByObservable=valuesByObservable
    }
    public static func fromTrajectory(_ trajectory: VivoTrajectory, identifier: String,
                                      configuration: VivoMDConfiguration,
                                      structureFingerprint: VivoFingerprint,
                                      observables: [VivoMolecularObservable]) throws -> Self {
        try trajectory.validate()
        guard trajectory.structureFingerprint == structureFingerprint else {
            throw VivoArtifactValidationError.incompatible("sampling trajectory differs from prepared structure")
        }
        let map=Array(0..<Int(trajectory.atomCount))
        var values=Array(repeating:[Double](),count:observables.count)
        for frame in trajectory.frames {
            for i in observables.indices {
                values[i].append(try observables[i].measure(positionsNM:frame.positionsNM,periodicCell:frame.periodicCell,
                    particleByAtom:map,potentialEnergyKJPerMol:frame.potentialEnergyKJPerMol,temperatureK:frame.temperatureK))
            }
        }
        return .init(identifier:identifier,sourceFingerprint:try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(trajectory)),
                     configuration:configuration,steps:trajectory.frames.map(\.step),timesPS:trajectory.frames.map(\.timePS),valuesByObservable:values)
    }
}
public struct VivoMolecularSamplingRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/molecular-sampling/v1"
    public var schema: String
    public var structureFingerprint: VivoFingerprint
    public var systemFingerprint: VivoFingerprint
    public var contextIdentifier: String
    public var observables: [VivoMolecularObservable]
    public var replicas: [VivoMolecularReplicaSeries]
    public var configuration: VivoMolecularSamplingConfiguration
    public init(structureFingerprint: VivoFingerprint, systemFingerprint: VivoFingerprint, contextIdentifier: String,
                observables: [VivoMolecularObservable], replicas: [VivoMolecularReplicaSeries],
                configuration: VivoMolecularSamplingConfiguration = .init()) {
        schema=Self.schema;self.structureFingerprint=structureFingerprint;self.systemFingerprint=systemFingerprint
        self.contextIdentifier=contextIdentifier;self.observables=observables;self.replicas=replicas;self.configuration=configuration
    }
}
public struct VivoMolecularObservableDiagnostics: Codable, Sendable, Equatable {
    public let identifier: String
    public let mean: Double
    public let meanStandardError: Double?
    public let rankNormalizedSplitRHat: Double?
    public let foldedRankNormalizedSplitRHat: Double?
    public let bulkEffectiveSamples: Double?
    public let tailEffectiveSamples: Double?
    public let autocorrelationResolved: Bool
    public let issues: [String]
}
public struct VivoMolecularSamplingResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/molecular-sampling-result/v1"
    public let schema: String
    public let request: VivoMolecularSamplingRequest
    public let requestFingerprint: VivoFingerprint
    public let retainedFramesPerReplica: Int
    public let diagnostics: [VivoMolecularObservableDiagnostics]
    public let converged: Bool
    public let issues: [String]
    public let interpretation: String
}

/// Rank-normalized/folded split-Rhat and multi-chain Geyer initial-positive,
/// initial-monotone ESS. Raw-scale ESS controls mean standard errors. Tail ESS
/// uses 5%/95% indicators. These diagnose DECLARED observables, not all molecular
/// conformations, unobserved rare events, force-field accuracy or protonation.
/// Algorithm reference: Vehtari et al., Bayesian Analysis 16 (2021), 667-718.
public enum VivoMolecularSampling {
    public static let interpretation="Replica-aware convergence of the declared molecular observables under one specified MD Hamiltonian and ensemble; not proof of global conformational exploration, accurate force fields, protonation equilibrium, reaction free energies or dynamical rates. Constant traces and unresolved autocorrelation tails do not certify convergence."

    public static func analyze(_ request: VivoMolecularSamplingRequest) throws -> VivoMolecularSamplingResult {
        let cfg=request.configuration
        try cfg.validate()
        guard request.schema==VivoMolecularSamplingRequest.schema,!request.contextIdentifier.isEmpty,
              !request.observables.isEmpty,request.observables.count<=64,
              request.replicas.count>=cfg.minimumReplicas,request.replicas.count<=64,
              Set(request.observables.map(\.identifier)).count==request.observables.count,
              Set(request.replicas.map(\.identifier)).count==request.replicas.count,
              Set(request.replicas.map(\.sourceFingerprint)).count==request.replicas.count,
              Set(request.replicas.map { $0.configuration.randomSeed }).count==request.replicas.count else {
            throw VivoArtifactValidationError.invalid("sampling schema, context, observables or independent replica identities/seeds")
        }
        for observable in request.observables { try observable.validate() }
        let count=request.replicas[0].timesPS.count
        let total=count.multipliedReportingOverflow(by:request.replicas.count)
        guard !total.overflow,total.partialValue<=cfg.maximumTotalSamples,count>=2 else {
            throw VivoArtifactValidationError.invalid("sampling frame count or capacity")
        }
        var referenceConfiguration=request.replicas[0].configuration
        referenceConfiguration.randomSeed=0
        let dt=request.replicas[0].timesPS[1]-request.replicas[0].timesPS[0]
        guard dt.isFinite,dt>0 else { throw VivoArtifactValidationError.invalid("sampling time interval") }
        for replica in request.replicas {
            try replica.configuration.validate()
            var common=replica.configuration; common.randomSeed=0
            guard !replica.identifier.isEmpty,common==referenceConfiguration,
                  common.thermostat == .langevinMiddle,common.frictionPerPS.map({$0>0}) == true,
                  replica.steps.count==count,replica.timesPS.count==count,
                  replica.valuesByObservable.count==request.observables.count,
                  replica.valuesByObservable.allSatisfy({$0.count==count && $0.allSatisfy(\.isFinite)}),
                  replica.timesPS.allSatisfy({$0.isFinite && $0>=0}) else {
                throw VivoArtifactValidationError.incompatible("replicas must share a Langevin NVT/NPT model and equal finite observation schedules")
            }
            for i in 1..<count {
                let delta=replica.timesPS[i]-replica.timesPS[i-1]
                guard replica.steps[i]>replica.steps[i-1],abs(delta-dt)<=max(1e-10,dt*1e-6) else {
                    throw VivoArtifactValidationError.invalid("molecular observations are not uniformly sampled")
                }
            }
        }
        let retained=max(0,count-cfg.discardFramesPerReplica),half=retained/2
        let fingerprint=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        guard retained>=cfg.minimumRetainedFramesPerReplica else {
            return .init(schema:VivoMolecularSamplingResult.schema,request:request,requestFingerprint:fingerprint,
                         retainedFramesPerReplica:retained,diagnostics:[],converged:false,
                         issues:["insufficient retained molecular frames"],interpretation:interpretation)
        }
        var diagnostics:[VivoMolecularObservableDiagnostics]=[],work=0
        for (index,observable) in request.observables.enumerated() {
            var chains:[[Double]]=[]
            for replica in request.replicas {
                let x=replica.valuesByObservable[index]
                chains.append(Array(x[cfg.discardFramesPerReplica..<(cfg.discardFramesPerReplica+half)]))
                chains.append(Array(x[(count-half)..<count]))
            }
            let all=chains.flatMap{$0},raw=try statistics(chains)
            let normalized=rankNormalize(chains)
            let rankStats=try statistics(normalized)
            let sorted=all.sorted(),median=quantile(sorted,0.5)
            let folded=chains.map{$0.map{abs($0-median)}}
            let foldedStats=try statistics(rankNormalize(folded))
            let bulk=try effectiveSamples(normalized,configuration:cfg,work:&work)
            let rawESS=try effectiveSamples(chains,configuration:cfg,work:&work)
            let low=quantile(sorted,0.05),high=quantile(sorted,0.95)
            let lowESS=try effectiveSamples(chains.map{$0.map{$0<=low ? 1.0:0.0}},configuration:cfg,work:&work)
            let highESS=try effectiveSamples(chains.map{$0.map{$0>=high ? 1.0:0.0}},configuration:cfg,work:&work)
            let tail:Double?
            if let l=lowESS.value,let h=highESS.value { tail=min(l,h) } else { tail=nil }
            let foldedRHat:Double? = folded.flatMap{$0}.allSatisfy({$0==folded[0][0]}) ? 1 : foldedStats.rHat
            let standardError=rawESS.value.map{sqrt(max(0,raw.variance)/$0)}
            let resolved=bulk.resolved && rawESS.resolved && lowESS.resolved && highESS.resolved
            let target=cfg.minimumEffectiveSamplesPerReplica*Double(request.replicas.count)
            var issues:[String]=[]
            if rankStats.rHat == nil || foldedRHat == nil { issues.append("constant or degenerate replica variance") }
            if let r=rankStats.rHat,let f=foldedRHat,max(r,f)>cfg.maximumRHat { issues.append("rank/folded split-Rhat exceeds tolerance") }
            if (bulk.value ?? 0)<target { issues.append("bulk effective sample count below target") }
            if (tail ?? 0)<target { issues.append("tail effective sample count below target") }
            if !resolved { issues.append("autocorrelation tail is unresolved within the lag budget") }
            if let se=standardError {
                if !se.isFinite || se>observable.maximumMeanStandardError { issues.append("mean standard error exceeds the observable precision") }
            } else { issues.append("mean uncertainty is not estimable") }
            diagnostics.append(.init(identifier:observable.identifier,mean:raw.mean,meanStandardError:standardError,
                rankNormalizedSplitRHat:rankStats.rHat,foldedRankNormalizedSplitRHat:foldedRHat,
                bulkEffectiveSamples:bulk.value,tailEffectiveSamples:tail,autocorrelationResolved:resolved,issues:issues))
        }
        let issues=diagnostics.flatMap{d in d.issues.map{"\(d.identifier): \($0)"}}
        return .init(schema:VivoMolecularSamplingResult.schema,request:request,requestFingerprint:fingerprint,
                     retainedFramesPerReplica:2*half,diagnostics:diagnostics,converged:issues.isEmpty,
                     issues:issues,interpretation:interpretation)
    }
    public static func validate(_ result: VivoMolecularSamplingResult) throws {
        guard result.schema==VivoMolecularSamplingResult.schema,result == (try analyze(result.request)) else {
            throw VivoArtifactValidationError.invalid("sampling result does not reconstruct from its measurements and criteria")
        }
    }

    private struct Statistics {
        let mean:Double
        let means:[Double]
        let within:Double
        let variance:Double
        let rHat:Double?
    }
    private static func statistics(_ chains:[[Double]]) throws -> Statistics {
        let n=chains[0].count,m=chains.count
        var means:[Double]=[],variances:[Double]=[]
        for x in chains {
            var mean=0.0,m2=0.0
            for (i,v) in x.enumerated() { let delta=v-mean;mean+=delta/Double(i+1);m2+=delta*(v-mean) }
            means.append(mean);variances.append(m2/Double(n-1))
        }
        let mean=means.reduce(0,+)/Double(m),within=variances.reduce(0,+)/Double(m)
        let between=Double(n)*means.reduce(0.0){$0+pow($1-mean,2)}/Double(m-1)
        let variance=Double(n-1)/Double(n)*within+between/Double(n)
        guard mean.isFinite,within.isFinite,variance.isFinite,within>=0,variance>=0 else {
            throw VivoArtifactValidationError.invalid("molecular sampling variance overflow")
        }
        let rHat:Double? = within>0 ? max(1,sqrt(variance/within)):nil
        return .init(mean:mean,means:means,within:within,variance:variance,rHat:rHat)
    }
    private static func effectiveSamples(_ chains:[[Double]],configuration cfg:VivoMolecularSamplingConfiguration,
                                         work:inout Int) throws -> (value:Double?,resolved:Bool) {
        let stats=try statistics(chains),n=chains[0].count,m=chains.count
        guard stats.within>0,stats.variance>0 else { return (nil,false) }
        let maximum=min(n-1,cfg.maximumAutocorrelationLag)
        func rho(_ lag:Int) throws -> Double {
            if lag==0 { return 1 }
            var covariance=0.0
            for chain in 0..<m {
                for i in 0..<(n-lag) {
                    covariance+=(chains[chain][i]-stats.means[chain])*(chains[chain][i+lag]-stats.means[chain])
                }
            }
            covariance/=Double(m*n)
            return 1-(stats.within-covariance)/stats.variance
        }
        var sum=0.0,previous=Double.infinity,resolved=false
        for lag in stride(from:0,to:maximum,by:2) {
            let amount=(2*n-lag-lag-1)*m
            guard amount<=cfg.maximumAutocorrelationWork-work else {
                throw VivoArtifactValidationError.invalid("molecular autocorrelation work budget exhausted")
            }
            work+=amount
            let pair=try rho(lag)+rho(lag+1)
            guard pair.isFinite else { throw VivoArtifactValidationError.invalid("nonfinite molecular autocorrelation") }
            if pair<=0 { resolved=true;break }
            previous=min(previous,pair);sum+=previous
        }
        let tau=max(1,-1+2*sum),ess=Double(m*n)/tau
        return (ess.isFinite && ess>0 ? ess:nil,resolved)
    }
    private static func quantile(_ sorted:[Double],_ probability:Double)->Double {
        let index=probability*Double(sorted.count-1),low=Int(index),high=min(sorted.count-1,low+1)
        return sorted[low]+(sorted[high]-sorted[low])*(index-Double(low))
    }
    private static func rankNormalize(_ chains:[[Double]])->[[Double]] {
        let flat=chains.flatMap{$0},n=flat.count,length=chains[0].count
        let order=flat.indices.sorted { flat[$0] == flat[$1] ? $0<$1 : flat[$0]<flat[$1] }
        var ranked=[Double](repeating:0,count:n),start=0
        while start<n {
            var end=start+1
            while end<n && flat[order[end]]==flat[order[start]] { end+=1 }
            let rank=0.5*Double(start+1+end)
            let p=(rank-0.375)/(Double(n)+0.25)
            // Bounded normal quantile inversion. No external statistics runtime.
            var low = -9.0,high = 9.0
            for _ in 0..<56 {
                let middle=0.5*(low+high),cdf=0.5*erfc(-middle/sqrt(2))
                if cdf<p { low=middle } else { high=middle }
            }
            let z=0.5*(low+high)
            for i in start..<end { ranked[order[i]]=z }
            start=end
        }
        return chains.indices.map{Array(ranked[($0*length)..<(($0+1)*length)])}
    }
}
