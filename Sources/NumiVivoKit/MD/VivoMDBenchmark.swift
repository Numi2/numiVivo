import Foundation

/// A comparison against supplied, separately generated reference observations.
/// Passing this gate does not authenticate a reference engine or establish an ensemble.
public enum VivoBenchmarkOutcome: String, Codable, Sendable { case passed, failed, inconclusive, unsupported }
public enum VivoMDBenchmarkPreparation: String, Codable, Sendable { case preserve, projectConstraints }

public struct VivoMDBenchmarkLimits: Codable, Sendable, Equatable {
    public var energyAbsolutePerParticleKJPerMol: Double = 0.002
    public var forceNormalizedRMS: Double = 0.001
    public var forceNormalizedMaximum: Double = 0.01
    public var forceNormalizationFloor: Double = 1
    public init() {}
    public func validate() throws {
        guard [energyAbsolutePerParticleKJPerMol, forceNormalizedRMS, forceNormalizedMaximum,
               forceNormalizationFloor].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoChemistryError.invalid("benchmark limits must be finite and positive")
        }
    }
}

public struct VivoMDBenchmarkReference: Codable, Sendable, Equatable {
    public var identifier: String
    public var geometry: VivoMDCandidateGeometry
    public var energyKJPerMol: Double
    public var forcesKJPerMolNM: [VivoVector3D]
    public init(identifier: String, geometry: VivoMDCandidateGeometry, energyKJPerMol: Double,
                forcesKJPerMolNM: [VivoVector3D]) {
        self.identifier=identifier; self.geometry=geometry; self.energyKJPerMol=energyKJPerMol
        self.forcesKJPerMolNM=forcesKJPerMolNM
    }
}

public struct VivoMDBenchmarkRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-benchmark-request/v1"
    public var schema: String = Self.schema
    public var identifier: String
    public var system: VivoClassicalSystem
    public var configuration: VivoMDConfiguration
    public var references: [VivoMDBenchmarkReference]
    /// Descriptive reference identity. The campaign runner must additionally
    /// verify source input, serialized reference-system and request bytes.
    public var referenceProvenance: [String: String]
    public var limits: VivoMDBenchmarkLimits
    public var dynamicsSteps: Int
    public var dynamicsPreparation: VivoMDBenchmarkPreparation?
    public init(identifier: String, system: VivoClassicalSystem, configuration: VivoMDConfiguration,
                references: [VivoMDBenchmarkReference], referenceProvenance: [String: String],
                limits: VivoMDBenchmarkLimits = .init(), dynamicsSteps: Int = 0,
                dynamicsPreparation: VivoMDBenchmarkPreparation? = nil) {
        self.identifier=identifier; self.system=system; self.configuration=configuration
        self.references=references; self.referenceProvenance=referenceProvenance
        self.limits=limits; self.dynamicsSteps=dynamicsSteps
        self.dynamicsPreparation=dynamicsPreparation
    }
    public func validate() throws {
        try limits.validate(); try configuration.validate(); try VivoClassicalSystemValidator.validate(system)
        guard schema == Self.schema, !identifier.isEmpty, !referenceProvenance.isEmpty,
              system.particles.count <= 1_000_000, (1...32).contains(references.count),
              (0...100_000).contains(dynamicsSteps), Set(references.map(\.identifier)).count == references.count else {
            throw VivoChemistryError.invalid("benchmark schema, identities or bounded work")
        }
        for ref in references {
            guard !ref.identifier.isEmpty, ref.energyKJPerMol.isFinite,
                  ref.geometry.particlePositionsNM.count == system.particles.count,
                  ref.geometry.particlePositionsNM.allSatisfy(\.isFinite), ref.geometry.periodicCell?.isValid != false,
                  ref.forcesKJPerMolNM.count == system.particles.count, ref.forcesKJPerMolNM.allSatisfy(\.isFinite),
                  ref.geometry.particlePositionsNM.allSatisfy({ v in [v.x,v.y,v.z].allSatisfy { Double(Float($0)) == $0 } }) else {
                throw VivoChemistryError.invalid("benchmark reference shape or exact FP32 coordinates")
            }
        }
    }
}

public struct VivoMDBenchmarkComparison: Codable, Sendable, Equatable {
    public let identifier: String
    public let energyErrorPerParticleKJPerMol: Double
    public let forceNormalizedRMS: Double
    public let forceNormalizedMaximum: Double
    public let geometryUnchanged: Bool
    public let outcome: VivoBenchmarkOutcome
    public static func compare(_ actual: VivoMDHamiltonianEvaluation, reference: VivoMDBenchmarkReference,
                               limits: VivoMDBenchmarkLimits) throws -> Self {
        try limits.validate()
        let n=reference.forcesKJPerMolNM.count
        guard n>0, actual.physicalParticleForcesKJPerMolNM.count==n,
              actual.energyKJPerMol.isFinite, reference.energyKJPerMol.isFinite,
              reference.forcesKJPerMolNM.allSatisfy(\.isFinite), actual.physicalParticleForcesKJPerMolNM.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("benchmark comparison shape or nonfinite values")
        }
        let scale=max(limits.forceNormalizationFloor, sqrt(reference.forcesKJPerMolNM.reduce(0) { $0+$1.squaredNorm }/Double(3*n)))
        let errors=zip(actual.physicalParticleForcesKJPerMolNM,reference.forcesKJPerMolNM).map { $0-$1 }
        let rms=sqrt(errors.reduce(0) { $0+$1.squaredNorm }/Double(3*n))/scale
        let maximum=errors.reduce(0.0) { max($0,abs($1.x),abs($1.y),abs($1.z)) }/scale
        let energy=abs(actual.energyKJPerMol-reference.energyKJPerMol)/Double(n)
        guard [scale,rms,maximum,energy].allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("benchmark reduction overflow") }
        let same=actual.evaluatedGeometry==reference.geometry
        return .init(identifier: reference.identifier, energyErrorPerParticleKJPerMol: energy,
            forceNormalizedRMS: rms, forceNormalizedMaximum: maximum, geometryUnchanged: same,
            outcome: same && energy<=limits.energyAbsolutePerParticleKJPerMol && rms<=limits.forceNormalizedRMS && maximum<=limits.forceNormalizedMaximum ? .passed:.failed)
    }
}

public struct VivoMDBenchmarkDynamics: Codable, Sendable, Equatable {
    public let requestedSteps: Int
    public let committedSteps: Int
    public let wallSeconds: Double
    public let rejected: VivoMDStepCertificate?
    public let start: VivoMDObservables
    public let preparedCheckpoint: VivoMDCheckpoint
    public let end: VivoMDObservables
    public let finalCheckpoint: VivoMDCheckpoint
    /// This short trajectory is execution evidence, not an equilibrium gate.
    public let ensembleOutcome: VivoBenchmarkOutcome = .inconclusive
}

public struct VivoMDBenchmarkReport: Codable, Sendable {
    public let schema: String = "numivivo.org/md-benchmark-report/v1"
    public let requestFingerprint: VivoFingerprint
    public let numericalContract: String
    public let identifier: String
    public let capability: VivoMDCapabilityReport
    public let deviceName: String?
    public let comparisons: [VivoMDBenchmarkComparison]
    public let evaluations: [VivoMDHamiltonianEvaluation]
    public let dynamics: VivoMDBenchmarkDynamics?
    public let executionError: String?
    public let outcome: VivoBenchmarkOutcome
}

public enum VivoMDBenchmark {
    public static func run(_ request: VivoMDBenchmarkRequest) async throws -> VivoMDBenchmarkReport {
        try request.validate(); try Task.checkCancellation()
        let first=request.references[0], system=request.system, config=request.configuration
        let initial=VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),
            positionsNM: first.geometry.particlePositionsNM, periodicCell:first.geometry.periodicCell)
        let capability=try VivoMDCapabilityAnalyzer.analyze(system:system,initialState:initial,configuration:config)
        let identity=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        guard capability.executable else {
            return .init(requestFingerprint:identity,numericalContract:VivoMDExecutionIdentity.current,
                identifier:request.identifier,capability:capability,deviceName:nil,comparisons:[],evaluations:[],dynamics:nil,executionError:nil,outcome:.unsupported)
        }
        let runtime=try await VivoMDMetalRuntime.make(system:system,initialState:initial,configuration:config)
        let before=try await runtime.checkpoint()
        var comparisons:[VivoMDBenchmarkComparison]=[], evaluations:[VivoMDHamiltonianEvaluation]=[]
        for ref in request.references {
            let value=try await runtime.evaluateHamiltonian(at:ref.geometry)
            evaluations.append(value)
            comparisons.append(try .compare(value,reference:ref,limits:request.limits))
        }
        let after=try await runtime.checkpoint()
        guard before==after else { throw VivoChemistryError.invalid("benchmark probes changed accepted state") }
        var dynamics:VivoMDBenchmarkDynamics?
        var executionError:String?
        if request.dynamicsSteps>0 {
          do {
            if request.dynamicsPreparation == .projectConstraints { _ = try await runtime.projectConstraints() }
            let prepared=try await runtime.thermalize(temperatureK:config.targetTemperatureK ?? 300,seed:config.randomSeed)
            let start=try await runtime.observables(), clock=ContinuousClock(), began=clock.now
            var committed=0,rejected:VivoMDStepCertificate?
            for _ in 0..<request.dynamicsSteps {
                try Task.checkCancellation()
                let step=try await runtime.step()
                if !step.committed { rejected=step;break };committed+=1
            }
            let duration=began.duration(to:clock.now).components
            let seconds=Double(duration.seconds)+Double(duration.attoseconds)/1e18
            let end=try await runtime.observables(), checkpoint=try await runtime.checkpoint()
            dynamics = .init(requestedSteps:request.dynamicsSteps,committedSteps:committed,wallSeconds:seconds,
                rejected:rejected,start:start,preparedCheckpoint:prepared,end:end,finalCheckpoint:checkpoint)
          } catch is CancellationError { throw CancellationError() }
          catch { executionError=String(describing:error) }
        }
        return .init(requestFingerprint:identity,numericalContract:VivoMDExecutionIdentity.current,
            identifier:request.identifier,capability:capability,deviceName:runtime.deviceName,
            comparisons:comparisons,evaluations:evaluations,dynamics:dynamics,executionError:executionError,
            outcome:comparisons.allSatisfy({$0.outcome == .passed}) && dynamics?.rejected == nil && executionError == nil ? .passed:.failed)
    }
}
