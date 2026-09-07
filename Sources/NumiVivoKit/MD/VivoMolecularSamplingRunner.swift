import Foundation

public struct VivoMolecularSamplingRunRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/molecular-sampling-run/v1"
    public var schema:String
    public var structure:VivoMolecularStructure
    public var system:VivoClassicalSystem
    public var initialStates:[VivoClassicalInitialState]
    public var replicaSeeds:[UInt64]
    public var md:VivoMDConfiguration
    public var contextIdentifier:String
    public var observables:[VivoMolecularObservable]
    public var convergence:VivoMolecularSamplingConfiguration
    public var minimization:VivoMDMinimizationConfiguration?
    public var equilibrationSteps:UInt64
    public var stepsPerBlock:UInt64
    public var sampleEvery:UInt64
    public var maximumBlocks:Int
    public var requiredConsecutivePasses:Int
    public var trajectoryChunkBytes:Int
    public init(structure:VivoMolecularStructure,system:VivoClassicalSystem,
                initialStates:[VivoClassicalInitialState],replicaSeeds:[UInt64],md:VivoMDConfiguration,
                contextIdentifier:String,observables:[VivoMolecularObservable],
                convergence:VivoMolecularSamplingConfiguration = .init(),
                minimization:VivoMDMinimizationConfiguration? = nil,
                equilibrationSteps:UInt64,stepsPerBlock:UInt64,sampleEvery:UInt64,
                maximumBlocks:Int,requiredConsecutivePasses:Int=2,trajectoryChunkBytes:Int=8*1024*1024) {
        schema=Self.schema;self.structure=structure;self.system=system;self.initialStates=initialStates
        self.replicaSeeds=replicaSeeds;self.md=md;self.contextIdentifier=contextIdentifier;self.observables=observables
        self.convergence=convergence;self.minimization=minimization;self.equilibrationSteps=equilibrationSteps
        self.stepsPerBlock=stepsPerBlock;self.sampleEvery=sampleEvery;self.maximumBlocks=maximumBlocks
        self.requiredConsecutivePasses=requiredConsecutivePasses;self.trajectoryChunkBytes=trajectoryChunkBytes
    }
    public func validate() throws {
        _ = try VivoStructureValidator.validate(structure)
        try VivoClassicalSystemValidator.validate(system,atomCount:UInt32(structure.atoms.count))
        try md.validate();try convergence.validate();try minimization?.validate()
        guard schema==Self.schema,system.structureFingerprint == (try VivoStructureCodec.fingerprint(structure)),
              !contextIdentifier.isEmpty,!observables.isEmpty,observables.count<=64,
              initialStates.count==replicaSeeds.count,replicaSeeds.count>=convergence.minimumReplicas,
              replicaSeeds.count<=64,Set(replicaSeeds).count==replicaSeeds.count,
              Set(observables.map(\.identifier)).count==observables.count,
              md.thermostat == .langevinMiddle,md.frictionPerPS.map({$0>0}) == true,
              md.targetTemperatureK != nil,stepsPerBlock>0,sampleEvery>0,
              stepsPerBlock%sampleEvery==0,stepsPerBlock/sampleEvery>=2,
              (1...10000).contains(maximumBlocks),(2...100).contains(requiredConsecutivePasses),
              requiredConsecutivePasses<=maximumBlocks,
              trajectoryChunkBytes>VivoMDTrajectoryChunkCodec.headerBytes,
              trajectoryChunkBytes<=VivoMDTrajectoryChunkCodec.maximumChunkBytes else {
            throw VivoArtifactValidationError.invalid("adaptive molecular sampling request, mapping, replicas or block schedule")
        }
        let blocks=stepsPerBlock.multipliedReportingOverflow(by:UInt64(maximumBlocks))
        guard !blocks.overflow,blocks.partialValue<=UInt64.max-equilibrationSteps else {
            throw VivoArtifactValidationError.invalid("molecular sampling step budget overflow")
        }
        let frames=(blocks.partialValue/sampleEvery).multipliedReportingOverflow(by:UInt64(replicaSeeds.count))
        let scalars=frames.partialValue.multipliedReportingOverflow(by:UInt64(observables.count))
        guard !frames.overflow,!scalars.overflow,scalars.partialValue<=UInt64(convergence.maximumTotalSamples) else {
            throw VivoArtifactValidationError.invalid("sampling schedule exceeds the scalar-observation capacity")
        }
        let id=try system.fingerprint(),map=try VivoMolecularSamplingRunner.atomMap(structure:structure,system:system)
        for initial in initialStates {
            try initial.validate(particleCount:system.particles.count)
            guard initial.systemFingerprint==id else { throw VivoArtifactValidationError.incompatible("replica initial state differs from the prepared classical system") }
            for observable in observables {
                try observable.validate()
                switch observable.kind {
                case .potentialEnergyKJPerMol,.temperatureK:break
                default:_ = try observable.measure(positionsNM:initial.positionsNM,periodicCell:initial.periodicCell,particleByAtom:map)
                }
            }
            let capabilities=try VivoMDCapabilityAnalyzer.analyze(system:system,initialState:initial,configuration:md)
            guard capabilities.executable else { throw VivoMDRuntimeError.unsupported(capabilities.blockers) }
        }
    }
}
public enum VivoMolecularSamplingRunStatus:String,Codable,Sendable {
    case converged,budgetExhausted,rejected,cancelled
}
public struct VivoMolecularSamplingRunReceipt:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/molecular-sampling-run-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let status:VivoMolecularSamplingRunStatus
    public let completedBlocks:Int
    public let consecutivePasses:Int
    public let checkpoint:VivoFingerprint
    public let trajectoryManifests:[VivoFingerprint]
    public let sampling:VivoMolecularSamplingResult?
    public let diagnostic:String?
    public let interpretation:String
}
struct VivoMolecularReplicaCursor:Codable,Sendable,Equatable {
    var series:VivoMolecularReplicaSeries
    var mdCheckpoint:VivoFingerprint?
    var trajectory:VivoFingerprint?
    var minimization:VivoFingerprint?
}
struct VivoMolecularSamplingCursor:Codable,Sendable,Equatable {
    static let schema="numivivo.org/molecular-sampling-checkpoint/v1"
    var schema:String
    var requestFingerprint:VivoFingerprint
    var completedBlocks:Int
    var consecutivePasses:Int
    var replicas:[VivoMolecularReplicaCursor]
    var diagnostics:[VivoFingerprint]
}

/// Sampling orchestration, not a second MD engine. A single replica's existing
/// Metal arena is resident at a time. FP32-exact checkpoints permit sequential
/// GPU reuse, while coordinates stream through the existing trajectory archive.
/// A complete cross-replica block is the durable publication boundary.
public enum VivoMolecularSamplingRunner {
    public static func checkpointReferenceName(requestFingerprint:VivoFingerprint) -> String {
        "molecular-sampling-"+requestFingerprint.hex+"-checkpoint"
    }
    private static func publish(_ cursor:VivoMolecularSamplingCursor,store:VivoArtifactStore) async throws -> VivoFingerprint {
        let artifact=try await store.put(data:VivoCanonicalJSON.encode(cursor),kind:"molecular-sampling-checkpoint",mediaType:"application/json")
        _ = try await store.setReference(checkpointReferenceName(requestFingerprint:cursor.requestFingerprint),to:artifact)
        return artifact.fingerprint
    }
    static func atomMap(structure:VivoMolecularStructure,system:VivoClassicalSystem) throws -> [Int] {
        var map=[Int](repeating:-1,count:structure.atoms.count)
        for particle in system.particles where particle.role == .atom {
            guard let atom=particle.atomIndex,Int(atom)<map.count,map[Int(atom)] == -1 else {
                throw VivoArtifactValidationError.invalid("sampling requires a unique structure-atom/physical-particle map")
            }
            map[Int(atom)]=Int(particle.index)
        }
        guard map.allSatisfy({$0>=0}) else { throw VivoArtifactValidationError.unresolved("sampling structure contains unmapped physical atoms") }
        return map
    }
    private static func put<T:Encodable>(_ value:T,kind:String,store:VivoArtifactStore) async throws -> VivoFingerprint {
        try await store.put(data:VivoCanonicalJSON.encode(value),kind:kind,mediaType:"application/json").fingerprint
    }
    private static func read<T:Decodable>(_ type:T.Type,id:VivoFingerprint,kind:String,store:VivoArtifactStore) async throws -> T {
        let descriptor=try await store.descriptor(for:id)
        guard descriptor.kind==kind else { throw VivoArtifactValidationError.incompatible("sampling artifact kind mismatch") }
        return try VivoCanonicalJSON.decode(type,from:await store.data(for:id,verify:true))
    }

    public static func run(_ request:VivoMolecularSamplingRunRequest,store:VivoArtifactStore,
                           resumeFrom:VivoFingerprint? = nil,
                           resumeReadLimits:VivoMolecularSamplingReadLimits = .init()) async throws -> VivoMolecularSamplingRunReceipt {
        try request.validate()
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        let systemID=try request.system.fingerprint(),structureID=try VivoStructureCodec.fingerprint(request.structure)
        let map=try atomMap(structure:request.structure,system:request.system)
        var cursor:VivoMolecularSamplingCursor
        var latest:VivoMolecularSamplingResult?
        var resumeValidations:[VivoMolecularSamplingReplicaRestart]=[]
        if let resumeFrom {
            let reader=try await VivoMolecularSamplingArchiveReader.open(store:store,checkpoint:resumeFrom,limits:resumeReadLimits)
            guard reader.cursor.requestFingerprint==requestID else {
                throw VivoArtifactValidationError.incompatible("sampling resume request identity differs")
            }
            resumeValidations=try await reader.validateForResume(request:request)
            cursor=reader.cursor;latest=reader.latest
        } else {
            _ = try await put(request,kind:"molecular-sampling-run",store:store)
            var replicas:[VivoMolecularReplicaCursor]=[]
            for (i,initial) in request.initialStates.enumerated() {
                let initialID=try await put(initial,kind:"classical-initial-state",store:store)
                var cfg=request.md;cfg.randomSeed=request.replicaSeeds[i]
                let series=VivoMolecularReplicaSeries(identifier:"replica-\(i)",sourceFingerprint:initialID,configuration:cfg,
                    steps:[],timesPS:[],valuesByObservable:Array(repeating:[],count:request.observables.count))
                replicas.append(.init(series:series,mdCheckpoint:nil,trajectory:nil,minimization:nil))
            }
            cursor = .init(schema:VivoMolecularSamplingCursor.schema,requestFingerprint:requestID,
                           completedBlocks:0,consecutivePasses:0,replicas:replicas,diagnostics:[])
        }
        var durable=try await publish(cursor,store:store)
        func receipt(_ status:VivoMolecularSamplingRunStatus,_ diagnostic:String? = nil)->VivoMolecularSamplingRunReceipt {
            .init(schema:VivoMolecularSamplingRunReceipt.schema,requestFingerprint:requestID,status:status,
                  completedBlocks:cursor.completedBlocks,consecutivePasses:cursor.consecutivePasses,checkpoint:durable,
                  trajectoryManifests:cursor.replicas.compactMap(\.trajectory),sampling:latest,diagnostic:diagnostic,
                  interpretation:VivoMolecularSampling.interpretation)
        }
        do {
            while cursor.completedBlocks<request.maximumBlocks && cursor.consecutivePasses<request.requiredConsecutivePasses {
                try Task.checkCancellation()
                var next=cursor
                for i in next.replicas.indices {
                    // The helper returns before the next arena is created. It
                    // never thermalizes or minimizes a resumed production state.
                    next.replicas[i]=try await advanceReplica(next.replicas[i],index:i,request:request,map:map,store:store,
                        validatedReplica:resumeValidations.isEmpty ? nil:resumeValidations[i])
                }
                resumeValidations=[]
                let analysis=VivoMolecularSamplingRequest(structureFingerprint:structureID,systemFingerprint:systemID,
                    contextIdentifier:request.contextIdentifier,observables:request.observables,
                    replicas:next.replicas.map(\.series),configuration:request.convergence)
                let result=try VivoMolecularSampling.analyze(analysis)
                next.completedBlocks+=1;next.consecutivePasses=result.converged ? next.consecutivePasses+1:0
                next.diagnostics.append(try await put(result,kind:"molecular-sampling-result",store:store))
                let checkpoint=try await publish(next,store:store)
                cursor=next;durable=checkpoint;latest=result
            }
            return receipt(cursor.consecutivePasses>=request.requiredConsecutivePasses ? .converged:.budgetExhausted)
        } catch is CancellationError {
            return receipt(.cancelled,"cancelled; resume from the last complete replica block")
        } catch {
            return receipt(.rejected,String(describing:error))
        }
    }
    private static func advanceReplica(_ previous:VivoMolecularReplicaCursor,index:Int,
                                       request:VivoMolecularSamplingRunRequest,map:[Int],store:VivoArtifactStore,
                                       validatedReplica:VivoMolecularSamplingReplicaRestart? = nil) async throws -> VivoMolecularReplicaCursor {
        var cursor=previous
        let cfg=cursor.series.configuration,runtime:VivoMDMetalRuntime
        if let validatedReplica {
            guard validatedReplica.checkpointFingerprint==cursor.mdCheckpoint,
                  validatedReplica.trajectoryValidation.manifestFingerprint==cursor.trajectory else {
                throw VivoArtifactValidationError.incompatible("sampling continuation differs from validated replica state")
            }
            runtime=try await .restore(system:request.system,configuration:cfg,checkpoint:validatedReplica.checkpoint)
        } else if let id=cursor.mdCheckpoint {
            let checkpoint=try await read(VivoMDCheckpoint.self,id:id,kind:"md-checkpoint",store:store)
            runtime=try await .restore(system:request.system,configuration:cfg,checkpoint:checkpoint)
        } else {
            runtime=try await .make(system:request.system,initialState:request.initialStates[index],configuration:cfg)
            if let settings=request.minimization {
                let minimized=try await runtime.minimize(settings)
                cursor.minimization=try await put(minimized,kind:"md-minimization",store:store)
                guard minimized.converged else { throw VivoArtifactValidationError.invalid("replica minimization did not converge") }
            }
            _ = try await runtime.thermalize(temperatureK:cfg.targetTemperatureK!,seed:cfg.randomSeed ^ 0x53414d504c455448)
            for _ in 0..<request.equilibrationSteps {
                try Task.checkCancellation()
                let step=try await runtime.step()
                guard step.committed else { throw VivoMDRuntimeError.candidateRejected(step.statusFlags) }
            }
        }
        let writer:VivoMDTrajectoryArchiveWriter
        if let validatedReplica {
            writer=try .resume(validated:validatedReplica.trajectoryValidation,targetChunkBytes:request.trajectoryChunkBytes)
        } else if let manifest=cursor.trajectory {
            writer=try await .resume(store:store,manifest:manifest,targetChunkBytes:request.trajectoryChunkBytes)
        } else {
            writer=try .init(store:store,systemFingerprint:runtime.systemFingerprint,
                configurationFingerprint:runtime.configurationFingerprint,particleCount:UInt32(request.system.particles.count),
                includeVelocities:false,targetChunkBytes:request.trajectoryChunkBytes)
        }
        for index in 1...request.stepsPerBlock {
            try Task.checkCancellation()
            let step=try await runtime.step()
            guard step.committed else { throw VivoMDRuntimeError.candidateRejected(step.statusFlags) }
            if index%request.sampleEvery==0 {
                let sample=try await runtime.sample(includeObservables:true)
                guard let obs=sample.observables else { throw VivoArtifactValidationError.invalid("MD sample omitted observables") }
                try await writer.append(sample.state)
                cursor.series.steps.append(sample.state.stepIndex);cursor.series.timesPS.append(sample.state.timePS)
                for j in request.observables.indices {
                    cursor.series.valuesByObservable[j].append(try request.observables[j].measure(
                        positionsNM:sample.state.positionsNM,periodicCell:sample.state.periodicCell,particleByAtom:map,
                        potentialEnergyKJPerMol:obs.potentialEnergyKJPerMol,temperatureK:obs.temperatureK))
                }
            }
        }
        cursor.mdCheckpoint=try await put(runtime.checkpoint(),kind:"md-checkpoint",store:store)
        cursor.trajectory=try await writer.snapshot().fingerprint
        cursor.series.sourceFingerprint=cursor.trajectory!
        return cursor
    }
}
