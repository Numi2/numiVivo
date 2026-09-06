import Foundation

public struct VivoQMMMFreeEnergyExecutionConfiguration: Codable, Sendable, Equatable {
    public var productionStepsPerCheckpointBlock:UInt64
    public var maximumTotalScalarSamples:Int
    public init(productionStepsPerCheckpointBlock:UInt64=5000,maximumTotalScalarSamples:Int=2_000_000) {
        self.productionStepsPerCheckpointBlock=productionStepsPerCheckpointBlock;self.maximumTotalScalarSamples=maximumTotalScalarSamples
    }
    public func validate(sampleEvery:UInt64) throws {
        guard productionStepsPerCheckpointBlock>0,sampleEvery>0,productionStepsPerCheckpointBlock%sampleEvery==0,
              maximumTotalScalarSamples>0,maximumTotalScalarSamples<=20_000_000 else {
            throw VivoChemistryError.invalid("free-energy checkpoint block or scalar capacity")
        }
    }
}

public struct VivoQMMMFreeEnergyExecutionRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-execution/v1"
    public var schema:String
    public var sampling:VivoQMMMFreeEnergyRunRequest
    public var systemFingerprint:VivoFingerprint
    public var baseProviderFingerprint:VivoFingerprint
    public var execution:VivoQMMMFreeEnergyExecutionConfiguration
    public init(sampling:VivoQMMMFreeEnergyRunRequest,systemFingerprint:VivoFingerprint,
                baseProviderFingerprint:VivoFingerprint,execution:VivoQMMMFreeEnergyExecutionConfiguration = .init()) {
        schema=Self.schema;self.sampling=sampling;self.systemFingerprint=systemFingerprint
        self.baseProviderFingerprint=baseProviderFingerprint;self.execution=execution
    }
}

public enum VivoQMMMFreeEnergyExecutionStatus:String,Codable,Sendable {
    case converged
    case analysisNotConverged
    case rejected
    case cancelled
}

public struct VivoQMMMFreeEnergyExecutionReceipt:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/qmmm-free-energy-execution-result/v1"
    public let schema:String
    public let requestFingerprint:VivoFingerprint
    public let status:VivoQMMMFreeEnergyExecutionStatus
    public let checkpoint:VivoFingerprint
    public let completedWindows:Int
    public let currentWindowProductionSteps:UInt64
    public let result:VivoQMMMActivationFreeEnergyResult?
    public let diagnostic:String?
}

private struct VivoQMMMFreeEnergyExecutionCursor:Codable,Sendable,Equatable {
    static let schema="numivivo.org/qmmm-free-energy-execution-checkpoint/v1"
    var schema:String
    var requestFingerprint:VivoFingerprint
    var completedWindows:Int
    var currentWindowProductionSteps:UInt64
    var currentMDCheckpoint:VivoFingerprint?
    var traces:[VivoQMMMUmbrellaTrace]
    var currentCoordinates:[Double]
    var currentEnergies:[Double]
}

/// Durable production orchestration over the existing Metal MD + composed BO
/// provider. A checkpoint is published only after a complete production block.
/// The checkpoint contains scalar CV/energy traces plus an exact accepted-state
/// MD checkpoint; it never serializes a closure or mutable SCF history.
public enum VivoQMMMFreeEnergyArchiveRunner {
    public static func checkpointReferenceName(requestFingerprint:VivoFingerprint)->String {
        "qmmm-free-energy/\(requestFingerprint.hex)/checkpoint"
    }
    private static func put<T:Encodable>(_ value:T,kind:String,store:VivoArtifactStore) async throws -> VivoFingerprint {
        try await store.put(data:VivoCanonicalJSON.encode(value),kind:kind,mediaType:"application/json").fingerprint
    }
    private static func read<T:Decodable>(_ type:T.Type,id:VivoFingerprint,kind:String,store:VivoArtifactStore) async throws -> T {
        let descriptor=try await store.descriptor(for:id)
        guard descriptor.kind==kind else { throw VivoChemistryError.invalid("free-energy artifact kind mismatch") }
        return try VivoCanonicalJSON.decode(type,from:await store.data(for:id,verify:true))
    }

    public static func run(_ request:VivoQMMMFreeEnergyExecutionRequest,system:VivoClassicalSystem,
                           baseProvider:VivoMDCandidateForceProvider,store:VivoArtifactStore,
                           resumeFrom:VivoFingerprint?=nil) async throws -> VivoQMMMFreeEnergyExecutionReceipt {
        let sampling=request.sampling
        try sampling.coordinate.validate();try sampling.analysis.validate();try sampling.dynamics.validate()
        try request.execution.validate(sampleEvery:sampling.sampleEvery)
        let systemID=try system.fingerprint()
        guard request.schema==VivoQMMMFreeEnergyExecutionRequest.schema,request.systemFingerprint==systemID,
              request.baseProviderFingerprint==baseProvider.fingerprint,baseProvider.retainedSystemFingerprint==systemID,
              sampling.windows.count>=2,sampling.windows.count==sampling.initialStates.count,
              sampling.windows.count==sampling.randomSeeds.count,Set(sampling.randomSeeds).count==sampling.randomSeeds.count,
              sampling.productionSteps>0,sampling.sampleEvery>0,sampling.productionSteps%sampling.sampleEvery==0,
              sampling.dynamics.ensemble == .nvt,sampling.dynamics.thermostat == .langevinMiddle,
              sampling.dynamics.targetTemperatureK != nil,sampling.dynamics.frictionPerPS.map({$0>0}) == true else {
            throw VivoChemistryError.invalid("artifact-backed free-energy execution identity or NVT schedule")
        }
        let sampleCount=sampling.productionSteps/sampling.sampleEvery
        let aggregate=sampleCount.multipliedReportingOverflow(by:UInt64(sampling.windows.count))
        guard !aggregate.overflow,aggregate.partialValue<=UInt64(request.execution.maximumTotalScalarSamples) else {
            throw VivoChemistryError.resourceLimit("free-energy scalar trace capacity")
        }
        for i in sampling.windows.indices {
            try sampling.windows[i].validate();try sampling.initialStates[i].validate(particleCount:system.particles.count)
            guard sampling.initialStates[i].systemFingerprint==systemID else { throw VivoChemistryError.invalid("free-energy initial state/system mismatch") }
        }
        let requestID=try await put(request,kind:"qmmm-free-energy-execution-request",store:store)
        var cursor:VivoQMMMFreeEnergyExecutionCursor
        if let resumeFrom {
            cursor=try await read(VivoQMMMFreeEnergyExecutionCursor.self,id:resumeFrom,kind:"qmmm-free-energy-execution-checkpoint",store:store)
            guard cursor.schema==VivoQMMMFreeEnergyExecutionCursor.schema,cursor.requestFingerprint==requestID,
                  cursor.completedWindows>=0,cursor.completedWindows<=sampling.windows.count,
                  cursor.currentWindowProductionSteps<=sampling.productionSteps,
                  cursor.currentWindowProductionSteps%sampling.sampleEvery==0,
                  cursor.traces.count==cursor.completedWindows,
                  cursor.currentCoordinates.count==Int(cursor.currentWindowProductionSteps/sampling.sampleEvery),
                  cursor.currentEnergies.count==cursor.currentCoordinates.count,
                  cursor.currentCoordinates.allSatisfy(\.isFinite),cursor.currentEnergies.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("free-energy resume cursor differs from the execution request")
            }
            for i in cursor.traces.indices {
                guard cursor.traces[i].window==sampling.windows[i],cursor.traces[i].randomSeed==sampling.randomSeeds[i],
                      cursor.traces[i].coordinateNM.count==Int(sampleCount),cursor.traces[i].potentialEnergyKJPerMol.count==Int(sampleCount) else {
                    throw VivoChemistryError.invalid("completed free-energy window trace differs from its request")
                }
            }
            if cursor.completedWindows<sampling.windows.count,cursor.currentWindowProductionSteps>0 {
                guard let id=cursor.currentMDCheckpoint else { throw VivoChemistryError.invalid("free-energy resume lacks exact MD checkpoint") }
                let checkpoint=try await read(VivoMDCheckpoint.self,id:id,kind:"md-checkpoint",store:store)
                var dynamics=sampling.dynamics;dynamics.randomSeed=sampling.randomSeeds[cursor.completedWindows]
                let biased=try VivoQMMMUmbrellaBias.provider(base:baseProvider,system:system,coordinate:sampling.coordinate,window:sampling.windows[cursor.completedWindows])
                let executionID=try VivoMDCandidateForceProvider.executionFingerprint(configuration:dynamics,provider:biased)
                guard checkpoint.systemFingerprint==systemID,checkpoint.configurationFingerprint==executionID,
                      checkpoint.acceptedStep==sampling.equilibrationSteps+cursor.currentWindowProductionSteps else {
                    throw VivoChemistryError.invalid("free-energy MD checkpoint clock or Hamiltonian identity mismatch")
                }
            }
        } else {
            cursor = .init(schema:VivoQMMMFreeEnergyExecutionCursor.schema,requestFingerprint:requestID,completedWindows:0,
                           currentWindowProductionSteps:0,currentMDCheckpoint:nil,traces:[],currentCoordinates:[],currentEnergies:[])
        }
        var durable=try await put(cursor,kind:"qmmm-free-energy-execution-checkpoint",store:store)
        func receipt(_ status:VivoQMMMFreeEnergyExecutionStatus,_ result:VivoQMMMActivationFreeEnergyResult?=nil,_ diagnostic:String?=nil)->VivoQMMMFreeEnergyExecutionReceipt {
            .init(schema:VivoQMMMFreeEnergyExecutionReceipt.schema,requestFingerprint:requestID,status:status,checkpoint:durable,
                  completedWindows:cursor.completedWindows,currentWindowProductionSteps:cursor.currentWindowProductionSteps,result:result,diagnostic:diagnostic)
        }
        do {
            while cursor.completedWindows<sampling.windows.count {
                try Task.checkCancellation()
                let index=cursor.completedWindows,window=sampling.windows[index]
                var dynamics=sampling.dynamics;dynamics.randomSeed=sampling.randomSeeds[index]
                let biased=try VivoQMMMUmbrellaBias.provider(base:baseProvider,system:system,coordinate:sampling.coordinate,window:window)
                let runtime:VivoMDMetalRuntime
                if let checkpointID=cursor.currentMDCheckpoint {
                    let checkpoint=try await read(VivoMDCheckpoint.self,id:checkpointID,kind:"md-checkpoint",store:store)
                    runtime=try await .restore(system:system,configuration:dynamics,checkpoint:checkpoint,forceProvider:biased)
                } else {
                    runtime=try await .make(system:system,initialState:sampling.initialStates[index],configuration:dynamics,forceProvider:biased)
                    _ = try await runtime.thermalize(temperatureK:dynamics.targetTemperatureK!,seed:sampling.randomSeeds[index]^0x554d4252454c4c41)
                    for _ in 0..<sampling.equilibrationSteps {
                        try Task.checkCancellation();guard try await runtime.step().committed else { throw VivoChemistryError.convergence("umbrella equilibration candidate rejected") }
                    }
                    let exact=try await runtime.checkpoint()
                    cursor.currentMDCheckpoint=try await put(exact,kind:"md-checkpoint",store:store)
                }
                let coordinate=try VivoQMMMResolvedCoordinate(source:sampling.coordinate,system:system)
                let remaining=sampling.productionSteps-cursor.currentWindowProductionSteps
                let block=min(remaining,request.execution.productionStepsPerCheckpointBlock)
                guard block%sampling.sampleEvery==0 else {
                    throw VivoChemistryError.invalid("final free-energy production block is not aligned to sampleEvery")
                }
                for step in 1...block {
                    try Task.checkCancellation();guard try await runtime.step().committed else { throw VivoChemistryError.convergence("umbrella production candidate rejected") }
                    if step%sampling.sampleEvery==0 {
                        let sample=try await runtime.sample(includeObservables:true)
                        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:sample.state.positionsNM,periodicCell:sample.state.periodicCell)
                        cursor.currentCoordinates.append(try coordinate.evaluate(geometry).valueNM)
                        guard let energy=sample.observables?.potentialEnergyKJPerMol else { throw VivoChemistryError.invalid("umbrella sample omitted potential energy") }
                        cursor.currentEnergies.append(energy)
                    }
                }
                cursor.currentWindowProductionSteps+=block
                let exact=try await runtime.checkpoint()
                cursor.currentMDCheckpoint=try await put(exact,kind:"md-checkpoint",store:store)
                if cursor.currentWindowProductionSteps==sampling.productionSteps {
                    cursor.traces.append(.init(window:window,randomSeed:sampling.randomSeeds[index],coordinateNM:cursor.currentCoordinates,
                                               potentialEnergyKJPerMol:cursor.currentEnergies))
                    cursor.completedWindows+=1;cursor.currentWindowProductionSteps=0;cursor.currentMDCheckpoint=nil
                    cursor.currentCoordinates=[];cursor.currentEnergies=[]
                }
                durable=try await put(cursor,kind:"qmmm-free-energy-execution-checkpoint",store:store)
            }
            let result=try VivoQMMMFreeEnergy.analyze(coordinate:sampling.coordinate,temperatureK:sampling.dynamics.targetTemperatureK!,
                                                       traces:cursor.traces,configuration:sampling.analysis)
            _ = try await put(result,kind:"qmmm-activation-free-energy",store:store)
            return receipt(result.converged ? .converged:.analysisNotConverged,result)
        } catch is CancellationError {
            return receipt(.cancelled,nil,"cancelled; resume from the last complete production block")
        } catch {
            return receipt(.rejected,nil,String(describing:error))
        }
    }
}
