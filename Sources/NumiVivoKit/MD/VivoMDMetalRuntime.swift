import Foundation
import Metal
import NumiVivoShaders

public enum VivoMDRuntimeError: Error, Sendable, CustomStringConvertible {
    case unsupported([String])
    case metal(String)
    public var description: String {
        switch self {
        case .unsupported(let values): return "MD execution is unsupported: " + values.joined(separator: "; ")
        case .metal(let value): return "MD Metal runtime: \(value)"
        }
    }
}

public struct VivoMDStepCertificate: Codable, Sendable, Equatable {
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let stepIndex: UInt64
    public let timeBeforePS: Double
    public let timeAfterPS: Double
    public let committed: Bool
    public let statusFlags: UInt32
    public let firstViolationParticle: UInt32?
    public let violationCount: UInt32
}

private struct VivoMDCommandABI {
    var particleCount: UInt32; var typeCount: UInt32; var electrostatics: UInt32; var periodic: UInt32
    var dtPS: Float; var cutoffNM: Float; var coulombPrefactor: Float; var reactionFieldK: Float
    var reactionFieldC: Float; var minimumDistanceNM: Float; var reserved0: Float = 0; var reserved1: Float = 0
    var cellA: SIMD4<Float>; var cellB: SIMD4<Float>; var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>; var reciprocalB: SIMD4<Float>; var reciprocalC: SIMD4<Float>
}
private struct VivoMDStatusABI { var flags: UInt32; var firstParticle: UInt32; var violationCount: UInt32; var reserved: UInt32 }

private enum VivoMDABIValidator {
    static func validate() throws {
        let expected: [(String, Int, Int)] = [
            ("Command", MemoryLayout<VivoMDCommandABI>.stride, 144),
            ("Status", MemoryLayout<VivoMDStatusABI>.stride, 16),
            ("Bond", MemoryLayout<VivoMDPackedSystem.Bond>.stride, 16),
            ("Angle", MemoryLayout<VivoMDPackedSystem.Angle>.stride, 32),
            ("Torsion", MemoryLayout<VivoMDPackedSystem.Torsion>.stride, 32),
            ("Incidence", MemoryLayout<VivoMDPackedSystem.Incidence>.stride, 8),
            ("PairException", MemoryLayout<VivoMDPackedSystem.PairException>.stride, 32)
        ]
        for (name, actual, required) in expected where actual != required {
            throw VivoMDRuntimeError.metal("Swift/Metal ABI mismatch for \(name): \(actual) != \(required)")
        }
    }
}

private final class VivoMDMetalArena: @unchecked Sendable {
    let particleCount: Int
    let positionA, positionB, velocityA, velocityB: MTLBuffer
    let forceEnergy, kinetic, status: MTLBuffer
    let dynamics, typeIndices, pairMatrix: MTLBuffer
    let bonds, bondOffsets, bondIncidence: MTLBuffer
    let angles, angleOffsets, angleIncidence: MTLBuffer
    let torsions, torsionOffsets, torsionIncidence: MTLBuffer
    let exceptions, exceptionOffsets, exceptionPartners, exceptionIndices: MTLBuffer
    let positionReadback, velocityReadback: MTLBuffer
    private(set) var acceptedIsA = true
    var acceptedPosition: MTLBuffer { acceptedIsA ? positionA : positionB }
    var candidatePosition: MTLBuffer { acceptedIsA ? positionB : positionA }
    var acceptedVelocity: MTLBuffer { acceptedIsA ? velocityA : velocityB }
    var candidateVelocity: MTLBuffer { acceptedIsA ? velocityB : velocityA }
    func commit() { acceptedIsA.toggle() }

    static func make(device: MTLDevice, queue: MTLCommandQueue, packed: VivoMDPackedSystem,
                     initial: VivoClassicalInitialState, velocities: [VivoVector3D]) async throws -> VivoMDMetalArena {
        try VivoMDABIValidator.validate()
        let count = Int(packed.particleCount)
        guard initial.positionsNM.count == count, velocities.count == count else {
            throw VivoMDRuntimeError.metal("initial position/velocity shape mismatch")
        }
        let positions = try initial.positionsNM.map { try vector($0, "position") }
        let velocityVectors = try velocities.map { try vector($0, "velocity") }
        guard let upload = queue.makeCommandBuffer(), let blit = upload.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("initial upload command unavailable")
        }
        var staging: [MTLBuffer] = []
        func immutable<T>(_ values: [T], _ label: String) throws -> MTLBuffer {
            let length = max(values.count * MemoryLayout<T>.stride, 16)
            guard length <= device.maxBufferLength,
                  let target = device.makeBuffer(length: length, options: .storageModePrivate) else {
                throw VivoMDRuntimeError.metal("cannot allocate \(label) (\(length) bytes)")
            }
            target.label = label
            if !values.isEmpty {
                let stage: MTLBuffer? = values.withUnsafeBufferPointer { buffer in
                    guard let base = buffer.baseAddress else { return nil }
                    return device.makeBuffer(bytes: base, length: values.count * MemoryLayout<T>.stride,
                                             options: .storageModeShared)
                }
                guard let stage else { throw VivoMDRuntimeError.metal("cannot stage \(label)") }
                staging.append(stage)
                blit.copy(from: stage, sourceOffset: 0, to: target, destinationOffset: 0,
                          size: values.count * MemoryLayout<T>.stride)
            }
            return target
        }
        let positionA=try immutable(positions,"NumiVivo.MD.positionA"), positionB=try immutable(positions,"NumiVivo.MD.positionB")
        let velocityA=try immutable(velocityVectors,"NumiVivo.MD.velocityA"), velocityB=try immutable(velocityVectors,"NumiVivo.MD.velocityB")
        let dynamics=try immutable(packed.particleDynamics,"NumiVivo.MD.dynamics")
        let typeIndices=try immutable(packed.particleTypeIndices,"NumiVivo.MD.typeIndices")
        let pairMatrix=try immutable(packed.typePairC12C6,"NumiVivo.MD.pairMatrix")
        let bonds=try immutable(packed.bonds,"NumiVivo.MD.bonds"), bondOffsets=try immutable(packed.bondOffsets,"NumiVivo.MD.bondOffsets"), bondIncidence=try immutable(packed.bondIncidence,"NumiVivo.MD.bondIncidence")
        let angles=try immutable(packed.angles,"NumiVivo.MD.angles"), angleOffsets=try immutable(packed.angleOffsets,"NumiVivo.MD.angleOffsets"), angleIncidence=try immutable(packed.angleIncidence,"NumiVivo.MD.angleIncidence")
        let torsions=try immutable(packed.torsions,"NumiVivo.MD.torsions"), torsionOffsets=try immutable(packed.torsionOffsets,"NumiVivo.MD.torsionOffsets"), torsionIncidence=try immutable(packed.torsionIncidence,"NumiVivo.MD.torsionIncidence")
        let exceptions=try immutable(packed.pairExceptions,"NumiVivo.MD.exceptions"), exceptionOffsets=try immutable(packed.exceptionOffsets,"NumiVivo.MD.exceptionOffsets"), exceptionPartners=try immutable(packed.exceptionPartners,"NumiVivo.MD.exceptionPartners"), exceptionIndices=try immutable(packed.exceptionIndices,"NumiVivo.MD.exceptionIndices")
        blit.endEncoding(); try await complete(upload); _ = staging
        let vectorBytes=max(count*MemoryLayout<SIMD4<Float>>.stride,16)
        guard let forceEnergy=device.makeBuffer(length:vectorBytes,options:.storageModePrivate),
              let kinetic=device.makeBuffer(length:max(count*4,16),options:.storageModePrivate),
              let status=device.makeBuffer(length:16,options:.storageModeShared),
              let positionReadback=device.makeBuffer(length:vectorBytes,options:.storageModeShared),
              let velocityReadback=device.makeBuffer(length:vectorBytes,options:.storageModeShared) else {
            throw VivoMDRuntimeError.metal("dynamic MD buffer allocation failed")
        }
        return .init(particleCount:count,positionA:positionA,positionB:positionB,velocityA:velocityA,velocityB:velocityB,
                     forceEnergy:forceEnergy,kinetic:kinetic,status:status,dynamics:dynamics,typeIndices:typeIndices,pairMatrix:pairMatrix,
                     bonds:bonds,bondOffsets:bondOffsets,bondIncidence:bondIncidence,angles:angles,angleOffsets:angleOffsets,angleIncidence:angleIncidence,
                     torsions:torsions,torsionOffsets:torsionOffsets,torsionIncidence:torsionIncidence,exceptions:exceptions,
                     exceptionOffsets:exceptionOffsets,exceptionPartners:exceptionPartners,exceptionIndices:exceptionIndices,
                     positionReadback:positionReadback,velocityReadback:velocityReadback)
    }

    private init(particleCount:Int, positionA:MTLBuffer,positionB:MTLBuffer,velocityA:MTLBuffer,velocityB:MTLBuffer,
                 forceEnergy:MTLBuffer,kinetic:MTLBuffer,status:MTLBuffer,dynamics:MTLBuffer,typeIndices:MTLBuffer,pairMatrix:MTLBuffer,
                 bonds:MTLBuffer,bondOffsets:MTLBuffer,bondIncidence:MTLBuffer,angles:MTLBuffer,angleOffsets:MTLBuffer,angleIncidence:MTLBuffer,
                 torsions:MTLBuffer,torsionOffsets:MTLBuffer,torsionIncidence:MTLBuffer,exceptions:MTLBuffer,exceptionOffsets:MTLBuffer,
                 exceptionPartners:MTLBuffer,exceptionIndices:MTLBuffer,positionReadback:MTLBuffer,velocityReadback:MTLBuffer) {
        self.particleCount=particleCount;self.positionA=positionA;self.positionB=positionB;self.velocityA=velocityA;self.velocityB=velocityB
        self.forceEnergy=forceEnergy;self.kinetic=kinetic;self.status=status;self.dynamics=dynamics;self.typeIndices=typeIndices;self.pairMatrix=pairMatrix
        self.bonds=bonds;self.bondOffsets=bondOffsets;self.bondIncidence=bondIncidence;self.angles=angles;self.angleOffsets=angleOffsets;self.angleIncidence=angleIncidence
        self.torsions=torsions;self.torsionOffsets=torsionOffsets;self.torsionIncidence=torsionIncidence;self.exceptions=exceptions
        self.exceptionOffsets=exceptionOffsets;self.exceptionPartners=exceptionPartners;self.exceptionIndices=exceptionIndices
        self.positionReadback=positionReadback;self.velocityReadback=velocityReadback
    }
    private static func vector(_ value:VivoVector3D,_ label:String)throws->SIMD4<Float>{
        guard value.isFinite,abs(value.x)<=Double(Float.greatestFiniteMagnitude),abs(value.y)<=Double(Float.greatestFiniteMagnitude),abs(value.z)<=Double(Float.greatestFiniteMagnitude) else { throw VivoMDRuntimeError.metal("initial \(label) is outside FP32") }
        return .init(Float(value.x),Float(value.y),Float(value.z),0)
    }
    private static func complete(_ command:MTLCommandBuffer)async throws{
        try await withCheckedThrowingContinuation{(continuation:CheckedContinuation<Void,Error>)in
            command.addCompletedHandler{ value in if let error=value.error{continuation.resume(throwing:VivoMDRuntimeError.metal(String(describing:error)))}else{continuation.resume(returning:())} };command.commit()
        }
    }
}

public actor VivoMDMetalRuntime {
    public nonisolated let system: VivoClassicalSystem
    public nonisolated let configuration: VivoMDConfiguration
    public nonisolated let systemFingerprint: VivoFingerprint
    public nonisolated let configurationFingerprint: VivoFingerprint
    public nonisolated let deviceName: String
    public nonisolated let deviceRegistryID: UInt64
    private let queue:MTLCommandQueue, packed:VivoMDPackedSystem, arena:VivoMDMetalArena
    private let pipelines:[NumiVivoKernel:NumiVivoPipeline]
    private let periodicCell:VivoPeriodicCell?
    private var acceptedStep:UInt64=0,acceptedTimePS:Double=0,inFlight=false

    public static func make(system:VivoClassicalSystem,initialState:VivoClassicalInitialState,configuration:VivoMDConfiguration,
                            initialVelocitiesNMPerPS:[VivoVector3D]?=nil,device requestedDevice:MTLDevice?=nil)async throws->VivoMDMetalRuntime{
        let report=try VivoMDCapabilityAnalyzer.analyze(system:system,initialState:initialState,configuration:configuration)
        var blockers=report.blockers
        if configuration.electrostatics == .pme { blockers.append("PME reciprocal-space execution is not installed yet") }
        if !system.constraints.isEmpty { blockers.append("constraint projection is not installed yet") }
        if configuration.ensemble != .nve || configuration.thermostat != .none || configuration.barostat != .none { blockers.append("current runtime stage executes NVE only") }
        if let cell=initialState.periodicCell { try validateMinimumImageCutoff(configuration.cutoffNM,cell:cell) }
        guard blockers.isEmpty else { throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted()) }
        let packed=try VivoMDSystemPacker.pack(system),device=try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory,let queue=device.makeCommandQueue() else { throw VivoMDRuntimeError.metal("Apple-silicon unified memory and command queue are required") }
        let catalog=try NumiVivoPipelineCatalog(device:device)
        let required:[NumiVivoKernel]=[.mdClearForce,.mdClearStatus,.mdBonded,.mdNonbondedDirect,.mdHalfKick,.mdDrift,.mdKinetic,.mdValidate]
        var pipelines:[NumiVivoKernel:NumiVivoPipeline]=[:];for kernel in required{pipelines[kernel]=try await catalog.pipeline(kernel)}
        let velocities=initialVelocitiesNMPerPS ?? [VivoVector3D](repeating:.zero,count:system.particles.count)
        let arena=try await VivoMDMetalArena.make(device:device,queue:queue,packed:packed,initial:initialState,velocities:velocities)
        queue.label="NumiVivo.MD.Queue"
        return try .init(system:system,configuration:configuration,packed:packed,queue:queue,arena:arena,
                         pipelines:pipelines,periodicCell:initialState.periodicCell,device:device)
    }
    private init(system:VivoClassicalSystem,configuration:VivoMDConfiguration,packed:VivoMDPackedSystem,queue:MTLCommandQueue,
                 arena:VivoMDMetalArena,pipelines:[NumiVivoKernel:NumiVivoPipeline],periodicCell:VivoPeriodicCell?,device:MTLDevice)throws{
        self.system=system;self.configuration=configuration;self.packed=packed;self.queue=queue;self.arena=arena;self.pipelines=pipelines;self.periodicCell=periodicCell
        systemFingerprint=packed.systemFingerprint;configurationFingerprint=try configuration.fingerprint();deviceName=device.name;deviceRegistryID=device.registryID
    }

    public func step()async throws->VivoMDStepCertificate{
        guard !inFlight else{throw VivoMDRuntimeError.metal("MD operation already in flight")};guard acceptedStep<UInt64.max else{throw VivoMDRuntimeError.metal("step index overflow")}
        inFlight=true;defer{inFlight=false};let before=acceptedTimePS,abi=try commandABI()
        guard let command=queue.makeCommandBuffer(),let blit=command.makeBlitCommandEncoder()else{throw VivoMDRuntimeError.metal("step command unavailable")}
        command.label="NumiVivo.MD.step.\(acceptedStep)";let bytes=arena.particleCount*MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from:arena.acceptedPosition,sourceOffset:0,to:arena.candidatePosition,destinationOffset:0,size:bytes)
        blit.copy(from:arena.acceptedVelocity,sourceOffset:0,to:arena.candidateVelocity,destinationOffset:0,size:bytes);blit.endEncoding()
        try encodeStatusClear(command);try encodeForces(command,position:arena.candidatePosition,abi:abi)
        try encode(.mdHalfKick,command:command,buffers:[arena.candidateVelocity,arena.forceEnergy,arena.dynamics],abi:abi)
        try encode(.mdDrift,command:command,buffers:[arena.candidatePosition,arena.candidateVelocity,arena.dynamics],abi:abi)
        try encodeForces(command,position:arena.candidatePosition,abi:abi)
        try encode(.mdHalfKick,command:command,buffers:[arena.candidateVelocity,arena.forceEnergy,arena.dynamics],abi:abi)
        try encode(.mdValidate,command:command,buffers:[arena.candidatePosition,arena.candidateVelocity,arena.forceEnergy,arena.status],abi:abi)
        try await complete(command);let status=arena.status.contents().assumingMemoryBound(to:VivoMDStatusABI.self).pointee
        let first=status.firstParticle==UInt32.max ? nil : status.firstParticle
        guard status.flags==0 else{return .init(systemFingerprint:systemFingerprint,configurationFingerprint:configurationFingerprint,deviceName:deviceName,deviceRegistryID:deviceRegistryID,stepIndex:acceptedStep,timeBeforePS:before,timeAfterPS:before,committed:false,statusFlags:status.flags,firstViolationParticle:first,violationCount:status.violationCount)}
        arena.commit();acceptedStep+=1;acceptedTimePS+=configuration.timeStepPS
        return .init(systemFingerprint:systemFingerprint,configurationFingerprint:configurationFingerprint,deviceName:deviceName,deviceRegistryID:deviceRegistryID,stepIndex:acceptedStep-1,timeBeforePS:before,timeAfterPS:acceptedTimePS,committed:true,statusFlags:0,firstViolationParticle:nil,violationCount:0)
    }

    public func snapshot()async throws->VivoMDStateSnapshot{
        guard !inFlight else{throw VivoMDRuntimeError.metal("cannot snapshot during MD operation")}
        guard let command=queue.makeCommandBuffer(),let blit=command.makeBlitCommandEncoder()else{throw VivoMDRuntimeError.metal("snapshot command unavailable")}
        let bytes=arena.particleCount*MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from:arena.acceptedPosition,sourceOffset:0,to:arena.positionReadback,destinationOffset:0,size:bytes);blit.copy(from:arena.acceptedVelocity,sourceOffset:0,to:arena.velocityReadback,destinationOffset:0,size:bytes);blit.endEncoding();try await complete(command)
        let pp=arena.positionReadback.contents().assumingMemoryBound(to:SIMD4<Float>.self),vv=arena.velocityReadback.contents().assumingMemoryBound(to:SIMD4<Float>.self)
        var p:[VivoVector3D]=[],v:[VivoVector3D]=[];p.reserveCapacity(arena.particleCount);v.reserveCapacity(arena.particleCount)
        for i in 0..<arena.particleCount{p.append(.init(Double(pp[i].x),Double(pp[i].y),Double(pp[i].z)));v.append(.init(Double(vv[i].x),Double(vv[i].y),Double(vv[i].z)))}
        let result=VivoMDStateSnapshot(systemFingerprint:systemFingerprint,configurationFingerprint:configurationFingerprint,stepIndex:acceptedStep,timePS:acceptedTimePS,positionsNM:p,velocitiesNMPerPS:v,periodicCell:periodicCell);try result.validate(particleCount:arena.particleCount);return result
    }

    private func encodeForces(_ command:MTLCommandBuffer,position:MTLBuffer,abi:VivoMDCommandABI)throws{
        try encode(.mdClearForce,command:command,buffers:[arena.forceEnergy],abi:abi)
        try encode(.mdBonded,command:command,buffers:[position,arena.forceEnergy,arena.bonds,arena.bondOffsets,arena.bondIncidence,arena.angles,arena.angleOffsets,arena.angleIncidence,arena.torsions,arena.torsionOffsets,arena.torsionIncidence,arena.status],abi:abi)
        try encode(.mdNonbondedDirect,command:command,buffers:[position,arena.forceEnergy,arena.dynamics,arena.typeIndices,arena.pairMatrix,arena.exceptions,arena.exceptionOffsets,arena.exceptionPartners,arena.exceptionIndices,arena.status],abi:abi)
    }
    private func encodeStatusClear(_ command:MTLCommandBuffer)throws{
        guard let pipeline=pipelines[.mdClearStatus],let encoder=command.makeComputeCommandEncoder()else{throw VivoMDRuntimeError.metal("status-clear encoder unavailable")};encoder.setComputePipelineState(pipeline.state);encoder.setBuffer(arena.status,offset:0,index:0);encoder.dispatchThreads(.init(width:1,height:1,depth:1),threadsPerThreadgroup:.init(width:1,height:1,depth:1));encoder.endEncoding()
    }
    private func encode(_ kernel:NumiVivoKernel,command:MTLCommandBuffer,buffers:[MTLBuffer],abi:VivoMDCommandABI)throws{
        guard let pipeline=pipelines[kernel],let encoder=command.makeComputeCommandEncoder()else{throw VivoMDRuntimeError.metal("encoder unavailable for \(kernel.rawValue)")};encoder.label=kernel.rawValue;encoder.setComputePipelineState(pipeline.state)
        for(i,b)in buffers.enumerated(){encoder.setBuffer(b,offset:0,index:i)};var value=abi;encoder.setBytes(&value,length:MemoryLayout<VivoMDCommandABI>.stride,index:buffers.count);encoder.dispatchThreads(pipeline.gridSize(for:arena.particleCount),threadsPerThreadgroup:pipeline.threadgroupSize(for:arena.particleCount));encoder.endEncoding()
    }
    private func commandABI()throws->VivoMDCommandABI{
        let mode:UInt32=configuration.electrostatics == .cutoff ? 0 : 1;let eps=configuration.relativeDielectric,cut=configuration.cutoffNM
        var krf:Float=0,crf:Float=0;if configuration.electrostatics == .reactionField{let er=configuration.reactionFieldDielectric;krf=Float((er-eps)/(2*er+eps)/pow(cut,3));crf=Float(3*er/(2*er+eps)/cut)}
        let reciprocal=try periodicCell.map(reciprocalCell);func f(_ x:VivoVector3D?)->SIMD4<Float>{guard let x else{return .zero};return .init(Float(x.x),Float(x.y),Float(x.z),0)}
        return .init(particleCount:packed.particleCount,typeCount:packed.typeCount,electrostatics:mode,periodic:periodicCell==nil ? 0:1,dtPS:Float(configuration.timeStepPS),cutoffNM:Float(cut),coulombPrefactor:Float(138.935456/eps),reactionFieldK:krf,reactionFieldC:crf,minimumDistanceNM:1e-5,cellA:f(periodicCell?.a),cellB:f(periodicCell?.b),cellC:f(periodicCell?.c),reciprocalA:f(reciprocal?.0),reciprocalB:f(reciprocal?.1),reciprocalC:f(reciprocal?.2))
    }
    private func reciprocalCell(_ c:VivoPeriodicCell)throws->(VivoVector3D,VivoVector3D,VivoVector3D){let d=c.a.dot(c.b.cross(c.c));guard d.isFinite,abs(d)>1e-18 else{throw VivoMDRuntimeError.metal("singular periodic cell")};return(c.b.cross(c.c)/d,c.c.cross(c.a)/d,c.a.cross(c.b)/d)}
    private static func validateMinimumImageCutoff(_ cutoff:Double,cell:VivoPeriodicCell)throws{let v=cell.a.dot(cell.b.cross(cell.c));let h=[abs(v)/cell.b.cross(cell.c).norm,abs(v)/cell.c.cross(cell.a).norm,abs(v)/cell.a.cross(cell.b).norm];guard h.allSatisfy({$0.isFinite&&$0>0}),cutoff<0.5*(h.min() ?? 0) else{throw VivoMDRuntimeError.unsupported(["cutoff must be less than half the shortest periodic face height for direct minimum-image execution"])} }
    private func complete(_ command:MTLCommandBuffer)async throws{try await withCheckedThrowingContinuation{(continuation:CheckedContinuation<Void,Error>)in command.addCompletedHandler{v in if let e=v.error{continuation.resume(throwing:VivoMDRuntimeError.metal(String(describing:e)))}else{continuation.resume(returning:())}};command.commit()}}
}
