import Foundation
import Metal

final class VivoMDGPUArena: @unchecked Sendable {
    let particleCount:Int
    let neighborCapacity:UInt32
    let positionA,positionB,velocityA,velocityB:MTLBuffer
    let positionScratch,velocityScratch:MTLBuffer
    let forceEnergy,kinetic,status:MTLBuffer
    let dynamics,typeIndices,pairMatrix:MTLBuffer
    let bonds,bondOffsets,bondIncidence:MTLBuffer
    let angles,angleOffsets,angleIncidence:MTLBuffer
    let torsions,torsionOffsets,torsionIncidence:MTLBuffer
    let constraints,constraintOffsets,constraintIncidence:MTLBuffer
    let exceptions,exceptionOffsets,exceptionPartners,exceptionIndices:MTLBuffer
    let neighborCounts,neighborIndices,neighborReferencePosition:MTLBuffer
    let positionReadback,velocityReadback,forceEnergyReadback,kineticReadback:MTLBuffer
    private(set)var acceptedIsA=true
    var acceptedPosition:MTLBuffer{acceptedIsA ? positionA:positionB}
    var candidatePosition:MTLBuffer{acceptedIsA ? positionB:positionA}
    var acceptedVelocity:MTLBuffer{acceptedIsA ? velocityA:velocityB}
    var candidateVelocity:MTLBuffer{acceptedIsA ? velocityB:velocityA}
    func commit(){acceptedIsA.toggle()}

    static func make(device:MTLDevice,queue:MTLCommandQueue,packed:VivoMDPackedSystem,
                     initial:VivoClassicalInitialState,velocities:[VivoVector3D],
                     configuration:VivoMDConfiguration)async throws->VivoMDGPUArena{
        try VivoMDMetalABI.validate();let count=Int(packed.particleCount)
        guard initial.positionsNM.count==count,velocities.count==count else{throw VivoMDRuntimeError.metal("initial position/velocity shape mismatch")}
        let positions=try initial.positionsNM.map{try vector($0,label:"position")},velocityVectors=try velocities.map{try vector($0,label:"velocity")}
        let constraintAdjacency=try makeConstraintIncidence(packed.constraints,particleCount:count)
        let maximumPossible=max(count-1,1),requested=Int(configuration.resolvedMaximumNeighborsPerParticle)
        let neighborCapacity=UInt32(max(1,min(requested,maximumPossible)))
        let neighborElements=UInt64(count)*UInt64(neighborCapacity)
        guard neighborElements<=UInt64(Int.max),neighborElements<=UInt64(device.maxBufferLength/MemoryLayout<UInt32>.stride) else{throw VivoMDRuntimeError.metal("Verlet neighbor-list capacity exceeds Metal buffer limits")}

        guard let upload=queue.makeCommandBuffer(),let blit=upload.makeBlitCommandEncoder()else{throw VivoMDRuntimeError.metal("initial upload command unavailable")}
        var staging:[MTLBuffer]=[]
        func immutable<T>(_ values:[T],label:String)throws->MTLBuffer{
            let length=max(values.count*MemoryLayout<T>.stride,16)
            guard length<=device.maxBufferLength,let target=device.makeBuffer(length:length,options:.storageModePrivate)else{throw VivoMDRuntimeError.metal("cannot allocate \(label) (\(length) bytes)")}
            target.label=label
            if !values.isEmpty{let stage:MTLBuffer?=values.withUnsafeBufferPointer{source in guard let base=source.baseAddress else{return nil};return device.makeBuffer(bytes:base,length:values.count*MemoryLayout<T>.stride,options:.storageModeShared)};guard let stage else{throw VivoMDRuntimeError.metal("cannot stage \(label)")};staging.append(stage);blit.copy(from:stage,sourceOffset:0,to:target,destinationOffset:0,size:values.count*MemoryLayout<T>.stride)}
            return target
        }
        let positionA=try immutable(positions,label:"NumiVivo.MD.positionA"),positionB=try immutable(positions,label:"NumiVivo.MD.positionB")
        let velocityA=try immutable(velocityVectors,label:"NumiVivo.MD.velocityA"),velocityB=try immutable(velocityVectors,label:"NumiVivo.MD.velocityB")
        let dynamics=try immutable(packed.particleDynamics,label:"NumiVivo.MD.dynamics"),typeIndices=try immutable(packed.particleTypeIndices,label:"NumiVivo.MD.typeIndices"),pairMatrix=try immutable(packed.typePairC12C6,label:"NumiVivo.MD.pairMatrix")
        let bonds=try immutable(packed.bonds,label:"NumiVivo.MD.bonds"),bondOffsets=try immutable(packed.bondOffsets,label:"NumiVivo.MD.bondOffsets"),bondIncidence=try immutable(packed.bondIncidence,label:"NumiVivo.MD.bondIncidence")
        let angles=try immutable(packed.angles,label:"NumiVivo.MD.angles"),angleOffsets=try immutable(packed.angleOffsets,label:"NumiVivo.MD.angleOffsets"),angleIncidence=try immutable(packed.angleIncidence,label:"NumiVivo.MD.angleIncidence")
        let torsions=try immutable(packed.torsions,label:"NumiVivo.MD.torsions"),torsionOffsets=try immutable(packed.torsionOffsets,label:"NumiVivo.MD.torsionOffsets"),torsionIncidence=try immutable(packed.torsionIncidence,label:"NumiVivo.MD.torsionIncidence")
        let constraints=try immutable(packed.constraints,label:"NumiVivo.MD.constraints"),constraintOffsets=try immutable(constraintAdjacency.offsets,label:"NumiVivo.MD.constraintOffsets"),constraintIncidence=try immutable(constraintAdjacency.entries,label:"NumiVivo.MD.constraintIncidence")
        let exceptions=try immutable(packed.pairExceptions,label:"NumiVivo.MD.exceptions"),exceptionOffsets=try immutable(packed.exceptionOffsets,label:"NumiVivo.MD.exceptionOffsets"),exceptionPartners=try immutable(packed.exceptionPartners,label:"NumiVivo.MD.exceptionPartners"),exceptionIndices=try immutable(packed.exceptionIndices,label:"NumiVivo.MD.exceptionIndices")
        blit.endEncoding();try await complete(upload);_=staging

        let vectorBytes=max(count*MemoryLayout<SIMD4<Float>>.stride,16),scalarBytes=max(count*MemoryLayout<Float>.stride,16)
        func buffer(_ length:Int,_ storage:MTLResourceOptions,_ label:String)throws->MTLBuffer{guard length<=device.maxBufferLength,let value=device.makeBuffer(length:length,options:storage)else{throw VivoMDRuntimeError.metal("cannot allocate \(label)")};value.label=label;return value}
        let neighborBytes=max(Int(neighborElements)*MemoryLayout<UInt32>.stride,16)
        return try .init(particleCount:count,neighborCapacity:neighborCapacity,positionA:positionA,positionB:positionB,velocityA:velocityA,velocityB:velocityB,
                         positionScratch:buffer(vectorBytes,.storageModePrivate,"NumiVivo.MD.positionScratch"),velocityScratch:buffer(vectorBytes,.storageModePrivate,"NumiVivo.MD.velocityScratch"),
                         forceEnergy:buffer(vectorBytes,.storageModePrivate,"NumiVivo.MD.forceEnergy"),kinetic:buffer(scalarBytes,.storageModePrivate,"NumiVivo.MD.kinetic"),status:buffer(16,.storageModeShared,"NumiVivo.MD.status"),
                         dynamics:dynamics,typeIndices:typeIndices,pairMatrix:pairMatrix,bonds:bonds,bondOffsets:bondOffsets,bondIncidence:bondIncidence,
                         angles:angles,angleOffsets:angleOffsets,angleIncidence:angleIncidence,torsions:torsions,torsionOffsets:torsionOffsets,torsionIncidence:torsionIncidence,
                         constraints:constraints,constraintOffsets:constraintOffsets,constraintIncidence:constraintIncidence,exceptions:exceptions,exceptionOffsets:exceptionOffsets,exceptionPartners:exceptionPartners,exceptionIndices:exceptionIndices,
                         neighborCounts:buffer(max(count*4,16),.storageModePrivate,"NumiVivo.MD.neighborCounts"),neighborIndices:buffer(neighborBytes,.storageModePrivate,"NumiVivo.MD.neighborIndices"),neighborReferencePosition:buffer(vectorBytes,.storageModePrivate,"NumiVivo.MD.neighborReferencePosition"),
                         positionReadback:buffer(vectorBytes,.storageModeShared,"NumiVivo.MD.positionReadback"),velocityReadback:buffer(vectorBytes,.storageModeShared,"NumiVivo.MD.velocityReadback"),forceEnergyReadback:buffer(vectorBytes,.storageModeShared,"NumiVivo.MD.forceReadback"),kineticReadback:buffer(scalarBytes,.storageModeShared,"NumiVivo.MD.kineticReadback"))
    }

    private init(particleCount:Int,neighborCapacity:UInt32,positionA:MTLBuffer,positionB:MTLBuffer,velocityA:MTLBuffer,velocityB:MTLBuffer,positionScratch:MTLBuffer,velocityScratch:MTLBuffer,forceEnergy:MTLBuffer,kinetic:MTLBuffer,status:MTLBuffer,dynamics:MTLBuffer,typeIndices:MTLBuffer,pairMatrix:MTLBuffer,bonds:MTLBuffer,bondOffsets:MTLBuffer,bondIncidence:MTLBuffer,angles:MTLBuffer,angleOffsets:MTLBuffer,angleIncidence:MTLBuffer,torsions:MTLBuffer,torsionOffsets:MTLBuffer,torsionIncidence:MTLBuffer,constraints:MTLBuffer,constraintOffsets:MTLBuffer,constraintIncidence:MTLBuffer,exceptions:MTLBuffer,exceptionOffsets:MTLBuffer,exceptionPartners:MTLBuffer,exceptionIndices:MTLBuffer,neighborCounts:MTLBuffer,neighborIndices:MTLBuffer,neighborReferencePosition:MTLBuffer,positionReadback:MTLBuffer,velocityReadback:MTLBuffer,forceEnergyReadback:MTLBuffer,kineticReadback:MTLBuffer){
        self.particleCount=particleCount;self.neighborCapacity=neighborCapacity;self.positionA=positionA;self.positionB=positionB;self.velocityA=velocityA;self.velocityB=velocityB;self.positionScratch=positionScratch;self.velocityScratch=velocityScratch;self.forceEnergy=forceEnergy;self.kinetic=kinetic;self.status=status;self.dynamics=dynamics;self.typeIndices=typeIndices;self.pairMatrix=pairMatrix;self.bonds=bonds;self.bondOffsets=bondOffsets;self.bondIncidence=bondIncidence;self.angles=angles;self.angleOffsets=angleOffsets;self.angleIncidence=angleIncidence;self.torsions=torsions;self.torsionOffsets=torsionOffsets;self.torsionIncidence=torsionIncidence;self.constraints=constraints;self.constraintOffsets=constraintOffsets;self.constraintIncidence=constraintIncidence;self.exceptions=exceptions;self.exceptionOffsets=exceptionOffsets;self.exceptionPartners=exceptionPartners;self.exceptionIndices=exceptionIndices;self.neighborCounts=neighborCounts;self.neighborIndices=neighborIndices;self.neighborReferencePosition=neighborReferencePosition;self.positionReadback=positionReadback;self.velocityReadback=velocityReadback;self.forceEnergyReadback=forceEnergyReadback;self.kineticReadback=kineticReadback
    }
    private static func makeConstraintIncidence(_ constraints:[VivoMDPackedSystem.Constraint],particleCount:Int)throws->(offsets:[UInt32],entries:[VivoMDPackedSystem.Incidence]){var buckets=[[VivoMDPackedSystem.Incidence]](repeating:[],count:particleCount);for(index,constraint)in constraints.enumerated(){guard let term=UInt32(exactly:index),constraint.atoms.x<UInt32(particleCount),constraint.atoms.y<UInt32(particleCount)else{throw VivoMDRuntimeError.metal("constraint incidence exceeds execution ABI")};buckets[Int(constraint.atoms.x)].append(.init(termIndex:term,localIndex:0));buckets[Int(constraint.atoms.y)].append(.init(termIndex:term,localIndex:1))};var offsets:[UInt32]=[0],entries:[VivoMDPackedSystem.Incidence]=[];for bucket in buckets{entries.append(contentsOf:bucket);guard let next=UInt32(exactly:entries.count)else{throw VivoMDRuntimeError.metal("constraint incidence exceeds UInt32")};offsets.append(next)};return(offsets,entries)}
    private static func vector(_ value:VivoVector3D,label:String)throws->SIMD4<Float>{guard value.isFinite,abs(value.x)<=Double(Float.greatestFiniteMagnitude),abs(value.y)<=Double(Float.greatestFiniteMagnitude),abs(value.z)<=Double(Float.greatestFiniteMagnitude)else{throw VivoMDRuntimeError.metal("initial \(label) is outside FP32")};return .init(Float(value.x),Float(value.y),Float(value.z),0)}
    private static func complete(_ command:MTLCommandBuffer)async throws{try await withCheckedThrowingContinuation{(continuation:CheckedContinuation<Void,Error>)in command.addCompletedHandler{value in if let error=value.error{continuation.resume(throwing:VivoMDRuntimeError.metal(String(describing:error)))}else{continuation.resume(returning:())}};command.commit()}}
}
