import Foundation
@preconcurrency import Metal

/// Low words follow the arena's existing transactional A/B ownership. Five
/// buffers cost 80 bytes per particle; no second state or commit owner is added.
final class VivoMDPositionCorrections: @unchecked Sendable {
    let a,b,scratch,reference,readback:MTLBuffer
    init(device:MTLDevice,positions:[VivoVector3D],exactCorrections:[VivoVector3D]?=nil) throws {
        guard !positions.isEmpty,positions.count<=device.maxBufferLength/16 else {
            throw VivoMDRuntimeError.metal("compensated coordinate allocation exceeds buffer limit")
        }
        func allocate(_ name:String)throws->MTLBuffer {
            guard let buffer=device.makeBuffer(length:positions.count*16,options:.storageModeShared) else {
                throw VivoMDRuntimeError.metal("compensated coordinate allocation: " + name)
            }
            buffer.label="NumiVivo.MD.positionCorrection."+name
            buffer.contents().initializeMemory(as:UInt8.self,repeating:0,count:buffer.length)
            return buffer
        }
        a=try allocate("A");b=try allocate("B");scratch=try allocate("scratch")
        reference=try allocate("reference");readback=try allocate("readback")
        for buffer in [a,b] {
            let pointer=buffer.contents().assumingMemoryBound(to:SIMD4<Float>.self)
            for (i,p) in positions.enumerated() {
                let low=exactCorrections?[i] ?? .init(p.x-Double(Float(p.x)),p.y-Double(Float(p.y)),p.z-Double(Float(p.z)))
                pointer[i] = .init(Float(low.x),Float(low.y),Float(low.z),0)
            }
        }
    }
    func forPosition(_ position:MTLBuffer,arena:VivoMDGPUArena)throws->MTLBuffer {
        if position === arena.positionA { return a }
        if position === arena.positionB { return b }
        if position === arena.positionScratch { return scratch }
        throw VivoMDRuntimeError.metal("unowned compensated position buffer")
    }
    func accepted(_ arena:VivoMDGPUArena)->MTLBuffer { arena.acceptedIsA ? a:b }
    func candidate(_ arena:VivoMDGPUArena)->MTLBuffer { arena.acceptedIsA ? b:a }
}
