import Foundation
@preconcurrency import Metal
import NumiVivoShaders

private struct VivoMultipoleMeshSource {
    var positionCharge: SIMD4<Float>
    var dipole: SIMD4<Float>
    var q0: SIMD4<Float>
    var q1: SIMD4<Float>
}
/// The existing PME FFT schedule is shared. This resource adds derivative-
/// consistent sixth-order multipole spread/influence/gather, not a second FFT.
/// Shared buffers are serialized with a lock; calls return after GPU completion.
/// CPU-side real/self/exclusion and electronic reductions retain FP64 authority.
final class VivoMultipolePMEResource: @unchecked Sendable {
    let mesh: VivoMultipoleMeshConfiguration
    let capacity: Int
    let queue: MTLCommandQueue
    let gridA: MTLBuffer,gridB: MTLBuffer,sourceBuffer: MTLBuffer,resultBuffer: MTLBuffer,modeBuffer: MTLBuffer
    let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let lock = NSLock()
    init(mesh: VivoMultipoleMeshConfiguration,capacity: Int,queue: MTLCommandQueue,gridA: MTLBuffer,gridB: MTLBuffer,
         source: MTLBuffer,result: MTLBuffer,modes: MTLBuffer,pipelines: [NumiVivoKernel: NumiVivoPipeline]) {
        self.mesh = mesh;self.capacity = capacity;self.queue = queue;self.gridA = gridA;self.gridB = gridB
        sourceBuffer = source;resultBuffer = result;modeBuffer = modes;self.pipelines = pipelines
    }
    static func make(mesh: VivoMultipoleMeshConfiguration,capacity: Int,device: MTLDevice,
                     catalog: NumiVivoPipelineCatalog,budget: VivoChemistryBudget) async throws -> Self {
        try mesh.validate();try budget.validate()
        guard capacity >= 0,capacity <= Int(UInt32.max),MemoryLayout<VivoMultipoleMeshSource>.stride == 64,
              let queue = device.makeCommandQueue() else { throw VivoChemistryError.resourceLimit("multipolar PME capacity, ABI or queue") }
        let g = mesh.gridDimensions,count = g[0]*g[1]*g[2]
        _ = try budget.elements([count,8]);_ = try budget.elements([max(1,capacity),16],simultaneousArrays: 2)
        let bytes = count*48+max(1,capacity)*128
        guard bytes <= budget.maximumBytes else { throw VivoChemistryError.resourceLimit("multipolar PME aggregate workspace") }
        func buffer(_ bytes: Int,_ label: String,_ options: MTLResourceOptions) throws -> MTLBuffer {
            guard bytes <= device.maxBufferLength,let b = device.makeBuffer(length: max(bytes,16),options: options) else {
                throw VivoChemistryError.resourceLimit("multipolar PME buffer allocation")
            }
            b.label = label;return b
        }
        let a = try buffer(count*8,"NumiVivo.PME.multipole.gridA",.storageModePrivate)
        let b = try buffer(count*8,"NumiVivo.PME.multipole.gridB",.storageModePrivate)
        let source = try buffer(max(1,capacity)*64,"NumiVivo.PME.multipole.sources",.storageModeShared)
        let result = try buffer(max(1,capacity)*64,"NumiVivo.PME.multipole.derivatives",.storageModeShared)
        let modes = try buffer(count*32,"NumiVivo.PME.multipole.strain",.storageModeShared)
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in [NumiVivoKernel.mdPMEClearGrid,.mdPMEBitReverse,.mdPMEFFTStage,.mdPMEScaleInverse,
                       .mdPMEMultipoleSpread,.mdPMEMultipoleInfluence,.mdPMEMultipoleGather] {
            pipelines[kernel] = try await catalog.pipeline(kernel)
        }
        return .init(mesh: mesh,capacity: capacity,queue: queue,gridA: a,gridB: b,source: source,result: result,modes: modes,pipelines: pipelines)
    }
    func evaluate(_ sources: [VivoCartesianMultipole],cell: VivoPeriodicCell,
                  configuration cfg: VivoPeriodicElectrostaticConfiguration,budget: VivoChemistryBudget) throws -> VivoReciprocalElectrostaticResult {
        lock.lock();defer { lock.unlock() }
        try cfg.validate();try budget.validate()
        guard cell.isValid,sources.count <= capacity else { throw VivoChemistryError.invalid("multipolar PME sources or cell") }
        let dimensions = mesh.gridDimensions,count = dimensions[0]*dimensions[1]*dimensions[2]
        guard zip(cfg.reciprocalHalfWidths,dimensions).allSatisfy({ $0.0 < $0.1/2 }),
              count*48+max(1,capacity)*128 <= budget.maximumBytes else {
            throw VivoChemistryError.resourceLimit("multipolar PME Fourier widths or aggregate workspace")
        }
        let lattice = [cell.a,cell.b,cell.c].map { $0/VivoAtomicUnits.bohrInNM }
        let determinant = lattice[0].dot(lattice[1].cross(lattice[2])),volume = abs(determinant)
        let reciprocal = [lattice[1].cross(lattice[2])/determinant,lattice[2].cross(lattice[0])/determinant,lattice[0].cross(lattice[1])/determinant]
        let gridVectors = (0..<3).map { reciprocal[$0]*Double(dimensions[$0]) }
        func xyz(_ p: VivoVector3D) -> [Double] { [p.x,p.y,p.z] }
        let transform = gridVectors.map(xyz)
        func f4(_ p: VivoVector3D) -> SIMD4<Float> { .init(Float(p.x),Float(p.y),Float(p.z),0) }
        var command = VivoPMEMetalCommand(particleCount: UInt32(sources.count),gridX: UInt32(dimensions[0]),gridY: UInt32(dimensions[1]),
            gridZ: UInt32(dimensions[2]),gridPointCount: UInt32(count),axis: 0,stage: 0,inverse: 0,
            betaPerNM: Float(cfg.alphaPerBohr),volumeNM3: Float(volume),coulombPrefactor: 1,inverseGridCount: 1/Float(count),
            reciprocalA: f4(reciprocal[0]),reciprocalB: f4(reciprocal[1]),reciprocalC: f4(reciprocal[2]))
        guard command.volumeNM3.isFinite,command.volumeNM3 > 0,command.betaPerNM.isFinite,command.betaPerNM > 0,
              [command.reciprocalA,command.reciprocalB,command.reciprocalC].allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite }) else {
            throw VivoChemistryError.invalid("multipolar PME cell outside FP32 range")
        }
        let input = sourceBuffer.contents().assumingMemoryBound(to: VivoMultipoleMeshSource.self)
        for (i,source) in sources.enumerated() {
            try source.validate()
            let fractions = reciprocal.map { let f = $0.dot(source.positionBohr);return f-floor(f) }
            let q = source.secondMomentsEBohr2,matrix = [[q[0],q[1],q[2]],[q[1],q[3],q[4]],[q[2],q[4],q[5]]]
            var transformed = [[Double]](repeating: [0,0,0],count: 3)
            for a in 0..<3 { for b in 0..<3 { for c in 0..<3 { for d in 0..<3 {
                transformed[a][b] += transform[a][c]*matrix[c][d]*transform[b][d]
            } } } }
            let u = (0..<3).map { fractions[$0]*Double(dimensions[$0]) },mu = gridVectors.map { $0.dot(source.dipoleEBohr) }
            let values = u+[source.chargeE]+mu+transformed.flatMap { $0 }
            guard values.allSatisfy({ $0.isFinite && Float($0).isFinite }) else { throw VivoChemistryError.invalid("multipolar PME source outside FP32 range") }
            input[i] = .init(positionCharge: .init(Float(u[0]),Float(u[1]),Float(u[2]),Float(source.chargeE)),
                dipole: .init(Float(mu[0]),Float(mu[1]),Float(mu[2]),0),
                q0: .init(Float(transformed[0][0]),Float(transformed[0][1]),Float(transformed[0][2]),Float(transformed[1][1])),
                q1: .init(Float(transformed[1][2]),Float(transformed[2][2]),0,0))
        }
        guard let buffer = queue.makeCommandBuffer() else { throw VivoChemistryError.resourceLimit("multipolar PME command buffer") }
        buffer.label = "NumiVivo.PME.multipole.reciprocal"
        try VivoPMEEngine.dispatch(pipelines: pipelines,.mdPMEClearGrid,commandBuffer: buffer,buffers: [gridA],command: &command,elements: count)
        if !sources.isEmpty {
            try VivoPMEEngine.dispatch(pipelines: pipelines,.mdPMEMultipoleSpread,commandBuffer: buffer,
                buffers: [sourceBuffer,gridA],command: &command,elements: sources.count)
        }
        var current = gridA,scratch = gridB
        try VivoPMEEngine.fft3D(pipelines: pipelines,commandBuffer: buffer,current: &current,scratch: &scratch,inverse: false,command: &command)
        guard let influence = pipelines[.mdPMEMultipoleInfluence],let encoder = buffer.makeComputeCommandEncoder() else {
            throw VivoChemistryError.resourceLimit("multipolar PME influence encoder")
        }
        encoder.setComputePipelineState(influence.state)
        encoder.setBuffer(current,offset: 0,index: 0);encoder.setBuffer(scratch,offset: 0,index: 1);encoder.setBuffer(modeBuffer,offset: 0,index: 2)
        var modes = SIMD4<UInt32>(UInt32(cfg.reciprocalHalfWidths[0]),UInt32(cfg.reciprocalHalfWidths[1]),UInt32(cfg.reciprocalHalfWidths[2]),0)
        encoder.setBytes(&modes,length: 16,index: 3);encoder.setBytes(&command,length: MemoryLayout<VivoPMEMetalCommand>.stride,index: 4)
        encoder.dispatchThreads(influence.gridSize(for: count),threadsPerThreadgroup: influence.threadgroupSize(for: count));encoder.endEncoding()
        swap(&current,&scratch)
        try VivoPMEEngine.fft3D(pipelines: pipelines,commandBuffer: buffer,current: &current,scratch: &scratch,inverse: true,command: &command)
        try VivoPMEEngine.dispatch(pipelines: pipelines,.mdPMEScaleInverse,commandBuffer: buffer,buffers: [current],command: &command,elements: count)
        if !sources.isEmpty {
            try VivoPMEEngine.dispatch(pipelines: pipelines,.mdPMEMultipoleGather,commandBuffer: buffer,
                buffers: [sourceBuffer,current,resultBuffer],command: &command,elements: sources.count)
        }
        buffer.commit();buffer.waitUntilCompleted()
        guard buffer.status == .completed,buffer.error == nil else { throw VivoChemistryError.convergence("multipolar PME GPU execution: \(String(describing: buffer.error))") }
        let values = modeBuffer.contents().assumingMemoryBound(to: Float.self)
        var energy = 0.0,strain = VivoQMMatrix(3,3)
        for i in 0..<count {
            let v = (0..<7).map { Double(values[i*8+$0]) }
            energy += v[0];strain[0,0] += v[1];strain[0,1] += v[2];strain[1,0] += v[2]
            strain[0,2] += v[3];strain[2,0] += v[3];strain[1,1] += v[4];strain[1,2] += v[5];strain[2,1] += v[5];strain[2,2] += v[6]
        }
        let output = resultBuffer.contents().assumingMemoryBound(to: Float.self)
        var forces: [VivoVector3D] = [],derivatives: [[Double]] = []
        for (i,source) in sources.enumerated() {
            let v = (0..<13).map { Double(output[i*16+$0]) }
            let force = gridVectors[0]*v[0]+gridVectors[1]*v[1]+gridVectors[2]*v[2]
            let lambdaMu = gridVectors[0]*v[4]+gridVectors[1]*v[5]+gridVectors[2]*v[6]
            let lambdaGrid = [[v[7],v[8]/2,v[9]/2],[v[8]/2,v[10],v[11]/2],[v[9]/2,v[11]/2,v[12]]]
            var lambdaLab = [[Double]](repeating: [0,0,0],count: 3)
            for a in 0..<3 { for b in 0..<3 { for c in 0..<3 { for d in 0..<3 {
                lambdaLab[a][b] += transform[c][a]*lambdaGrid[c][d]*transform[d][b]
            } } } }
            derivatives.append([v[3],lambdaMu.x,lambdaMu.y,lambdaMu.z,lambdaLab[0][0],2*lambdaLab[0][1],
                2*lambdaLab[0][2],lambdaLab[1][1],2*lambdaLab[1][2],lambdaLab[2][2]])
            forces.append(force)
            let mu = xyz(source.dipoleEBohr),lm = xyz(lambdaMu),q = source.secondMomentsEBohr2
            let matrix = [[q[0],q[1],q[2]],[q[1],q[3],q[4]],[q[2],q[4],q[5]]]
            for a in 0..<3 { for b in 0..<3 {
                strain[a,b] -= lm[a]*mu[b]
                for c in 0..<3 { strain[a,b] -= 2*lambdaLab[a][c]*matrix[b][c] }
            } }
        }
        guard energy.isFinite,forces.allSatisfy(\.isFinite),derivatives.joined().allSatisfy(\.isFinite),strain.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite multipolar PME energy or derivatives")
        }
        let modeCount = cfg.reciprocalHalfWidths.reduce(1) { $0*(2*$1+1) }-1
        return .init(energyHartree: energy,forcesHartreePerBohr: forces,momentDerivatives: derivatives,
            affineStrainDerivativeHartree: strain,modeCount: modeCount)
    }
}

public extension VivoReciprocalElectrostaticOperator {
    static func metalPME(configuration: VivoMultipoleMeshConfiguration,maximumSources: Int,
                         device: MTLDevice? = nil,budget: VivoChemistryBudget = .init()) async throws -> Self {
        guard let device = device ?? MTLCreateSystemDefaultDevice() else { throw VivoChemistryError.unsupported("Metal multipolar PME device unavailable") }
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let resource = try await VivoMultipolePMEResource.make(mesh: configuration,capacity: maximumSources,device: device,catalog: catalog,budget: budget)
        return .init(configuration: configuration,evaluate: { sources,cell,cfg,budget in
            try resource.evaluate(sources,cell: cell,configuration: cfg,budget: budget)
        })
    }
}
