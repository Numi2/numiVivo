import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MultipolePMETests {
    private func sources() -> [VivoCartesianMultipole] {
        [.init(positionBohr: .init(0.37,1.23,2.19),chargeE: 0.8,dipoleEBohr: .init(0.13,-0.2,0.08),secondMomentsEBohr2: [0.09,0.03,-0.04,0.12,0.02,-0.05]),
         .init(positionBohr: .init(5.31,4.11,3.24),chargeE: -0.8,dipoleEBohr: .init(-0.1,0.06,-0.04),secondMomentsEBohr2: [-0.03,0.02,0.01,0.07,0.02,0.04])]
    }
    private func cell() -> VivoPeriodicCell {
        .init(a: .init(15,0,0)*VivoAtomicUnits.bohrInNM,b: .init(2,16,0)*VivoAtomicUnits.bohrInNM,c: .init(1,0.5,17)*VivoAtomicUnits.bohrInNM)
    }
    private func configuration() -> VivoPeriodicElectrostaticConfiguration {
        .init(alphaPerBohr: 0.35,realCutoffBohr: 12,reciprocalHalfWidths: [5,5,5])
    }
    @Test func multipolarMeshUsesSharedFFTAndConvergesToDirectEwald() async throws {
        let source = sources(),box = cell(),cfg = configuration()
        let reference = try VivoPeriodicElectrostatics.evaluate(sources: source,cell: box,configuration: cfg)
        var errors: [Double] = []
        for size in [16,64] {
            let mesh = try await VivoReciprocalElectrostaticOperator.metalPME(configuration: .init(gridDimensions: [size,size,size]),maximumSources: source.count)
            let result = try VivoPeriodicElectrostatics.evaluate(sources: source,cell: box,configuration: cfg,reciprocalOperator: mesh)
            errors.append(abs(result.energyHartree-reference.energyHartree))
            if size == 64 {
                #expect(errors.last! < 2e-6)
                for (a,b) in zip(result.forcesHartreePerBohr,reference.forcesHartreePerBohr) { #expect((a-b).norm < 3e-6) }
                for (a,b) in zip(result.momentDerivatives.joined(),reference.momentDerivatives.joined()) { #expect(abs(a-b) < 5e-6) }
                #expect(try result.affineStrainDerivativeHartree.adding(reference.affineStrainDerivativeHartree,scale: -1).frobeniusNorm < 1e-5)
            }
        }
        #expect(errors[1] < errors[0]+2e-7)
    }
    @Test func reciprocalMeshForceAndStrainDifferentiateItsOwnEnergy() async throws {
        let source = sources(),box = cell(),cfg = configuration()
        let mesh = try await VivoReciprocalElectrostaticOperator.metalPME(configuration: .init(gridDimensions: [32,32,32]),maximumSources: source.count)
        let result = try VivoPeriodicElectrostatics.evaluate(sources: source,cell: box,configuration: cfg,reciprocalOperator: mesh)
        let h = 0.01
        var plus = source,minus = source
        plus[0].positionBohr.x += h;minus[0].positionBohr.x -= h
        let fp = try VivoPeriodicElectrostatics.evaluate(sources: plus,cell: box,configuration: cfg,reciprocalOperator: mesh)
        let fm = try VivoPeriodicElectrostatics.evaluate(sources: minus,cell: box,configuration: cfg,reciprocalOperator: mesh)
        #expect(abs(-(fp.energyHartree-fm.energyHartree)/(2*h)-result.forcesHartreePerBohr[0].x) < 5e-6)
        func affine(_ e: Double) throws -> Double {
            var s = source
            for i in s.indices { s[i].positionBohr.x += e*s[i].positionBohr.y }
            func move(_ v: VivoVector3D) -> VivoVector3D { .init(v.x+e*v.y,v.y,v.z) }
            let b = VivoPeriodicCell(a: move(box.a),b: move(box.b),c: move(box.c))
            return try VivoPeriodicElectrostatics.evaluate(sources: s,cell: b,configuration: cfg,reciprocalOperator: mesh).energyHartree
        }
        #expect(abs(try (affine(0.001)-affine(-0.001))/0.002-result.affineStrainDerivativeHartree[0,1]) < 2e-5)
    }
    @Test func fixedClassicalMeshSurvivesCellMovesAndRejectsInvalidWidths() async throws {
        let first = try VivoPMEPlan.make(cell: cell(),cutoffNM: 0.2,tolerance: 1e-5,targetGridSpacingNM: 0.1,fixedGridDimensions: [16,32,16])
        let box = cell(),changed = VivoPeriodicCell(a: box.a*1.05,b: box.b*1.05,c: box.c*1.05)
        let second = try VivoPMEPlan.make(cell: changed,cutoffNM: 0.2,tolerance: 1e-5,targetGridSpacingNM: 0.1,fixedGridDimensions: [16,32,16])
        #expect(first.gridX == second.gridX && first.gridY == second.gridY && first.gridZ == second.gridZ)
        #expect(first.cellVolumeNM3 != second.cellVolumeNM3)
        let mesh = try await VivoReciprocalElectrostaticOperator.metalPME(configuration: .init(gridDimensions: [16,16,16]),maximumSources: 2)
        var cfg = configuration();cfg.reciprocalHalfWidths = [8,5,5]
        #expect(throws: (any Error).self) { _ = try VivoPeriodicElectrostatics.evaluate(sources: sources(),cell: cell(),configuration: cfg,reciprocalOperator: mesh) }
    }
}
