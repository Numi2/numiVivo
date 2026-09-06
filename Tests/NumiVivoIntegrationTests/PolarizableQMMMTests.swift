import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PolarizableQMMMTests {
    private func fixture() -> (VivoElectronicSystem,VivoPeriodicCell,VivoPeriodicQMMMConfiguration) {
        let source=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0.1,0.2,0.3)),.init(atomicNumber:1,positionBohr:.init(1.5,0.1,-0.1))],
            alphaElectrons:1,betaElectrons:1,pointCharges:[
                .init(chargeE:0.2,positionBohr:.init(4.5,2,1),classicalParticleIndex:0),
                .init(chargeE:-0.2,positionBohr:.init(-3,4.1,-1),classicalParticleIndex:1),
                .init(chargeE:0,positionBohr:.init(5.8,2.2,1.2),classicalParticleIndex:2),
                .init(chargeE:0,positionBohr:.init(4.6,3.4,1.3),classicalParticleIndex:3)])
        let cell=VivoPeriodicCell(a:.init(1.2,0,0),b:.init(0.1,1.1,0),c:.init(0,0.1,1))
        let unit=pow(VivoAtomicUnits.bohrInNM,3)
        let polarization=VivoInducedDipoleConfiguration(sites:[
            .init(particleIndex:0,principalPolarizabilitiesNM3:.init(0.8,1.1,1.4)*unit,zParticle:2,xParticle:3),
            .init(particleIndex:1,principalPolarizabilitiesNM3:.init(1,1,1)*unit)],
            pairs:[.init(firstParticle:0,secondParticle:1,screeningPerNM3:0.003/unit,permanentFieldScale:0.8,mutualInductionScale:0.9)],
            parameterProvenance:"synthetic derivative fixture; not force-field parameters")
        var cfg=VivoPeriodicQMMMConfiguration(ewald:.init(alphaPerBohr:0.3,realCutoffBohr:22,reciprocalHalfWidths:[8,8,8]),
            exactNearRadiusBohr:3,exactNearSwitchOffBohr:8,polarization:polarization)
        cfg.scf.energyToleranceHartree=1e-12;cfg.scf.densityTolerance=1e-10;cfg.scf.commutatorTolerance=1e-10
        return (source,cell,cfg)
    }
    @Test func mutualDensityResponseAndAnisotropicFramesDifferentiateOneEnergy() throws {
        let (source,cell,cfg)=fixture(),basis=VivoGaussianBasis.hydrogenSTO3G(nucleusIndices:[0,1])
        let charges=[0.12,-0.12,0.0,0.0],h=1e-4
        func evaluate(_ system: VivoElectronicSystem,_ cell: VivoPeriodicCell) throws -> VivoPeriodicHartreeFockResult {
            try VivoPeriodicHartreeFock.evaluate(system:system,basis:basis,cell:cell,configuration:cfg,physicalPermanentChargesE:charges)
        }
        let result=try evaluate(source,cell),response=try #require(result.polarization)
        #expect(response.maximumResponseResidual<1e-10 && response.minimumResponsePivot>0)
        #expect(response.physicalPermanentChargesE==charges)
        var maximum=0.0
        for (center,axis) in [(0,0),(1,2),(2,0),(3,1),(4,1),(5,2)] {
            var plus=source,minus=source
            if center<2 { plus.nuclei[center].positionBohr[axis]+=h;minus.nuclei[center].positionBohr[axis]-=h }
            else { plus.pointCharges[center-2].positionBohr[axis]+=h;minus.pointCharges[center-2].positionBohr[axis]-=h }
            let numerical=try -(evaluate(plus,cell).reference.energyHartree-evaluate(minus,cell).reference.energyHartree)/(2*h)
            let actual=center<2 ? result.nucleusForcesHartreePerBohr[center]:result.pointChargeForcesHartreePerBohr[center-2]
            maximum=max(maximum,abs(numerical-[actual.x,actual.y,actual.z][axis]))
        }
        #expect(maximum<1e-7)
        func deform(_ p: VivoVector3D,_ sign: Double) -> VivoVector3D { .init(p.x+sign*h*p.y,p.y,p.z) }
        var plus=source,minus=source
        for i in source.nuclei.indices { plus.nuclei[i].positionBohr.x+=h*source.nuclei[i].positionBohr.y;minus.nuclei[i].positionBohr.x-=h*source.nuclei[i].positionBohr.y }
        for i in source.pointCharges.indices { plus.pointCharges[i].positionBohr.x+=h*source.pointCharges[i].positionBohr.y;minus.pointCharges[i].positionBohr.x-=h*source.pointCharges[i].positionBohr.y }
        let cp=VivoPeriodicCell(a:deform(cell.a,1),b:deform(cell.b,1),c:deform(cell.c,1))
        let cm=VivoPeriodicCell(a:deform(cell.a,-1),b:deform(cell.b,-1),c:deform(cell.c,-1))
        let strain=try (evaluate(plus,cp).reference.energyHartree-evaluate(minus,cm).reference.energyHartree)/(2*h)
        #expect(abs(strain-result.affineStrainDerivativeHartree[0,1])<1e-7)
        if let path=ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] {
            let root=URL(fileURLWithPath:path);try FileManager.default.createDirectory(at:root,withIntermediateDirectories:true)
            try VivoCanonicalJSON.encode(result).write(to:root.appendingPathComponent("mutual-polarizable-hf.json"))
        }
    }
    @Test func classicalInductionIncludesDampingAndRejectsUnstableResponse() throws {
        let (source,cell,cfg)=fixture(),model=try #require(cfg.polarization)
        let result=try VivoInducedDipoles.evaluate(pointCharges:source.pointCharges,cell:cell,configuration:model,ewald:cfg.ewald)
        #expect(result.response.maximumResponseResidual<1e-10)
        let h=1e-4
        var plus=source.pointCharges,minus=plus;plus[2].positionBohr.y+=h;minus[2].positionBohr.y-=h
        let fd=try -(VivoInducedDipoles.evaluate(pointCharges:plus,cell:cell,configuration:model,ewald:cfg.ewald).energyHartree
            - VivoInducedDipoles.evaluate(pointCharges:minus,cell:cell,configuration:model,ewald:cfg.ewald).energyHartree)/(2*h)
        #expect(abs(fd-result.centerForcesHartreePerBohr[2].y)<1e-8)
        var bad=model;bad.minimumResponsePivot=1e9
        #expect(throws:(any Error).self) { _ = try VivoInducedDipoles.evaluate(pointCharges:source.pointCharges,cell:cell,configuration:bad,ewald:cfg.ewald) }
        var missing=source.pointCharges;missing.removeLast()
        #expect(throws:(any Error).self) { _ = try VivoInducedDipoles.evaluate(pointCharges:missing,cell:cell,configuration:model,ewald:cfg.ewald) }
    }
    @Test func canonicalPolarizationCannotBeSilentlyIgnoredByMetalMD() async throws {
        let (source,_,cfg)=fixture(),model=try #require(cfg.polarization)
        let fp=try VivoCanonicalJSON.fingerprint(Data("polarization-runtime-fixture".utf8))
        let particles=source.pointCharges.enumerated().map { i,q in VivoClassicalParticle(index:UInt32(i),atomIndex:UInt32(i),
            typeIdentifier:"X",massDa:12,chargeE:q.chargeE,sigmaNM:0,epsilonKJPerMol:0) }
        let system=VivoClassicalSystem(identifier:"polarization-runtime-fixture",structureFingerprint:fp,particles:particles,polarization:model)
        let cell=VivoPeriodicCell(a:.init(1.2,0,0),b:.init(0,1.2,0),c:.init(0,0,1.2))
        let positions=source.pointCharges.map { VivoVector3D($0.positionBohr.x,$0.positionBohr.y,$0.positionBohr.z)*VivoAtomicUnits.bohrInNM }
        let initial=VivoClassicalInitialState(systemFingerprint:try system.fingerprint(),positionsNM:positions,periodicCell:cell)
        let dynamics=VivoMDConfiguration(timeStepPS:1e-5,cutoffNM:0.4,neighborSkinNM:0.1,ensemble:.nve,
            thermostat:.none,targetTemperatureK:nil,frictionPerPS:nil,pmeGridDimensions:[32,32,32])
        #expect(try !VivoMDCapabilityAnalyzer.analyze(system:system,initialState:initial,configuration:dynamics).executable)
        let provider=try VivoMDCandidateForceProvider.inducedDipoles(system:system,ewald:cfg.ewald)
        #expect(try VivoMDCapabilityAnalyzer.analyze(system:system,initialState:initial,configuration:dynamics,forceProvider:provider).executable)
        let runtime=try await VivoMDMetalRuntime.make(system:system,initialState:initial,configuration:dynamics,forceProvider:provider)
        let before=try await runtime.checkpoint()
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:before.positionsNM,periodicCell:before.periodicCell)
        let probe=try await runtime.evaluateHamiltonian(at:geometry)
        #expect(probe.normalizedForceResidual != nil && probe.energyKJPerMol.isFinite)
        #expect(try await runtime.checkpoint()==before)
        #expect(try await runtime.step().committed)
    }
    @Test func nonlinearSiteStrainIncludesPeriodicParentImages() throws {
        let graph=try VivoVirtualSiteGraph(particleCount:4,physicalParticles:[0,1,2],definitions:[
            .init(siteParticle:3,parentParticles:[0,1,2],rule:.localCoordinates(originWeights:[1,0,0],xWeights:[-1,1,0],yWeights:[-1,0,1],positionNM:.init(0.02,0.01,0.03)),provenance:"synthetic")])
        let cell=VivoPeriodicCell(a:.init(2,0,0),b:.init(0,2,0),c:.init(0,0,2))
        let positions:[VivoVector3D]=[.init(1.95,1.95,0),.init(0.08,1.98,0),.init(1.97,0.08,0.05),.zero]
        let state=try graph.construct(positionsNM:positions,periodicCell:cell)
        let target=VivoVector3D(1,1,1),siteForce=(state.positionsNM[3]-target) * -1
        let extra=try graph.affineStrainCorrection(rawForces:[.zero,.zero,.zero,siteForce],state:state,periodicCell:cell)
        let h=1e-5
        func energy(_ t:Double) throws -> Double {
            func deform(_ r:VivoVector3D) -> VivoVector3D { .init(r.x+t*r.y,r.y,r.z) }
            let box=VivoPeriodicCell(a:deform(cell.a),b:deform(cell.b),c:deform(cell.c))
            let output=try graph.construct(positionsNM:positions.map(deform),periodicCell:box)
            return 0.5*(output.positionsNM[3]-target).squaredNorm
        }
        let actual = -siteForce.x*state.positionsNM[3].y+extra[0,1]
        #expect(abs(try (energy(h)-energy(-h))/(2*h)-actual)<1e-7)
    }
}
