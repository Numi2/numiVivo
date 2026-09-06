import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct AdaptiveQMMMHamiltonianTests {
    private func fixture() throws -> (VivoMolecularStructureDocument,VivoClassicalSystem,VivoMDCandidateGeometry,VivoAdaptiveQMMMConfiguration) {
        let h=try #require(VivoElement.from(symbol:"H"))
        let atoms=(0..<4).map { VivoMolecularAtom(index:UInt32($0),name:"H\($0)",element:h) }
        let positions:[VivoVector3D]=[.init(0.5,0.5,0.5),.init(0.574,0.5,0.5),.init(0.85,0.5,0.5),.init(0.924,0.5,0.5)]
        let rounded=positions.map { VivoVector3D(Double(Float($0.x)),Double(Float($0.y)),Double(Float($0.z))) }
        let cell=VivoPeriodicCell(a:.init(2,0,0),b:.init(0,2,0),c:.init(0,0,2))
        let document=try VivoMolecularStructureDocument(structure:.init(identifier:"adaptive-h2-fixture",atoms:atoms,
            bonds:[.init(atomA:0,atomB:1),.init(atomA:2,atomB:3)],conformers:[.init(positionsNM:rounded)],periodicCell:cell))
        let particles=atoms.map { VivoClassicalParticle(index:$0.index,atomIndex:$0.index,typeIdentifier:"H",massDa:1,chargeE:0,sigmaNM:0,epsilonKJPerMol:0) }
        let system=VivoClassicalSystem(identifier:"adaptive-h2-fixture",structureFingerprint:document.structureFingerprint,particles:particles,
            bonds:[.init(a:0,b:1,lengthNM:0.074,forceConstant:100),.init(a:2,b:3,lengthNM:0.074,forceConstant:100)])
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:rounded,periodicCell:cell)
        let config=VivoAdaptiveQMMMConfiguration(fixedCore:.init(qmAtomIndices:[0,1],alphaElectrons:1,betaElectrons:1),
            centerAtomIndices:[0,1],molecules:[.init(identifier:"solvent-H2",atomIndices:[2,3],promotionReferenceKJPerMol:-2900,
                compensationCoefficientsKJPerMol:[0.3,-0.1],parameterProvenance:"synthetic reference; not equilibrium calibration")],innerRadiusNM:0.2,outerRadiusNM:0.4)
        return (document,system,geometry,config)
    }
    @Test func energyWeightsIncludeTransitionForcesMovingCoreAndCompensation() throws {
        let (document,system,geometry,cfg)=try fixture()
        let plan=try VivoAdaptiveQMMMPlan.compile(document:document,system:system,configuration:cfg)
        func evaluate(_ geometry:VivoMDCandidateGeometry) throws -> VivoAdaptiveHamiltonianResult {
            let weights=try plan.weights(document:document,system:system,geometry:geometry)
            let values=weights.partitions.map { partition in
                let constant=partition.selectedMoleculeIndices.isEmpty ? 1.0:3.0
                return VivoMDHamiltonianEvaluation(systemFingerprint:plan.systemFingerprint,configurationFingerprint:plan.structureFingerprint,
                    evaluatedGeometry:geometry,energyKJPerMol:constant+partition.energyReferenceKJPerMol,
                    physicalParticleForcesKJPerMolNM:Array(repeating:.zero,count:4),normalizedForceResidual:0)
            }
            return try VivoAdaptiveHamiltonian.combine(weights:weights,particleCount:4,evaluations:values)
        }
        let result=try evaluate(geometry),h=1e-6
        #expect(result.weightState.partitions.count==2 && result.weightState.switchingMoleculeIndices==[0])
        for i in [0,2] {
            var plus=geometry.particlePositionsNM,minus=plus;plus[i].x+=h;minus[i].x-=h
            let p=try VivoMDCandidateGeometry(particlePositionsNM:plus,periodicCell:geometry.periodicCell)
            let m=try VivoMDCandidateGeometry(particlePositionsNM:minus,periodicCell:geometry.periodicCell)
            let fd=try -(evaluate(p).energyKJPerMol-evaluate(m).energyKJPerMol)/(2*h)
            #expect(abs(fd-result.physicalParticleForcesKJPerMolNM[i].x)<1e-6)
        }
        #expect(result.physicalParticleForcesKJPerMolNM.reduce(.zero,+).norm<1e-10)
        var limited=cfg;limited.maximumSwitchingMolecules=0
        let small=try VivoAdaptiveQMMMPlan.compile(document:document,system:system,configuration:limited)
        #expect(throws:(any Error).self) { _ = try small.weights(document:document,system:system,geometry:geometry) }
        var constrained=system;constrained.constraints=[.init(a:2,b:3,distanceNM:0.074)]
        #expect(throws:(any Error).self) { _ = try VivoAdaptiveQMMMPlan.compile(document:document,system:constrained,configuration:cfg) }
    }
    @Test func completeAdaptivePartitionsRunThroughExistingMetalHamiltonians() async throws {
        let (document,system,geometry,cfg)=try fixture()
        let dynamics=VivoMDConfiguration(timeStepPS:1e-5,cutoffNM:0.6,neighborSkinNM:0.1,ensemble:.nve,
            thermostat:.none,targetTemperatureK:nil,frictionPerPS:nil,pmeGridDimensions:[32,32,32])
        let basis=VivoAdaptiveBasisSpecification(atomBasis:.hydrogenSTO3G(nucleusIndices:[0,1,2,3]))
        let electronic=VivoQMMMDynamicsElectronicConfiguration(basis:.hydrogenSTO3G(nucleusIndices:[0,1]))
        let setup=try await VivoAdaptiveQMMMDynamics.make(document:document,system:system,initialGeometry:geometry,
            dynamics:dynamics,adaptive:cfg,basis:basis,electronic:electronic)
        let provider=setup.provider,base=setup.retainedBaselineSystem
        let initial=VivoClassicalInitialState(systemFingerprint:try base.fingerprint(),positionsNM:geometry.particlePositionsNM,periodicCell:geometry.periodicCell)
        let runtime=try await VivoMDMetalRuntime.make(system:base,initialState:initial,configuration:dynamics,forceProvider:provider)
        let before=try await runtime.checkpoint(),first=try await runtime.evaluateHamiltonian(at:geometry)
        #expect(first.energyKJPerMol.isFinite && first.normalizedForceResidual != nil)
        #expect(try await runtime.checkpoint()==before)
        // Re-evaluating after cache reuse must not mutate a reference offset or state.
        let repeated=try await runtime.evaluateHamiltonian(at:geometry)
        #expect(first==repeated)
        #expect(try await runtime.step().committed)
        let after=try await runtime.checkpoint();#expect(after.acceptedStep==1)
        let restored=try await VivoMDMetalRuntime.restore(system:base,configuration:dynamics,checkpoint:after,forceProvider:provider)
        #expect(try await restored.checkpoint()==after)
    }
}
