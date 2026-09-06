import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct PeriodicQMMMHamiltonianTests {
    private func identity(_ label: String) throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(Data(label.utf8)) }
    private func fixture() throws -> (VivoMolecularStructureDocument,VivoClassicalSystem) {
        let c = try #require(VivoElement.from(symbol: "C"))
        let atoms = (0..<4).map { VivoMolecularAtom(index: UInt32($0),name: "C\($0)",element: c) }
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "ownership",atoms: atoms,
            bonds: [.init(atomA: 0,atomB: 1),.init(atomA: 1,atomB: 2),.init(atomA: 2,atomB: 3)]))
        let charges = [0.1,-0.1,0.25,-0.25]
        let system = VivoClassicalSystem(identifier: "ownership",structureFingerprint: document.structureFingerprint,
            particles: atoms.map { .init(index: $0.index,atomIndex: $0.index,typeIdentifier: "C",massDa: 12,
                chargeE: charges[Int($0.index)],sigmaNM: 0.3,epsilonKJPerMol: 0.2) },
            bonds: [.init(a: 0,b: 1,lengthNM: 0.15,forceConstant: 100),.init(a: 1,b: 2,lengthNM: 0.15,forceConstant: 100),
                    .init(a: 2,b: 3,lengthNM: 0.15,forceConstant: 100)],
            angles: [.init(a: 0,b: 1,c: 2,angleRadians: 1.9,forceConstant: 10),.init(a: 1,b: 2,c: 3,angleRadians: 1.9,forceConstant: 10)],
            torsions: [.init(a: 0,b: 1,c: 2,d: 3,periodicity: 3,phaseRadians: 0,barrierKJPerMol: 2)],
            nonbondedExceptions: [.init(a: 0,b: 2,coulombScale: 0.5,lennardJonesScale: 0.3)])
        return (document,system)
    }
    private func configuration() -> VivoQMMMHamiltonianConfiguration {
        .init(region: .init(qmAtomIndices: [1,0],alphaElectrons: 7,betaElectrons: 6),boundary: .periodicElectrostatic,
            boundaryChargeTransfers: [.init(donorParticle: 2,recipients: [.init(particleIndex: 3,weight: 1)])])
    }
    @Test func retainedTermCompilerOwnsEveryTermAndSeparatesChargeChannels() throws {
        let (document,system) = try fixture(),before = system
        let plan = try VivoQMMMHamiltonianPlan.compile(document: document,system: system,configuration: configuration())
        #expect(system == before)
        #expect(plan.bondOwners == [.electronic,.retainedClassical,.retainedClassical])
        #expect(plan.angleOwners == [.electronic,.retainedClassical])
        #expect(plan.torsionOwners == [.retainedClassical])
        #expect(plan.retainedSystem.bonds.map(\.a) == [1,2])
        #expect(plan.retainedSystem.particles.map(\.chargeE) == [0,0,0.25,-0.25])
        #expect(plan.embeddingChargesE == [0,0,0,0])
        #expect(plan.retainedSystem.nonbondedExceptions.contains { $0.a == 0 && $0.b == 1 && $0.lennardJonesScale == 0 })
        #expect(plan.retainedSystem.nonbondedExceptions.contains(system.nonbondedExceptions[0]))
        #expect(plan.expectedQMChargeE == 0 && plan.totalMMChargeE == 0)
        let packed = try VivoMDSystemPacker.pack(plan.retainedSystem)
        #expect(packed.bonds.count == 2 && packed.angles.count == 1 && packed.torsions.count == 1)
        let roundtrip = try VivoCanonicalJSON.decode(VivoQMMMHamiltonianPlan.self,from: VivoCanonicalJSON.encode(plan))
        try roundtrip.validate(document: document,source: system)
    }
    @Test func chargedUnmappedBoundaryAndUndeclaredConstraintsFail() throws {
        let (document,source) = try fixture()
        var cfg = configuration();cfg.boundaryChargeTransfers = []
        #expect(throws: (any Error).self) { _ = try VivoQMMMHamiltonianPlan.compile(document: document,system: source,configuration: cfg) }
        cfg = configuration();cfg.region.alphaElectrons -= 1
        #expect(throws: (any Error).self) { _ = try VivoQMMMHamiltonianPlan.compile(document: document,system: source,configuration: cfg) }
        var constrained = source;constrained.constraints = [.init(a: 0,b: 1,distanceNM: 0.15)]
        cfg = configuration()
        #expect(throws: (any Error).self) { _ = try VivoQMMMHamiltonianPlan.compile(document: document,system: constrained,configuration: cfg) }
        cfg.constraintPolicy = .preserveExplicitManifold
        let result = try VivoQMMMHamiltonianPlan.compile(document: document,system: constrained,configuration: cfg)
        #expect(result.retainedSystem.constraints == constrained.constraints)
        cfg.maximumGeneratedExceptions = 0
        #expect(throws: (any Error).self) { _ = try VivoQMMMHamiltonianPlan.compile(document: document,system: constrained,configuration: cfg) }
    }
    @Test func allMMOwnershipPreservesPhysicalForceField() throws {
        let (document,system) = try fixture()
        let plan = try VivoQMMMHamiltonianPlan.compile(document: document,system: system,
            configuration: .init(region: .init(qmAtomIndices: [],alphaElectrons: 0,betaElectrons: 0),boundary: .periodicElectrostatic))
        #expect(plan.retainedSystem.particles == system.particles && plan.retainedSystem.bonds == system.bonds)
        #expect(plan.retainedSystem.angles == system.angles && plan.retainedSystem.torsions == system.torsions)
        #expect(plan.retainedSystem.nonbondedExceptions == system.nonbondedExceptions)
        #expect(plan.embeddingChargesE == system.particles.map(\.chargeE))
    }
    @Test func nonlinearAndNestedSitesHaveConsistentJacobianVelocitiesAndForces() throws {
        let graph = try VivoVirtualSiteGraph(particleCount: 6,physicalParticles: [0,1,2],definitions: [
            .init(siteParticle: 3,parentParticles: [0,1,2],rule: .outOfPlane(weight12: 0.2,weight13: 0.3,crossWeightPerNM: 0.4),provenance: "fixture"),
            .init(siteParticle: 4,parentParticles: [0,1,2],rule: .localCoordinates(originWeights: [1,0,0],xWeights: [-1,1,0],yWeights: [-1,0,1],positionNM: .init(0.1,0.2,-0.1)),provenance: "fixture"),
            .init(siteParticle: 5,parentParticles: [3,4],rule: .linear(weights: [0.4,0.6]),provenance: "fixture")])
        let positions: [VivoVector3D] = [.init(0.2,0.1,0.3),.init(0.4,0.2,0.2),.init(0.1,0.4,0.5),.zero,.zero,.zero]
        let state = try graph.construct(positionsNM: positions),direction = state.positionsNM[5]-state.positionsNM[0]
        var forces = [VivoVector3D](repeating: .zero,count: 6);forces[5] = direction * -1;forces[0] = direction
        let mapped = try graph.redistribute(rawForces: forces,state: state)
        func displaced(_ v: VivoVector3D,_ axis: Int,_ amount: Double) -> VivoVector3D {
            var v = v;if axis == 0 { v.x += amount } else if axis == 1 { v.y += amount } else { v.z += amount };return v
        }
        let h = 1e-6
        for i in 0..<3 { for axis in 0..<3 {
            var plus = positions,minus = positions
            plus[i] = displaced(plus[i],axis,h);minus[i] = displaced(minus[i],axis,-h)
            let p = try graph.construct(positionsNM: plus).positionsNM,m = try graph.construct(positionsNM: minus).positionsNM
            let fd = -0.5*((p[5]-p[0]).squaredNorm-(m[5]-m[0]).squaredNorm)/(2*h)
            #expect(abs(fd-[mapped[i].x,mapped[i].y,mapped[i].z][axis]) < 1e-8)
        } }
        #expect(mapped[3...5].allSatisfy { $0 == .zero })
        #expect(mapped.reduce(.zero,+).norm < 1e-12)
        #expect(zip(state.positionsNM,mapped).reduce(VivoVector3D.zero) { $0+$1.0.cross($1.1) }.norm < 1e-12)
        let velocities: [VivoVector3D] = [.init(0.1,0.2,0.3),.init(-0.2,0.1,0.2),.init(0.3,0.1,-0.1),.zero,.zero,.zero]
        let result = try graph.velocities(velocities,state: state)
        let moved = try graph.construct(positionsNM: zip(positions,velocities).map { $0+$1*h }).positionsNM
        #expect(((moved[5]-state.positionsNM[5])/h-result[5]).norm < 1e-6)
        var degenerate = positions;degenerate[2] = degenerate[1]
        #expect(throws: (any Error).self) { _ = try graph.construct(positionsNM: degenerate) }
    }
    @Test func multipolarEwaldDifferentiatesSourcesAndLattice() throws {
        let cell = VivoPeriodicCell(a: .init(1,0,0),b: .init(0.25,0.95,0),c: .init(0.1,0.15,1.1))
        let sources = [VivoCartesianMultipole(positionBohr: .init(1,2,3),chargeE: 0.8,dipoleEBohr: .init(0.12,-0.2,0.13),secondMomentsEBohr2: [0.15,0.05,-0.02,0.18,0.04,0.12]),
            .init(positionBohr: .init(5,4,3.5),chargeE: -0.8,dipoleEBohr: .init(-0.1,0.03,0.2),secondMomentsEBohr2: [-0.12,0.02,0.04,0.02,-0.03,0.06])]
        var cfg = VivoPeriodicElectrostaticConfiguration(alphaPerBohr: 0.35,realCutoffBohr: 22,reciprocalHalfWidths: [9,9,9],chargeConvention: .uniformNeutralizingBackground)
        let actual = try VivoPeriodicElectrostatics.evaluate(sources: sources,cell: cell,configuration: cfg),h = 1e-5
        var plus = sources,minus = sources;plus[0].positionBohr.x += h;minus[0].positionBohr.x -= h
        let force = try -(VivoPeriodicElectrostatics.evaluate(sources: plus,cell: cell,configuration: cfg).energyHartree
            - VivoPeriodicElectrostatics.evaluate(sources: minus,cell: cell,configuration: cfg).energyHartree)/(2*h)
        #expect(abs(force-actual.forcesHartreePerBohr[0].x) < 1e-8)
        plus = sources;minus = sources;plus[1].secondMomentsEBohr2[1] += h;minus[1].secondMomentsEBohr2[1] -= h
        let moment = try (VivoPeriodicElectrostatics.evaluate(sources: plus,cell: cell,configuration: cfg).energyHartree
            - VivoPeriodicElectrostatics.evaluate(sources: minus,cell: cell,configuration: cfg).energyHartree)/(2*h)
        #expect(abs(moment-actual.momentDerivatives[1][5]) < 1e-8)
        func strain(_ p: VivoVector3D,_ sign: Double) -> VivoVector3D { .init(p.x+sign*h*p.y,p.y,p.z) }
        let cp = VivoPeriodicCell(a: strain(cell.a,1),b: strain(cell.b,1),c: strain(cell.c,1))
        let cm = VivoPeriodicCell(a: strain(cell.a,-1),b: strain(cell.b,-1),c: strain(cell.c,-1))
        plus = sources;minus = sources
        for i in sources.indices { plus[i].positionBohr = strain(sources[i].positionBohr,1);minus[i].positionBohr = strain(sources[i].positionBohr,-1) }
        let derivative = try (VivoPeriodicElectrostatics.evaluate(sources: plus,cell: cp,configuration: cfg).energyHartree
            - VivoPeriodicElectrostatics.evaluate(sources: minus,cell: cm,configuration: cfg).energyHartree)/(2*h)
        #expect(abs(derivative-actual.affineStrainDerivativeHartree[0,1]) < 1e-8)
        cfg.alphaPerBohr = 0.25
        let alternate = try VivoPeriodicElectrostatics.evaluate(sources: sources,cell: cell,configuration: cfg)
        #expect(abs(alternate.energyHartree-actual.energyHartree) < 1e-8)
        plus = sources;plus[1].positionBohr = plus[1].positionBohr+cell.b/VivoAtomicUnits.bohrInNM
        #expect(abs(try VivoPeriodicElectrostatics.evaluate(sources: plus,cell: cell,configuration: cfg).energyHartree-alternate.energyHartree) < 1e-12)
        #expect(actual.forcesHartreePerBohr.reduce(.zero,+).norm < 1e-12)
    }
    @Test func periodicHFIncludesDensityResponsePulayAndNearSwitchDerivatives() throws {
        let system = VivoElectronicSystem(nuclei: [.init(atomicNumber: 1,positionBohr: .init(0.1,0.2,0.3)),.init(atomicNumber: 1,positionBohr: .init(1.5,0.1,-0.1))],
            alphaElectrons: 1,betaElectrons: 1,pointCharges: [.init(chargeE: 0.2,positionBohr: .init(4.5,2,1)),.init(chargeE: -0.2,positionBohr: .init(-3,4.1,-1))])
        let basis = VivoGaussianBasis.hydrogenSTO3G(nucleusIndices: [0,1]),cell = VivoPeriodicCell(a: .init(1.2,0,0),b: .init(0.1,1.1,0),c: .init(0,0.1,1))
        var cfg = VivoPeriodicQMMMConfiguration(ewald: .init(alphaPerBohr: 0.3,realCutoffBohr: 22,reciprocalHalfWidths: [8,8,8]),exactNearRadiusBohr: 3,exactNearSwitchOffBohr: 8)
        cfg.scf.energyToleranceHartree = 1e-12;cfg.scf.densityTolerance = 1e-10;cfg.scf.commutatorTolerance = 1e-10
        let actual = try VivoPeriodicHartreeFock.evaluate(system: system,basis: basis,cell: cell,configuration: cfg),h = 1e-4
        var maximum = 0.0
        for center in 0..<4 { for axis in 0..<3 {
            var plus = system,minus = system
            if center < 2 { plus.nuclei[center].positionBohr[axis] += h;minus.nuclei[center].positionBohr[axis] -= h }
            else { plus.pointCharges[center-2].positionBohr[axis] += h;minus.pointCharges[center-2].positionBohr[axis] -= h }
            let eplus = try VivoPeriodicHartreeFock.evaluate(system: plus,basis: basis,cell: cell,configuration: cfg).reference.energyHartree
            let eminus = try VivoPeriodicHartreeFock.evaluate(system: minus,basis: basis,cell: cell,configuration: cfg).reference.energyHartree
            let force = center < 2 ? actual.nucleusForcesHartreePerBohr[center] : actual.pointChargeForcesHartreePerBohr[center-2]
            maximum = max(maximum,abs((eplus-eminus)/(2*h)+[force.x,force.y,force.z][axis]))
        } }
        #expect(maximum < 1e-7)
        #expect((actual.nucleusForcesHartreePerBohr+actual.pointChargeForcesHartreePerBohr).reduce(.zero,+).norm < 1e-9)
        if let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] {
            let root = URL(fileURLWithPath: path);try FileManager.default.createDirectory(at: root,withIntermediateDirectories: true)
            try VivoCanonicalJSON.encode(actual).write(to: root.appendingPathComponent("periodic-hf-analytic.json"),options: .atomic)
        }
    }
    @Test func actualMetalRuntimeCommitsBOForcesAndRejectsMissingProviderOnRestore() async throws {
        let h = try #require(VivoElement.from(symbol: "H")),he = try #require(VivoElement.from(symbol: "He"))
        let atoms = [VivoMolecularAtom(index: 0,name: "H1",element: h),.init(index: 1,name: "H2",element: h),
                     .init(index: 2,name: "M1",element: he),.init(index: 3,name: "M2",element: he)]
        let cell = VivoPeriodicCell(a: .init(2,0,0),b: .init(0,2,0),c: .init(0,0,2))
        let positions: [VivoVector3D] = [.init(0.5,0.5,0.5),.init(0.574,0.5,0.5),.init(0.7,0.7,0.7),.init(1.1,0.8,1)]
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "native-bo",atoms: atoms,
            bonds: [.init(atomA: 0,atomB: 1)],conformers: [.init(positionsNM: positions)],periodicCell: cell))
        let charges = [0.0,0.0,0.2,-0.2]
        let system = VivoClassicalSystem(identifier: "native-bo",structureFingerprint: document.structureFingerprint,
            particles: atoms.map { .init(index: $0.index,atomIndex: $0.index,typeIdentifier: $0.element.symbol,massDa: 2,
                chargeE: charges[Int($0.index)],sigmaNM: 0,epsilonKJPerMol: 0) },
            bonds: [.init(a: 0,b: 1,lengthNM: 0.074,forceConstant: 100)])
        let plan = try VivoQMMMHamiltonianPlan.compile(document: document,system: system,
            configuration: .init(region: .init(qmAtomIndices: [0,1],alphaElectrons: 1,betaElectrons: 1),boundary: .periodicElectrostatic))
        let provider = try VivoMDCandidateForceProvider.hartreeFock(document: document,sourceSystem: system,plan: plan,
            configuration: .init(basis: .hydrogenSTO3G(nucleusIndices: [0,1])))
        let dynamics = VivoMDConfiguration(timeStepPS: 0.00001,cutoffNM: 0.6,neighborSkinNM: 0.1,electrostatics: .pme,
            ensemble: .nve,thermostat: .none,targetTemperatureK: nil,frictionPerPS: nil)
        let initial = VivoClassicalInitialState(systemFingerprint: try plan.retainedSystem.fingerprint(),positionsNM: positions,periodicCell: cell)
        let runtime = try await VivoMDMetalRuntime.make(system: plan.retainedSystem,initialState: initial,configuration: dynamics,forceProvider: provider)
        let before = try await runtime.checkpoint()
        #expect(before.configurationFingerprint != (try dynamics.fingerprint()))
        let step = try await runtime.step();#expect(step.committed)
        let after = try await runtime.checkpoint();#expect(after.acceptedStep == 1 && after.positionsNM != before.positionsNM)
        let observables = try await runtime.observables();#expect(observables.potentialEnergyKJPerMol < -1000)
        let restored = try await VivoMDMetalRuntime.restore(system: plan.retainedSystem,configuration: dynamics,checkpoint: after,forceProvider: provider)
        #expect(try await restored.checkpoint() == after)
        do {
            _ = try await VivoMDMetalRuntime.restore(system: plan.retainedSystem,configuration: dynamics,checkpoint: after)
            Issue.record("BO checkpoint resumed without its Hamiltonian provider")
        } catch {}
        let rejecting = try VivoMDCandidateForceProvider(fingerprint: identity("rejecting-bo"),retainedSystemFingerprint: plan.retainedSystem.fingerprint(),
            boundary: .periodicElectrostatic,supportsCellMoves: false,maximumAcceptedResidual: 1e-8,molecularConnectivitySystem: system,
            evaluate: { _ in throw VivoChemistryError.convergence("injected electronic failure") })
        let failed = try await VivoMDMetalRuntime.make(system: plan.retainedSystem,initialState: initial,configuration: dynamics,forceProvider: rejecting)
        let failureBefore = try await failed.checkpoint()
        do { _ = try await failed.step();Issue.record("failed electronic candidate committed") } catch {}
        #expect(try await failed.checkpoint() == failureBefore)
    }
}
