import Foundation
import Testing
@preconcurrency import Metal
@testable import NumiVivoKit

@Suite(.serialized) struct PrecisionMetalSamplingTests {
    private func specification(k: Double) throws -> VivoNuclearMetalSpecification {
        let fp = try PrecisionSamplingFixtures.id("native-dimer-structure")
        let system = VivoClassicalSystem(identifier: "native-dimer-\(k)",structureFingerprint: fp,
            particles: [.init(index: 0,atomIndex: 1,typeIdentifier: "A",massDa: 12,chargeE: 0,sigmaNM: 0,epsilonKJPerMol: 0),
                        .init(index: 1,atomIndex: 0,typeIdentifier: "B",massDa: 1,chargeE: 0,sigmaNM: 0,epsilonKJPerMol: 0)],
            bonds: [.init(a: 0,b: 1,lengthNM: 1,forceConstant: k)])
        let dynamics = VivoMDConfiguration(timeStepPS: 0.001,cutoffNM: 2,neighborSkinNM: 0,electrostatics: .cutoff,
            ensemble: .nve,thermostat: .none,targetTemperatureK: nil,frictionPerPS: nil,neighborListEnabled: false)
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),positionsNM: PrecisionSamplingFixtures.pairPositions(1))
        return .init(system: system,initialState: initial,dynamics: dynamics)
    }
    private func record<T: Encodable>(_ value: T,_ name: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let root = URL(fileURLWithPath: path); try FileManager.default.createDirectory(at: root,withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: root.appendingPathComponent(name+".json"),options: .atomic)
    }
    @Test func batchedMetalInferenceMatchesFP64EnergyAndAnalyticForce() async throws {
        _ = try #require(MTLCreateSystemDefaultDevice())
        let (model,_,_,_) = try await PrecisionSamplingFixtures.trainedPair()
        let backend = try VivoReactiveEnergyForceBackend.metal(model: model,maximumBatchSize: 16)
        let frames = [0.86,0.94,1.08,1.14].map(PrecisionSamplingFixtures.pairPositions)
        let gpu = try await backend.predict(frames)
        #expect(gpu.count == frames.count)
        for i in frames.indices {
            let cpu = try VivoReactiveDeltaSurrogate.predict(model,positionsNM: frames[i])
            #expect(abs(cpu.deltaEnergyKJPerMol-gpu[i].deltaEnergyKJPerMol) < 0.001)
            for atom in frames[i].indices { #expect((cpu.deltaForcesKJPerMolNM[atom]-gpu[i].deltaForcesKJPerMolNM[atom]).norm < 0.02) }
            #expect(cpu.eligible == gpu[i].eligible)
        }
        try record(gpu,"reactive-metal-predictions")
    }
    @Test func realMetalPotentialFeedsRingWorkflowAndResume() async throws {
        let spec = try specification(k: 40), potential = try await spec.make()
        let q = PrecisionSamplingFixtures.pairPositions(1.1), e = try await potential.checked(q)
        #expect(abs(e.energyKJPerMol-0.2) < 1e-5)
        #expect(potential.definition.atomIndices == [1,0] && potential.definition.coordinateEvaluation == .projectedFP32)
        let c = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 2,timeStepPS: 0.0005,integrationSteps: 2)
        let cp = try VivoRingPolymerCheckpoint(definition: spec.definition(),configuration: c,seed: 17,beadPositionsNM: [q,q])
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: PrecisionSamplingFixtures.id("precision-workflow"))
        let operation = try registry.definition("vivo.platform.ring-polymer-sample").operation
        let input = ["request":try VivoCanonicalJSON.encode(VivoRingPolymerWorkflowRequest(potential: spec,checkpoint: cp,sweeps: 3))]
        let out = try await operation.execute(.object([:]),input,.init())
        try operation.validateOutputs(.object([:]),input,out,.init())
        let run = try VivoPlatformOperations.input(VivoRingPolymerRun.self,"run",out)
        #expect(run.observations.count == 3 && run.end.sweep == 3)
        let resumed = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: run.end,sweeps: 2)
        try resumed.validate()
        #expect(resumed.end.sweep == 5)
        try record(run,"ring-native-workflow")
        var forged = out; forged["checkpoint"] = try VivoCanonicalJSON.encode(cp)
        #expect(throws: (any Error).self) { try operation.validateOutputs(.object([:]),input,forged,.init()) }
    }
    @Test func nativeLabelFitAndCorrectedSurrogateWorkflow() async throws {
        let a = try specification(k: 40), b = try specification(k: 10)
        var geometries: [VivoReactiveLabelGeometry] = []
        for i in 0...12 { geometries.append(.init(identifier: "train-\(i)",sourceGroup: "train-\(i)",positionsNM: PrecisionSamplingFixtures.pairPositions(0.8+Double(i)/30))) }
        for (i,r) in [0.85,1.05,1.15].enumerated() { geometries.append(.init(identifier: "held-\(i)",sourceGroup: "held",positionsNM: PrecisionSamplingFixtures.pairPositions(r))) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: PrecisionSamplingFixtures.id("precision-chain"))
        let labelOp = try registry.definition("vivo.platform.reactive-surrogate-label").operation
        let input = ["request":try VivoCanonicalJSON.encode(VivoReactiveLabelWorkflowRequest(authority: a,baseline: b,geometries: geometries))]
        let out = try await labelOp.execute(.object([:]),input,.init())
        try labelOp.validateOutputs(.object([:]),input,out,.init())
        let labels = try VivoPlatformOperations.input([VivoReactiveTrainingLabel].self,"labels",out)
        let cfg = VivoReactiveSurrogateConfiguration(features: [.init(atomA: 0,atomB: 1,scaleNM: 0.2)],maximumCenters: 10,committeeSize: 3,
            maximumHeldOutEnergyErrorKJPerMol: 0.02,maximumHeldOutForceErrorKJPerMolNM: 0.1)
        let training = VivoReactiveTrainingRequest(authority: try a.definition(),baseline: try b.definition(),labels: labels,heldOutGroups: ["held"],configuration: cfg)
        let fit = try registry.definition("vivo.platform.reactive-surrogate-train").operation
        let fitInput = ["request":try VivoCanonicalJSON.encode(training)]
        let fitted = try await fit.execute(.object([:]),fitInput,.init())
        try fit.validateOutputs(.object([:]),fitInput,fitted,.init())
        let model = try VivoPlatformOperations.input(VivoReactiveSurrogateModel.self,"model",fitted)
        #expect(model.payload.qualification.passed)
        let cp = try VivoReactiveSamplingCheckpoint(model: model,configuration: .init(temperatureK: 300,timeStepPS: 0.001,integrationSteps: 3),
            chainIdentifier: "native-production-epoch-1",seed: 4,positionsNM: PrecisionSamplingFixtures.pairPositions(1),backendProfile: "metal-fp32-rbf-analytic-derivative/v1")
        let sampling = VivoReactiveSamplingWorkflowRequest(authority: a,baseline: b,model: model,checkpoint: cp,sweeps: 4)
        let op = try registry.definition("vivo.platform.reactive-surrogate-sample").operation
        let runInput = ["request":try VivoCanonicalJSON.encode(sampling)]
        let runOut = try await op.execute(.object([:]),runInput,.init())
        try op.validateOutputs(.object([:]),runInput,runOut,.init())
        let run = try VivoPlatformOperations.input(VivoReactiveSamplingRun.self,"run",runOut)
        #expect(run.authorityLabels.count == 4 && run.authorityEvaluations == 5)
        #expect(run.baselineEvaluations > run.authorityEvaluations)
        try record(run,"reactive-native-corrected-sampling")
    }
    @Test func genuineEmbeddedElectronicAuthorityEvaluatesEveryBead() async throws {
        let h = try #require(VivoElement.from(symbol: "H")), he = try #require(VivoElement.from(symbol: "He"))
        let atoms = [VivoMolecularAtom(index: 0,name: "H1",element: h),.init(index: 1,name: "H2",element: h),
                     .init(index: 2,name: "MM1",element: he),.init(index: 3,name: "MM2",element: he)]
        let q: [VivoVector3D] = [.init(0.5,0.5,0.5),.init(0.574,0.5,0.5),.init(0.7,0.7,0.7),.init(1.1,0.8,1)]
        let document = try VivoMolecularStructureDocument(structure: .init(identifier: "nuclear-embedded-H2",
            atoms: atoms,bonds: [.init(atomA: 0,atomB: 1)],conformers: [.init(positionsNM: q)]))
        let charges = [0.0,0.0,0.2,-0.2]
        let system = VivoClassicalSystem(identifier: "nuclear-embedded-H2",structureFingerprint: document.structureFingerprint,
            particles: atoms.map { .init(index: $0.index,atomIndex: $0.index,typeIdentifier: $0.element.symbol,massDa: 2,
                chargeE: charges[Int($0.index)],sigmaNM: 0,epsilonKJPerMol: 0) },
            bonds: [.init(a: 0,b: 1,lengthNM: 0.074,forceConstant: 100)])
        let plan = try VivoQMMMHamiltonianPlan.compile(document: document,system: system,
            configuration: .init(region: .init(qmAtomIndices: [0,1],alphaElectrons: 1,betaElectrons: 1),boundary: .finiteCluster))
        let electronic = VivoQMMMDynamicsElectronicConfiguration(basis: .hydrogenSTO3G(nucleusIndices: [0,1]),
            budget: .init(maximumBytes: 16*1024*1024))
        let dynamics = VivoMDConfiguration(timeStepPS: 0.00001,cutoffNM: 2,neighborSkinNM: 0,electrostatics: .cutoff,
            ensemble: .nve,thermostat: .none,targetTemperatureK: nil,frictionPerPS: nil,neighborListEnabled: false)
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),positionsNM: q)
        let spec = VivoNuclearMetalSpecification(system: system,initialState: initial,dynamics: dynamics,
            electronic: .fixed(document: document,plan: plan,electronic: electronic))
        let potential = try await spec.make()
        let baseline = try await potential.checked(q)
        #expect(baseline.energyKJPerMol < -1000 && baseline.normalizedConvergenceResidual <= 1)
        let cfg = VivoRingPolymerConfiguration(temperatureK: 300,beadCount: 2,timeStepPS: 0.00001,integrationSteps: 1)
        let cp = try VivoRingPolymerCheckpoint(definition: potential.definition,configuration: cfg,seed: 31,beadPositionsNM: [q,q])
        let run = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: cp,sweeps: 2)
        try run.validate()
        #expect(run.forceEvaluations == 6 && run.observations.count == 2)
        #expect(run.observations.allSatisfy { $0.meanPotentialEnergyKJPerMol < -1000 })
        try record(run,"ring-native-QMMM-authority")
    }

}
