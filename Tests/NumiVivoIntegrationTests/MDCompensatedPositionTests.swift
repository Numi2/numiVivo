import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MDCompensatedPositionTests {
    private func model(_ count:Int=3)throws->VivoClassicalSystem {
        .init(identifier:"compensated-coordinate-tests",structureFingerprint:try VivoCanonicalJSON.fingerprint(Data("compensated-tests".utf8)),
            particles:(0..<count).map { .init(index:UInt32($0),atomIndex:UInt32($0),typeIdentifier:"C",massDa:12,chargeE:0,sigmaNM:0,epsilonKJPerMol:0) })
    }
    private func config(_ iterations:UInt32=128)->VivoMDConfiguration {
        .init(timeStepPS:0.001,cutoffNM:1,electrostatics:.cutoff,ensemble:.nve,thermostat:.none,
            targetTemperatureK:nil,frictionPerPS:nil,maximumConstraintIterations:iterations,neighborListEnabled:false,positionPrecision:.compensated)
    }
    @Test func tinyMovementsAndExactWordsSurviveRestart() async throws {
        let model=try model(1),configuration=config()
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:[.init(1000,0,0)])
        let runtime=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:configuration,initialVelocitiesNMPerPS:[.init(1e-6,0,0)])
        for _ in 0..<16 { #expect(try await runtime.step().committed) }
        let checkpoint=try await runtime.checkpoint()
        #expect(abs(checkpoint.positionsNM[0].x-1000-16e-9)<1e-12)
        #expect(checkpoint.positionCorrectionsNM![0].x != 0)
        let restored=try await VivoMDMetalRuntime.restore(system:model,configuration:configuration,checkpoint:checkpoint)
        #expect(try await restored.checkpoint()==checkpoint)
        for _ in 0..<16 { #expect(try await runtime.step().committed);#expect(try await restored.step().committed) }
        #expect(try await runtime.checkpoint()==restored.checkpoint())
        // A low word too small for the Double coordinate view must still survive.
        var tiny=checkpoint
        tiny.positionHighNM=[.init(1000,0,0)];tiny.positionCorrectionsNM=[.init(Double(Float(1e-20)),0,0)]
        tiny.positionsNM=[.init(1000,0,0)]
        let exact=try await VivoMDMetalRuntime.restore(system:model,configuration:configuration,checkpoint:tiny)
        #expect(try await exact.checkpoint()==tiny)
        var missing=checkpoint;missing.positionCorrectionsNM=nil
        #expect(throws:Error.self) { try missing.validate(particleCount:1) }
        var noncanonical=tiny;noncanonical.positionHighNM=[.init(999,0,0)];noncanonical.positionCorrectionsNM=[.init(1,0,0)]
        #expect(throws:Error.self) { try noncanonical.validate(particleCount:1) }
        var wrong=configuration;wrong.positionPrecision = .fp32
        await #expect(throws:Error.self) { try await VivoMDMetalRuntime.restore(system:model,configuration:wrong,checkpoint:checkpoint) }
    }
    @Test(arguments:[false,true]) func strictConstraintsAtLargeCoordinatesAndPeriodicBoundary(_ periodic:Bool) async throws {
        var model=try model()
        model.constraints=[.init(a:0,b:1,distanceNM:0.09572),.init(a:0,b:2,distanceNM:0.09572),.init(a:1,b:2,distanceNM:0.15139)]
        let positions:[VivoVector3D]=periodic ? [.init(0.01,4,4),.init(7.91,4,4),.init(0.035,4.096,4)]
            : [.init(1000,1000,1000),.init(999.9,1000,1000),.init(1000.025,1000.096,1000)]
        let cell:VivoPeriodicCell?=periodic ? .init(a:.init(8,0,0),b:.init(0,8,0),c:.init(0,0,8)):nil
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:positions,periodicCell:cell)
        let runtime=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:config())
        _ = try await runtime.projectConstraints()
        _ = try await runtime.thermalize(temperatureK:300,seed:77)
        for _ in 0..<16 { #expect(try await runtime.step().committed) }
        let state=try await runtime.checkpoint()
        for q in model.constraints {
            var d=state.positionsNM[Int(q.a)]-state.positionsNM[Int(q.b)]
            if periodic { d = .init(d.x-8*(d.x/8).rounded(),d.y-8*(d.y/8).rounded(),d.z-8*(d.z/8).rounded()) }
            #expect(abs(d.norm-q.distanceNM)/q.distanceNM<1e-6)
        }
        let restarted=try await VivoMDMetalRuntime.restore(system:model,configuration:config(),checkpoint:state)
        #expect(try await restarted.checkpoint()==state)
        #expect(try await runtime.step().committed);#expect(try await restarted.step().committed)
        #expect(try await runtime.checkpoint()==restarted.checkpoint())
    }
    @Test func failedProjectionAndCancelledStepRetainBothWords() async throws {
        var model=try model()
        model.constraints=[.init(a:0,b:1,distanceNM:0.125),.init(a:1,b:2,distanceNM:0.125),.init(a:0,b:2,distanceNM:0.1875)]
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:[.init(1000,0,0),.init(1000.15625,0,0),.init(1000.15625,0.15625,0)])
        let runtime=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:config(1))
        let before=try await runtime.checkpoint()
        await #expect(throws:Error.self) { try await runtime.projectConstraints() }
        #expect(try await runtime.checkpoint()==before)
        let cancelled=Task { withUnsafeCurrentTask { $0?.cancel() };return try await runtime.step() }
        await #expect(throws:Error.self) { try await cancelled.value }
        #expect(try await runtime.checkpoint()==before)
        await #expect(throws:Error.self) { try await runtime.minimize() }
        #expect(try await runtime.checkpoint()==before)
    }
    @Test func unsupportedArchiveAndPrecisionChangesFailBeforeExecution() throws {
        let model=try model(),configuration=config()
        var npt=configuration;npt.ensemble = .npt;npt.barostat = .monteCarloIsotropic;npt.targetPressureBar=1
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:[.zero,.init(0.2,0,0),.init(0,0.2,0)])
        #expect(!VivoMDExecutionPreflight.blockers(system:model,initial:initial,configuration:npt).isEmpty)
        let stage=VivoMDProtocolStage(identifier:"sample",kind:.dynamics,configuration:configuration,steps:1,sampleEvery:1)
        #expect(throws:Error.self) { try stage.validate() }
    }
    @Test func benchmarkObservationSeriesIsBoundedAndUsesAcceptedSteps() async throws {
        let model=try model(1)
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:[.init(1,0,0)],periodicCell:nil)
        let reference=VivoMDBenchmarkReference(identifier:"zero-force",geometry:geometry,energyKJPerMol:0,forcesKJPerMolNM:[.zero])
        var request=VivoMDBenchmarkRequest(identifier:"observations",system:model,configuration:config(),references:[reference],
            referenceProvenance:["reference":"one free classical atom"],dynamicsSteps:8,dynamicsObserveEvery:2)
        let report=try await VivoMDBenchmark.run(request)
        #expect(report.outcome == .passed)
        #expect(report.dynamics?.observations?.map(\.stepIndex)==[0,2,4,6,8])
        request.dynamicsObserveEvery=0
        #expect(throws:Error.self) { try request.validate() }
        request.dynamicsObserveEvery=1;request.dynamicsSteps=1001
        #expect(throws:Error.self) { try request.validate() }
    }
    @Test(arguments:[VivoMDPositionPrecision.fp32,.compensated])
    func freeRigidRotorRetainsEnergyAndAngularMomentum(_ precision:VivoMDPositionPrecision) async throws {
        var model=try model(2),configuration=config(32)
        model.constraints=[.init(a:0,b:1,distanceNM:0.125)]
        configuration.positionPrecision=precision
        let initial=VivoClassicalInitialState(systemFingerprint:try model.fingerprint(),positionsNM:[.init(-0.0625,0,0),.init(0.0625,0,0)])
        let runtime=try await VivoMDMetalRuntime.make(system:model,initialState:initial,configuration:configuration,
            initialVelocitiesNMPerPS:[.init(0,1,0),.init(0,-1,0)])
        let energy=try await runtime.observables().totalEnergyKJPerMol
        var maximumError=0.0
        for i in 0..<300 {
            #expect(try await runtime.step().committed)
            if i%10==9 { maximumError=max(maximumError,abs(try await runtime.observables().totalEnergyKJPerMol/energy-1)) }
        }
        let state=try await runtime.checkpoint()
        let angular=zip(state.positionsNM,state.velocitiesNMPerPS).reduce(0.0) { $0+12*($1.0.x*$1.1.y-$1.0.y*$1.1.x) }
        let limit=precision == .compensated ? 2e-4:5e-3
        #expect(maximumError<limit)
        #expect(abs(angular/(-1.5)-1)<limit)
    }
}
