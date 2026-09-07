import Foundation
import Testing
@testable import NumiVivoKit

/// Host resource admission only: no Metal runtime or grid buffer is allocated.
struct NuclearMetalResourceTests {
    private func specification(grid: [UInt32]?,spacing: Double,maximumMeshPoints: Int) throws -> VivoNuclearMetalSpecification {
        let system = VivoClassicalSystem(identifier: "nuclear-mesh-reservation",
            structureFingerprint: try VivoCanonicalJSON.fingerprint(Data("nuclear-mesh-reservation".utf8)),particles: [
                .init(index: 0,atomIndex: 0,typeIdentifier: "X",massDa: 12,chargeE: 0,sigmaNM: 0,epsilonKJPerMol: 0)])
        let cell = VivoPeriodicCell(a: .init(2,0,0),b: .init(0,2,0),c: .init(0,0,2))
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),positionsNM: [.init(0.5,0.5,0.5)],periodicCell: cell)
        let dynamics = VivoMDConfiguration(timeStepPS: 0.001,cutoffNM: 0.75,neighborSkinNM: 0.125,
            electrostatics: .pme,pmeGridSpacingNM: spacing,ensemble: .nve,thermostat: .none,
            targetTemperatureK: nil,frictionPerPS: nil,pmeGridDimensions: grid)
        return .init(system: system,initialState: initial,dynamics: dynamics,maximumMeshPoints: maximumMeshPoints)
    }

    @Test func largerExplicitMeshControlsBothPointAndByteLimits() throws {
        let grid: [UInt32] = [128,64,32], points = 128*64*32
        // Spacing alone chooses 4^3. It must not mask the explicitly requested
        // 128*64*32 mesh during either point-count or memory admission.
        let implicit = try specification(grid: nil,spacing: 0.5,maximumMeshPoints: 64)
        let limitedBytes = VivoChemistryBudget(maximumBytes: 32*1024*1024)
        _ = try implicit.reservationBytes(budget: limitedBytes)
        let tooManyPoints = try specification(grid: grid,spacing: 0.5,maximumMeshPoints: points-1)
        #expect(throws: (any Error).self) { _ = try tooManyPoints.reservationBytes(budget: .init()) }
        let explicit = try specification(grid: grid,spacing: 0.5,maximumMeshPoints: points)
        #expect(throws: (any Error).self) { _ = try explicit.reservationBytes(budget: limitedBytes) }
        let reserved = try explicit.reservationBytes(budget: .init(maximumBytes: 128*1024*1024))
        // Even the two actual complex FFT buffers alone exceed the old small
        // reservation; the admission API deliberately reserves more scratch.
        #expect(reserved > points*2*MemoryLayout<SIMD2<Float>>.stride)
        #expect(reserved > limitedBytes.maximumBytes)
        #expect(try explicit.reservationBytes(budget: .init(maximumBytes: reserved)) == reserved)
        #expect(throws: (any Error).self) {
            _ = try explicit.reservationBytes(budget: .init(maximumBytes: reserved-1))
        }
    }

    @Test func smallerExplicitMeshDoesNotReserveIgnoredSpacingDerivedAxes() throws {
        let grid: [UInt32] = [4,4,4], budget = VivoChemistryBudget(maximumBytes: 1024*1024)
        let coarse = try specification(grid: grid,spacing: 0.5,maximumMeshPoints: 64)
        let expected = try coarse.reservationBytes(budget: budget)
        // The finer spacings exceed the existing derived-axis guard; none
        // changes the explicit four-point mesh in PMEEngine.
        for spacing in [0x1p-6,0x1p-20,1e-30] {
            let fixed = try specification(grid: grid,spacing: spacing,maximumMeshPoints: 64)
            #expect(try fixed.reservationBytes(budget: budget) == expected)
            let implicit = try specification(grid: nil,spacing: spacing,maximumMeshPoints: 64)
            #expect(throws: (any Error).self) { _ = try implicit.reservationBytes(budget: budget) }
        }
        let belowSmallCap = try specification(grid: grid,spacing: 0.5,maximumMeshPoints: 63)
        #expect(throws: (any Error).self) { _ = try belowSmallCap.reservationBytes(budget: budget) }
    }

    private func validateWorkflow(_ specification: VivoNuclearMetalSpecification,maximumMeshPoints: Int,
                                  budget: VivoChemistryBudget,resume: Bool) throws {
        let identifier = resume ? "vivo.platform.md-continue" : "vivo.platform.md-start"
        let definitions = VivoPlatformMDOperations.definitions(implementationFingerprint:
            try VivoCanonicalJSON.fingerprint(Data("host-mesh-reservation-validation".utf8)))
        let definition = try #require(definitions.first { $0.operation.identifier == identifier })
        let configuration = VivoWorkflowMDConfiguration(dynamics: specification.dynamics,steps: 1,
            maximumMeshPoints: maximumMeshPoints)
        let cfg = try VivoJSONValue.decode(data: VivoCanonicalJSON.encode(configuration))
        let systemID = try specification.system.fingerprint(), configurationID = try specification.dynamics.fingerprint()
        let start = VivoMDCheckpoint(systemFingerprint: systemID,configurationFingerprint: configurationID,
            acceptedStep: 0,timePS: 0,positionsNM: specification.initialState.positionsNM,
            velocitiesNMPerPS: [.zero],periodicCell: specification.initialState.periodicCell)
        var end = start; end.acceptedStep = 1; end.timePS = specification.dynamics.timeStepPS
        let observables = VivoMDObservables(systemFingerprint: systemID,configurationFingerprint: configurationID,
            stepIndex: end.acceptedStep,timePS: end.timePS,potentialEnergyKJPerMol: 0,kineticEnergyKJPerMol: 0,
            temperatureK: nil,degreesOfFreedom: 3)
        // Fully synthetic cached-output fixture. Calling only validateOutputs
        // exercises the same prepare/admission owner as execute, without a GPU.
        // The two callers must agree on the one-atom reservation.
        let reservation = try specification.reservationBytes(budget: .init(maximumBytes: 128*1024*1024))
        let result = VivoWorkflowMDStageResult(configuration: configuration,start: start,end: end,
            observables: observables,deviceName: "synthetic host validation fixture",deviceRegistryID: 0,
            allocationReservationBytes: reservation,validationScope: definition.validationScope)
        var input = ["system":try VivoCanonicalJSON.encode(specification.system)]
        if resume { input["checkpoint"] = try VivoCanonicalJSON.encode(start) }
        else { input["state"] = try VivoCanonicalJSON.encode(specification.initialState) }
        let output = ["checkpoint":try VivoCanonicalJSON.encode(end),"result":try VivoCanonicalJSON.encode(result)]
        try definition.operation.validateOutputs(cfg,input,output,budget)
    }

    @Test func workflowStartAndResumeAdmitTheSameExplicitMeshAsNuclearSampling() throws {
        let points = 128*64*32
        let large = try specification(grid: [128,64,32],spacing: 0.5,maximumMeshPoints: points)
        let tiny = try specification(grid: [4,4,4],spacing: 1e-30,maximumMeshPoints: 64)
        for resume in [false,true] {
            try validateWorkflow(large,maximumMeshPoints: points,budget: .init(maximumBytes: 128*1024*1024),resume: resume)
            #expect(throws: (any Error).self) {
                try validateWorkflow(large,maximumMeshPoints: points-1,budget: .init(maximumBytes: 128*1024*1024),resume: resume)
            }
            #expect(throws: (any Error).self) {
                try validateWorkflow(large,maximumMeshPoints: points,budget: .init(maximumBytes: 32*1024*1024),resume: resume)
            }
            try validateWorkflow(tiny,maximumMeshPoints: 64,budget: .init(maximumBytes: 1024*1024),resume: resume)
        }
    }

    @Test func commandPlanningHonorsFixedGridWithExtremelyFineValidSpacing() throws {
        let coarse = try specification(grid: [4,4,4],spacing: 0.5,maximumMeshPoints: 64)
        let fine = try specification(grid: [4,4,4],spacing: 1e-30,maximumMeshPoints: 64)
        try fine.dynamics.validate() // 1e-30 is positive and representable in FP32.
        let packed = try VivoMDSystemPacker.pack(fine.system)
        let reference = try VivoMDMetalABI.command(packed: packed,configuration: coarse.dynamics,
            cell: coarse.initialState.periodicCell,stepIndex: 0)
        let command = try VivoMDMetalABI.command(packed: packed,configuration: fine.dynamics,
            cell: fine.initialState.periodicCell,stepIndex: 0)
        #expect(command.electrostatics == 2 && command.periodic == 1)
        #expect(command.reactionFieldK.isFinite && command.reactionFieldK > 0)
        #expect(command.reactionFieldK == reference.reactionFieldK)
        #expect(command.reactionFieldC == reference.reactionFieldC)
    }

    @Test func directPlannerRejectsUnrepresentableDerivedAxesWithoutTrapping() throws {
        let cell = VivoPeriodicCell(a: .init(2,0,0),b: .init(0,2,0),c: .init(0,0,2))
        func plan(_ cell: VivoPeriodicCell,spacing: Double) throws -> VivoPMEPlan {
            try .make(cell: cell,cutoffNM: 0.75,tolerance: 1e-5,targetGridSpacingNM: spacing)
        }
        // Finite >Int range and infinite quotients previously reached Int(...).
        for spacing in [1e-30,Double.leastNonzeroMagnitude] {
            #expect(throws: (any Error).self) { _ = try plan(cell,spacing: spacing) }
        }
        let boundary = VivoPeriodicCell(a: .init(0x1p31,0,0),b: .init(0,4,0),c: .init(0,0,4))
        let largestAxis = try plan(boundary,spacing: 1)
        #expect(largestAxis.gridX == UInt32(1)<<31)
        #expect(largestAxis.gridPointCount == UInt64(1)<<35)
        // The planner's UInt64 count is separate from the engine's tighter ABI.
        // This finite UInt32 request needs the unrepresentable next power 2^32.
        let roundingOverflow = VivoPeriodicCell(a: .init(Double(UInt32.max),0,0),b: .init(0,4,0),c: .init(0,0,4))
        #expect(throws: (any Error).self) { _ = try plan(roundingOverflow,spacing: 1) }
        let countOverflow = VivoPeriodicCell(a: .init(0x1p22,0,0),b: .init(0,0x1p22,0),c: .init(0,0,0x1p22))
        #expect(throws: (any Error).self) { _ = try plan(countOverflow,spacing: 1) }
        let fixed = try VivoPMEPlan.make(cell: cell,cutoffNM: 0.75,tolerance: 1e-5,
            targetGridSpacingNM: 1e-30,fixedGridDimensions: [4,4,4])
        #expect(fixed.gridPointCount == 64)
    }
}
