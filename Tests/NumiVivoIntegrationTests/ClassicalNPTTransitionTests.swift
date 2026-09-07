import Foundation
import Darwin
@preconcurrency import Metal
import Testing
import NumiVivoKit

/// Conditional NPT proposal-score and state-transition qualification only.
/// These gamma-zero, zero-force fixtures do not qualify an equilibrium volume
/// distribution, thermostat, random-number distribution or interacting NPT.
@Suite(.serialized) struct ClassicalNPTTransitionTests {
    // Fixed before native measurements. Preserve failed observations rather than
    // selecting a new seed or widening a limit after seeing the outcome.
    private struct Limits: Encodable {
        let proposalCount = 128
        let absoluteScoreError = 1e-9
        let relativeCellAndVolumeError = 1e-12
        let maximumPositionErrorNM = 1e-5
        let maximumCenterErrorNM = 5e-6
        let maximumInternalLengthErrorNM = 6.25e-7
        let maximumInternalVectorChangeNM = 2e-6
        let maximumVelocityNMPerPS = 5e-5
        let absolutePotentialEnergyKJPerMol = 1e-12
        let maximumKineticEnergyKJPerMol = 1e-6
    }
    private let limits = Limits()
    private let seed: UInt64 = 0x4e50545452414e53
    private let initialTimePS = 0.375
    private let timeStepPS = 0.125
    private let temperatureK = 300.0
    private let pressureBar = 150.0
    private let maximumLogVolumeStep = 0.25
    private let initialSideNM = 2.0
    private let dimerLengthNM = 0.125
    // SI definitions, converted independently of the native acceptance helper.
    private let boltzmannJPerK = 1.380649e-23
    private let avogadroPerMol = 6.02214076e23

    private struct Fixture: Encodable {
        let name: String
        let massesDa: [Double]
        let moleculeParticles: [[Int]]
        let initialPositionsNM: [VivoVector3D]
        let system: VivoClassicalSystem
    }
    private struct RunDefinition: Encodable {
        let fixture: Fixture
        let configuration: VivoMDConfiguration
        let limits: Limits
        let initialTimePS: Double
        let initialVelocitiesNMPerPS: [VivoVector3D]
        let boltzmannJPerK: Double
        let avogadroPerMol: Double
        let pressureVolumeKJPerMolPerBarNM3: Double
        let scope: String
    }
    private struct Observation: Encodable {
        let proposalIndex: Int
        let before: VivoMDCheckpoint
        let certificate: VivoMDStepCertificate
        let after: VivoMDCheckpoint
        let observables: VivoMDObservables
    }
    private struct AttemptFailure: Encodable {
        let proposalIndex: Int
        let before: VivoMDCheckpoint
        let certificate: VivoMDStepCertificate?
        let error: String
    }
    private struct Comparison: Encodable {
        let proposalIndex: Int
        let moleculeCount: Int
        let proposedVolumeNM3: Double
        let pressureWorkKJPerMol: Double
        let expectedLogAcceptanceRatio: Double
        let expectedLogAcceptanceRatioUsingAtomCount: Double
        let scoreAbsoluteError: Double
        let analyticalProposedCell: VivoPeriodicCell
        let analyticalProposedPositionsNM: [VivoVector3D]
        let analyticalProposedCentersNM: [VivoVector3D]
        let expectedRetainedCell: VivoPeriodicCell
        let cellAbsoluteErrorNM: Double
        let expectedRetainedPositionsNM: [VivoVector3D]
        let expectedRetainedCentersNM: [VivoVector3D]
        let measuredCentersNM: [VivoVector3D]
        let maximumPositionErrorNM: Double
        let maximumCenterErrorNM: Double
        let maximumInternalLengthErrorNM: Double
        let maximumInternalVectorChangeNM: Double
        let maximumVelocityNMPerPS: Double
        let maximumVelocityChangeNMPerPS: Double
        let positionsExactlyUnchanged: Bool
        let velocitiesExactlyUnchanged: Bool
    }
    private struct RunSummary: Encodable {
        let fixtureName: String
        let proposalCount: Int
        let acceptedProposalCount: Int
        let rejectedProposalCount: Int
        let finalCheckpoint: VivoMDCheckpoint
        let minimumVolumeNM3: Double
        let maximumVolumeNM3: Double
        let maximumScoreAbsoluteError: Double
        let maximumPositionErrorNM: Double
        let maximumCenterErrorNM: Double
        let maximumInternalLengthErrorNM: Double
        let maximumVelocityNMPerPS: Double
    }

    private func cell(side: Double) -> VivoPeriodicCell {
        .init(a: .init(side, 0, 0), b: .init(0, side, 0), c: .init(0, 0, side))
    }
    private func configuration() -> VivoMDConfiguration {
        .init(timeStepPS: timeStepPS, cutoffNM: 0.125, neighborSkinNM: 0.0625,
              electrostatics: .cutoff, ensemble: .npt, thermostat: .langevinMiddle,
              targetTemperatureK: temperatureK, frictionPerPS: 0,
              barostat: .monteCarloIsotropic, targetPressureBar: pressureBar,
              barostatInterval: 1, barostatMaximumLogVolumeStep: maximumLogVolumeStep,
              constraintTolerance: 1e-6, maximumConstraintIterations: 32,
              neighborListEnabled: false, maximumNeighborsPerParticle: 15, randomSeed: seed)
    }
    // Independent cubic-cell geometry; the fixture has no half-cell ties.
    private func image(_ displacement: VivoVector3D, side: Double) -> VivoVector3D {
        func reduce(_ x: Double) -> Double { x - side * floor(x / side + 0.5) }
        return .init(reduce(displacement.x), reduce(displacement.y), reduce(displacement.z))
    }
    private func wrap(_ position: VivoVector3D, side: Double) -> VivoVector3D {
        func reduce(_ x: Double) -> Double { x - side * floor(x / side) }
        return .init(reduce(position.x), reduce(position.y), reduce(position.z))
    }
    private func makeFixture(dimers: Bool) throws -> Fixture {
        var positions: [VivoVector3D] = [], masses: [Double] = []
        var molecules: [[Int]] = [], constraints: [VivoDistanceConstraint] = []
        for x in [0.03125, 1.03125] { for y in [0.5, 1.5] { for z in [0.5, 1.5] {
            let center = VivoVector3D(x, y, z), first = positions.count
            if dimers {
                // Masses 1 and 3 give these exact binary offsets from their COM.
                // Four molecules straddle the periodic boundary along x.
                positions.append(wrap(center + .init(-0.09375, 0, 0), side: initialSideNM))
                positions.append(wrap(center + .init(0.03125, 0, 0), side: initialSideNM))
                masses.append(contentsOf: [1, 3]); molecules.append([first, first + 1])
                constraints.append(.init(a: UInt32(first), b: UInt32(first + 1), distanceNM: dimerLengthNM))
            } else {
                positions.append(center); masses.append(4); molecules.append([first])
            }
        } } }
        let name = dimers ? "eight-unequal-mass-rigid-dimers" : "eight-monatomic-particles"
        let identity = try VivoCanonicalJSON.fingerprint(Data("classical-npt-transition-v1".utf8))
        let particles = masses.enumerated().map { index, mass in
            VivoClassicalParticle(index: UInt32(index), atomIndex: UInt32(index), typeIdentifier: "ideal-zero-force",
                massDa: mass, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
        }
        let system = VivoClassicalSystem(identifier: name, structureFingerprint: identity,
            particles: particles, constraints: constraints)
        return .init(name: name, massesDa: masses, moleculeParticles: molecules,
                     initialPositionsNM: positions, system: system)
    }
    private func unwrapped(_ positions: [VivoVector3D], members: [Int], side: Double) -> [VivoVector3D] {
        let anchor = positions[members[0]]
        return members.map { anchor + image(positions[$0] - anchor, side: side) }
    }
    private func center(_ coordinates: [VivoVector3D], members: [Int], masses: [Double]) -> VivoVector3D {
        var weighted = VivoVector3D.zero, mass = 0.0
        for (offset, particle) in members.enumerated() {
            weighted = weighted + coordinates[offset] * masses[particle]
            mass += masses[particle]
        }
        return weighted / mass
    }
    private func compare(_ row: Observation, fixture: Fixture) throws -> Comparison {
        let move = try #require(row.certificate.barostat)
        let oldCell = try #require(row.before.periodicCell), newCell = try #require(row.after.periodicCell)
        let oldSide = oldCell.a.x, oldVolume = oldCell.volumeNM3
        let delta = move.logVolumeDelta, scale = exp(delta / 3), proposedVolume = oldVolume * exp(delta)
        let pressureWork = pressureBar * 1e5 * (proposedVolume - oldVolume) * 1e-27 * avogadroPerMol / 1000
        let thermalEnergy = boltzmannJPerK * avogadroPerMol / 1000 * temperatureK
        let score = -pressureWork / thermalEnergy + Double(fixture.moleculeParticles.count + 1) * delta
        let atomScore = -pressureWork / thermalEnergy + Double(fixture.massesDa.count + 1) * delta
        let proposedCell = VivoPeriodicCell(a: oldCell.a * scale, b: oldCell.b * scale, c: oldCell.c * scale)
        let retainedScale = move.accepted ? scale : 1
        let expectedCell = VivoPeriodicCell(a: oldCell.a * retainedScale,
            b: oldCell.b * retainedScale, c: oldCell.c * retainedScale)
        let expectedSide = oldSide * retainedScale
        var expectedPositions = row.before.positionsNM
        var proposedPositions = row.before.positionsNM, proposedCenters: [VivoVector3D] = []
        var expectedCenters: [VivoVector3D] = [], measuredCenters: [VivoVector3D] = []
        var centerError = 0.0, lengthError = 0.0, internalChange = 0.0
        for members in fixture.moleculeParticles {
            let beforeCoordinates = unwrapped(row.before.positionsNM, members: members, side: oldSide)
            let oldCenter = center(beforeCoordinates, members: members, masses: fixture.massesDa)
            for (offset, particle) in members.enumerated() {
                proposedPositions[particle] = wrap(beforeCoordinates[offset] + oldCenter * (scale - 1), side: oldSide * scale)
                expectedPositions[particle] = wrap(beforeCoordinates[offset] + oldCenter * (retainedScale - 1), side: expectedSide)
            }
            proposedCenters.append(wrap(oldCenter * scale, side: oldSide * scale))
            let expectedCenter = wrap(oldCenter * retainedScale, side: expectedSide)
            let afterCoordinates = unwrapped(row.after.positionsNM, members: members, side: newCell.a.x)
            let measuredCenter = wrap(center(afterCoordinates, members: members, masses: fixture.massesDa), side: newCell.a.x)
            expectedCenters.append(expectedCenter); measuredCenters.append(measuredCenter)
            centerError = max(centerError, image(measuredCenter - expectedCenter, side: expectedSide).norm)
            if members.count == 2 {
                let beforeVector = beforeCoordinates[1] - beforeCoordinates[0]
                let afterVector = afterCoordinates[1] - afterCoordinates[0]
                lengthError = max(lengthError, abs(afterVector.norm - dimerLengthNM))
                internalChange = max(internalChange, (afterVector - beforeVector).norm)
            }
        }
        let positionError = zip(row.after.positionsNM, expectedPositions).map {
            image($0 - $1, side: expectedSide).norm
        }.max() ?? 0
        let velocity = row.after.velocitiesNMPerPS.map(\.norm).max() ?? 0
        let velocityChange = zip(row.after.velocitiesNMPerPS, row.before.velocitiesNMPerPS).map { ($0 - $1).norm }.max() ?? 0
        let cellError = max((newCell.a - expectedCell.a).norm,
            max((newCell.b - expectedCell.b).norm, (newCell.c - expectedCell.c).norm))
        return .init(proposalIndex: row.proposalIndex, moleculeCount: fixture.moleculeParticles.count,
            proposedVolumeNM3: proposedVolume, pressureWorkKJPerMol: pressureWork,
            expectedLogAcceptanceRatio: score, expectedLogAcceptanceRatioUsingAtomCount: atomScore,
            scoreAbsoluteError: abs(move.logAcceptanceRatio - score), analyticalProposedCell: proposedCell,
            analyticalProposedPositionsNM: proposedPositions, analyticalProposedCentersNM: proposedCenters,
            expectedRetainedCell: expectedCell,
            cellAbsoluteErrorNM: cellError, expectedRetainedPositionsNM: expectedPositions,
            expectedRetainedCentersNM: expectedCenters, measuredCentersNM: measuredCenters,
            maximumPositionErrorNM: positionError, maximumCenterErrorNM: centerError,
            maximumInternalLengthErrorNM: lengthError, maximumInternalVectorChangeNM: internalChange,
            maximumVelocityNMPerPS: velocity, maximumVelocityChangeNMPerPS: velocityChange,
            positionsExactlyUnchanged: row.before.positionsNM == row.after.positionsNM,
            velocitiesExactlyUnchanged: row.before.velocitiesNMPerPS == row.after.velocitiesNMPerPS)
    }

    @Test func idealMolecularVolumeProposalsHaveAnalyticalScoresAndPreserveAcceptedStateBoundaries() async throws {
        let evidence = try NPTTransitionEvidenceSink()
        let device = try VivoMetalDeviceSelector.productionDevice()
        try #require(device.hasUnifiedMemory)
        let cfg = configuration()
        var allRuns: [[Observation]] = []
        for dimers in [false, true] {
            let fixture = try makeFixture(dimers: dimers)
            let zeros = [VivoVector3D](repeating: .zero, count: fixture.massesDa.count)
            try evidence.record(RunDefinition(fixture: fixture, configuration: cfg, limits: limits,
                initialTimePS: initialTimePS, initialVelocitiesNMPerPS: zeros,
                boltzmannJPerK: boltzmannJPerK, avogadroPerMol: avogadroPerMol,
                pressureVolumeKJPerMolPerBarNM3: 1e5 * 1e-27 * avogadroPerMol / 1000,
                scope: "conditional ideal NPT proposal score and state transition; native pre-dynamics/post-step states and analytical proposed geometry; internal native proposal is not exposed; no equilibrium or RNG-distribution qualification"),
                name: "definition-\(fixture.name)")
            let initial = VivoClassicalInitialState(systemFingerprint: try fixture.system.fingerprint(),
                positionsNM: fixture.initialPositionsNM, periodicCell: cell(side: initialSideNM), sourceTimePS: initialTimePS)
            let runtime = try await VivoMDMetalRuntime.make(system: fixture.system, initialState: initial,
                configuration: cfg, initialVelocitiesNMPerPS: zeros, device: device)
            var rows: [Observation] = []
            for index in 0..<limits.proposalCount {
                try Task.checkCancellation()
                let before = try await runtime.checkpoint()
                var certificate: VivoMDStepCertificate?
                do {
                    let step = try await runtime.step(); certificate = step
                    let after = try await runtime.checkpoint(), observables = try await runtime.observables()
                    let row = Observation(proposalIndex: index, before: before, certificate: step,
                        after: after, observables: observables)
                    // Retain complete certificates and physical observations before
                    // applying any score, geometry or acceptance assertion.
                    try evidence.record(row, name: "observation-\(fixture.name)-\(index)")
                    rows.append(row)
                } catch {
                    try evidence.record(AttemptFailure(proposalIndex: index, before: before,
                        certificate: certificate, error: String(describing: error)), name: "failure-\(fixture.name)-\(index)")
                    throw error
                }
            }
            allRuns.append(rows)
            let comparisons = try rows.map { try compare($0, fixture: fixture) }
            try evidence.record(comparisons, name: "analytical-comparisons-\(fixture.name)")
            let accepted = rows.filter { $0.certificate.barostat?.accepted == true }.count
            let rejected = rows.filter { $0.certificate.barostat?.accepted == false }.count
            let final = try #require(rows.last?.after)
            let volumes = rows.compactMap { $0.after.periodicCell?.volumeNM3 }
            try evidence.record(RunSummary(fixtureName: fixture.name, proposalCount: rows.count,
                acceptedProposalCount: accepted, rejectedProposalCount: rejected, finalCheckpoint: final,
                minimumVolumeNM3: volumes.min() ?? 0, maximumVolumeNM3: volumes.max() ?? 0,
                maximumScoreAbsoluteError: comparisons.map(\.scoreAbsoluteError).max() ?? 0,
                maximumPositionErrorNM: comparisons.map(\.maximumPositionErrorNM).max() ?? 0,
                maximumCenterErrorNM: comparisons.map(\.maximumCenterErrorNM).max() ?? 0,
                maximumInternalLengthErrorNM: comparisons.map(\.maximumInternalLengthErrorNM).max() ?? 0,
                maximumVelocityNMPerPS: comparisons.map(\.maximumVelocityNMPerPS).max() ?? 0), name: "summary-\(fixture.name)")
            #expect(rows.count == limits.proposalCount)
            #expect(accepted > 0 && rejected > 0)
            for (row, comparison) in zip(rows, comparisons) {
                let move = try #require(row.certificate.barostat)
                let beforeCell = try #require(row.before.periodicCell), afterCell = try #require(row.after.periodicCell)
                #expect(row.certificate.committed && row.certificate.statusFlags == 0)
                #expect(row.certificate.firstViolationParticle == nil && row.certificate.violationCount == 0)
                #expect(row.after.systemFingerprint == row.before.systemFingerprint &&
                    row.after.configurationFingerprint == row.before.configurationFingerprint)
                #expect(row.certificate.systemFingerprint == row.before.systemFingerprint &&
                    row.certificate.configurationFingerprint == row.before.configurationFingerprint)
                #expect(move.attempted && move.rejectionReason == nil)
                #expect(abs(move.logVolumeDelta) <= maximumLogVolumeStep)
                #expect(row.before.acceptedStep == UInt64(row.proposalIndex))
                #expect(row.after.acceptedStep == row.before.acceptedStep + 1)
                #expect(row.certificate.stepIndex == row.before.acceptedStep)
                #expect(move.stepIndex == row.after.acceptedStep)
                // Binary timestep/start time make exact clock checks intentional.
                #expect(row.before.timePS == initialTimePS + Double(row.proposalIndex) * timeStepPS)
                #expect(row.after.timePS == row.before.timePS + timeStepPS)
                #expect(row.certificate.timeBeforePS == row.before.timePS && row.certificate.timeAfterPS == row.after.timePS)
                #expect(row.observables.stepIndex == row.after.acceptedStep && row.observables.timePS == row.after.timePS)
                #expect(move.volumeBeforeNM3 == beforeCell.volumeNM3 && move.volumeAfterNM3 == afterCell.volumeNM3)
                #expect(comparison.scoreAbsoluteError <= limits.absoluteScoreError)
                #expect(comparison.cellAbsoluteErrorNM <= limits.relativeCellAndVolumeError * max(1, comparison.expectedRetainedCell.a.norm))
                if move.accepted {
                    #expect(abs(afterCell.volumeNM3 - comparison.proposedVolumeNM3) <=
                        limits.relativeCellAndVolumeError * max(1, comparison.proposedVolumeNM3))
                } else {
                    #expect(afterCell == beforeCell)
                }
                // Favorable proposals have Metropolis probability one. Other
                // decisions are observed; RNG-distribution quality is separate.
                if comparison.expectedLogAcceptanceRatio >= limits.absoluteScoreError { #expect(move.accepted) }
                // Even at gamma=0, native FP32 wrapping and constraint projection
                // still execute before the pressure proposal. Exact bit equality
                // is recorded above; the full-step physical test allows this floor.
                #expect(comparison.maximumPositionErrorNM <= limits.maximumPositionErrorNM)
                #expect(comparison.maximumCenterErrorNM <= limits.maximumCenterErrorNM)
                #expect(comparison.maximumInternalLengthErrorNM <= limits.maximumInternalLengthErrorNM)
                #expect(comparison.maximumInternalVectorChangeNM <= limits.maximumInternalVectorChangeNM)
                #expect(comparison.maximumVelocityNMPerPS <= limits.maximumVelocityNMPerPS)
                #expect(comparison.maximumVelocityChangeNMPerPS <= limits.maximumVelocityNMPerPS)
                #expect(abs(row.observables.potentialEnergyKJPerMol) <= limits.absolutePotentialEnergyKJPerMol)
                #expect(row.observables.kineticEnergyKJPerMol <= limits.maximumKineticEnergyKJPerMol)
            }
        }
        // Both systems contain eight independent molecular centers and have the
        // same zero potential. A molecular barostat must therefore produce the
        // same scalar volume/acceptance trace for the same seed, despite 8 vs 16 atoms.
        try #require(allRuns.count == 2 && allRuns.allSatisfy { $0.count == limits.proposalCount })
        for (atoms, dimers) in zip(allRuns[0], allRuns[1]) {
            let a = try #require(atoms.certificate.barostat), d = try #require(dimers.certificate.barostat)
            #expect(a == d)
            #expect(atoms.after.periodicCell == dimers.after.periodicCell)
        }
    }
}

private struct NPTTransitionEvidenceSink {
    let root: URL?

    init() throws {
        guard let directory = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else {
            root = nil; return
        }
        guard !directory.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoArtifactValidationError.invalid("NUMIVIVO_TEST_ARTIFACTS must name a nonempty directory")
        }
        let parent = URL(fileURLWithPath: directory, isDirectory: true).standardizedFileURL
        try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: true)
        let directoryURL = parent.appendingPathComponent("classical-npt-transition-\(UUID().uuidString)")
        guard Darwin.mkdir(directoryURL.path, 0o700) == 0 else {
            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: [NSFilePathErrorKey: directoryURL.path])
        }
        root = directoryURL
        print("NUMIVIVO_NPT_EVIDENCE_ROOT=\(directoryURL.path)")
    }
    func record<T: Encodable>(_ value: T, name: String) throws {
        let bytes = try VivoCanonicalJSON.encode(value)
        if let root { try bytes.write(to: root.appendingPathComponent(name + ".json"), options: .withoutOverwriting) }
        print("NUMIVIVO_NPT_OBSERVATION[\(name)]=\(String(decoding: bytes, as: UTF8.self))")
    }
}
