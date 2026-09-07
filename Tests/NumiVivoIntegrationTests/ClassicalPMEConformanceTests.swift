import Foundation
import Darwin
@preconcurrency import Metal
import Testing
import NumiVivoKit

/// Native force/energy conformance only. No integration, thermostat, barostat,
/// ensemble sampling, timing or performance qualification follows from this suite.
@Suite(.serialized) struct ClassicalPMEConformanceTests {
    // Preregistered before hardware measurements. Do not enlarge these limits
    // in response to a failed run; retain the observations and diagnose the owner.
    private struct Limits: Encodable {
        let fineEnergyAbsoluteKJPerMol = 0.05
        let fineForceNormalizedRMS = 1e-3
        let fineForceNormalizedMaximumComponent = 3e-3
        let forceRMSDenominatorFloorKJPerMolNM = 1.0
        let oracleEnergyAbsoluteKJPerMol = 1e-5
        let oracleForceMaximumComponentKJPerMolNM = 1e-4
        let aggregateEnergyFloorKJPerMol = 1e-3
        let aggregateForceNormalizedFloor = 1e-5
    }
    private let limits = Limits()
    private let cutoffNM = 0.75
    private let expandedOracleCutoffNM = 1.125
    private let pmeTolerance = 1e-7
    private let coulombPrefactorKJNMPerMol = 138.935456

    private struct Fixture: Encodable, Sendable {
        let name: String
        let chargesE: [Double]
        let geometry: VivoMDCandidateGeometry
        let primaryPair01Scale: Double?
    }
    private struct ForceError: Encodable {
        let referenceComponentRMSKJPerMolNM: Double
        let normalizationKJPerMolNM: Double
        let componentRMSAbsoluteKJPerMolNM: Double
        let normalizedRMS: Double
        let maximumComponentAbsoluteKJPerMolNM: Double
        let maximumComponentNormalized: Double
        let componentErrorsKJPerMolNM: [VivoVector3D]
    }
    private struct ConvertedOracle: Encodable {
        let energyKJPerMol: Double
        let forcesKJPerMolNM: [VivoVector3D]
        let realEnergyKJPerMol: Double
        let reciprocalEnergyKJPerMol: Double
        let selfEnergyKJPerMol: Double
        let backgroundEnergyKJPerMol: Double
        let exceptionEnergyKJPerMol: Double
        let reciprocalModes: Int
        let imagePairEvaluations: Int
    }
    private struct OracleDifference: Encodable {
        let energyAbsoluteKJPerMol: Double
        let force: ForceError
    }
    private struct OracleRecord: Encodable {
        let fixture: Fixture
        let referenceConfiguration: VivoPeriodicElectrostaticConfiguration
        let referenceBetaPerNM: Double
        let nativeFP32BetaPerNM: Double
        let bohrInNM: Double
        let hartreeInKJPerMol: Double
        let mdCoulombPrefactorKJNMPerMol: Double
        let atomicUnitPrefactorNormalization: Double
        let energyConversionKJPerMolPerHartree: Double
        let forceConversionKJPerMolNMPerHartreeBohr: Double
        let modes14ExpandedCutoff: ConvertedOracle
        let modes18NativeCutoff: ConvertedOracle
        let referenceModes18ExpandedCutoff: ConvertedOracle
        let reciprocalConvergence: OracleDifference
        let realCutoffConvergence: OracleDifference
        let analyticBackgroundEnergyKJPerMol: Double
        let analyticPrimaryExceptionEnergyKJPerMol: Double
        let minimumImagePairDistancesNM: [Double]
    }
    private struct ProbeRecord: Encodable {
        let fixture: Fixture
        let gridDimensions: [UInt32]
        let configuration: VivoMDConfiguration
        let deviceName: String
        let deviceRegistryID: UInt64
        let acceptedCheckpointBefore: VivoMDCheckpoint
        let acceptedCheckpointAfter: VivoMDCheckpoint
        let acceptedCheckpointUnchanged: Bool
        let evaluatedGeometryMatchesOracle: Bool
        let native: VivoMDHamiltonianEvaluation
        let reference: ConvertedOracle
        let energyErrorKJPerMol: Double
        let forceError: ForceError
    }
    private struct Aggregate: Encodable {
        let grid: UInt32
        let probeCount: Int
        let energyRMSAbsoluteKJPerMol: Double
        let forceNormalizedRMS: Double
    }

    private func cell(_ side: Double = 2) -> VivoPeriodicCell {
        .init(a: .init(side, 0, 0), b: .init(0, side, 0), c: .init(0, 0, side))
    }
    private func basePositions() -> [VivoVector3D] {
        // Odd multiples of 1/256 nm are exactly FP32, off mesh nodes at every
        // declared grid. Pair 0-1 crosses the boundary; pair 0-2 is ~0.7443 nm.
        [[43, 73, 107], [477, 101, 131], [233, 85, 115], [309, 351, 409]].map {
            .init(Double($0[0]) / 256, Double($0[1]) / 256, Double($0[2]) / 256)
        }
    }
    private func fixtures() throws -> [Fixture] {
        let base = basePositions(), neutral = [1.0, -1, 0.5, -0.5], charged = [1.0, -0.5, 0.5, -0.25]
        let displaced = [base[0], base[1] + .init(-4.0 / 256, 6.0 / 256, 2.0 / 256),
                         base[2] + .init(0, 4.0 / 256, -6.0 / 256),
                         base[3] + .init(-8.0 / 256, -2.0 / 256, 10.0 / 256)]
        func fixture(_ name: String, _ positions: [VivoVector3D], charges: [Double] = [1, -1, 0.5, -0.5],
                     scale: Double? = nil, box: VivoPeriodicCell? = nil) throws -> Fixture {
            .init(name: name, chargesE: charges,
                  geometry: try .init(particlePositionsNM: positions, periodicCell: box ?? cell()),
                  primaryPair01Scale: scale)
        }
        return try [
            fixture("neutral-base", base, charges: neutral),
            fixture("neutral-displaced", displaced),
            fixture("neutral-subcell-translation", base.map { $0 + .init(6.0 / 256, 10.0 / 256, 14.0 / 256) }),
            fixture("neutral-excluded-primary", base, scale: 0),
            fixture("neutral-half-scaled-primary", base, scale: 0.5),
            fixture("charged-background", base, charges: charged),
            fixture("charged-half-scaled-primary", base, charges: charged, scale: 0.5),
            fixture("neutral-lattice-translation", base.map { $0 + .init(2, -2, 4) }),
            fixture("neutral-affine-cell", base.map { $0 * (17.0 / 16) }, box: cell(2.125))
        ]
    }
    private func configuration(grid: UInt32) -> VivoMDConfiguration {
        .init(timeStepPS: 1e-5, cutoffNM: cutoffNM, neighborSkinNM: 0.125,
              electrostatics: .pme, relativeDielectric: 1, pmeTolerance: pmeTolerance,
              pmeGridSpacingNM: 2 / Double(grid), ensemble: .nve, thermostat: .none,
              targetTemperatureK: nil, frictionPerPS: nil, barostat: .none,
              neighborListEnabled: true, maximumNeighborsPerParticle: 3,
              pmeGridDimensions: [grid, grid, grid])
    }
    private func system(_ fixture: Fixture) throws -> VivoClassicalSystem {
        let identity = try VivoCanonicalJSON.fingerprint(Data("classical-pme-conformance-v1".utf8))
        let particles = fixture.chargesE.enumerated().map { index, charge in
            VivoClassicalParticle(index: UInt32(index), atomIndex: UInt32(index), typeIdentifier: "charge-only",
                massDa: 12, chargeE: charge, sigmaNM: 0, epsilonKJPerMol: 0)
        }
        let exceptions = fixture.primaryPair01Scale.map {
            [VivoNonbondedException(a: 0, b: 1, coulombScale: $0, lennardJonesScale: 0)]
        } ?? []
        return .init(identifier: fixture.name, structureFingerprint: identity,
                     particles: particles, nonbondedExceptions: exceptions)
    }
    private func components(_ vector: VivoVector3D) -> [Double] { [vector.x, vector.y, vector.z] }
    private func forceError(_ actual: [VivoVector3D], _ reference: [VivoVector3D]) throws -> ForceError {
        try #require(actual.count == reference.count && !reference.isEmpty)
        let errors = zip(actual, reference).map { $0 - $1 }
        let actualComponents = errors.flatMap(components), referenceComponents = reference.flatMap(components)
        let referenceRMS = sqrt(referenceComponents.reduce(0) { $0 + $1 * $1 } / Double(referenceComponents.count))
        let absoluteRMS = sqrt(actualComponents.reduce(0) { $0 + $1 * $1 } / Double(actualComponents.count))
        let normalization = max(referenceRMS, limits.forceRMSDenominatorFloorKJPerMolNM)
        let maximum = actualComponents.map(abs).max() ?? 0
        return .init(referenceComponentRMSKJPerMolNM: referenceRMS, normalizationKJPerMolNM: normalization,
            componentRMSAbsoluteKJPerMolNM: absoluteRMS, normalizedRMS: absoluteRMS / normalization,
            maximumComponentAbsoluteKJPerMolNM: maximum, maximumComponentNormalized: maximum / normalization,
            componentErrorsKJPerMolNM: errors)
    }
    private func oracle(_ fixture: Fixture) throws -> OracleRecord {
        let box = try #require(fixture.geometry.periodicCell)
        let cfg = configuration(grid: 64)
        let plan = try VivoPMEPlan.make(cell: box, cutoffNM: cutoffNM, tolerance: pmeTolerance,
            targetGridSpacingNM: cfg.resolvedPMEGridSpacingNM, fixedGridDimensions: cfg.pmeGridDimensions)
        let bohrNM = VivoAtomicUnits.bohrInNM, hartreeKJ = VivoAtomicUnits.hartreeInKJPerMol
        // Atomic electrostatics uses e, Bohr and Hartree. Normalize its Coulomb
        // constant explicitly to the public classical MD value (not a fitted
        // factor and not the kernel's rounded FP32 beta or prefactor).
        let normalization = coulombPrefactorKJNMPerMol / (hartreeKJ * bohrNM)
        let energyConversion = hartreeKJ * normalization
        let forceConversion = hartreeKJ / bohrNM * normalization
        let sources = zip(fixture.geometry.particlePositionsNM, fixture.chargesE).map {
            VivoCartesianMultipole(positionBohr: $0 / bohrNM, chargeE: $1)
        }
        let scales = fixture.primaryPair01Scale.map { [VivoMultipolePairScale(first: 0, second: 1, scale: $0)] } ?? []
        let totalCharge = fixture.chargesE.reduce(0, +)
        let referenceConfiguration = VivoPeriodicElectrostaticConfiguration(alphaPerBohr: plan.ewaldBetaPerNM * bohrNM,
            realCutoffBohr: expandedOracleCutoffNM / bohrNM, reciprocalHalfWidths: [18, 18, 18],
            chargeConvention: totalCharge == 0 ? .requireNeutral : .uniformNeutralizingBackground)
        func evaluate(modes: Int, realCutoffNM: Double) throws -> ConvertedOracle {
            try Task.checkCancellation()
            var configuration = referenceConfiguration
            configuration.reciprocalHalfWidths = [modes, modes, modes]
            configuration.realCutoffBohr = realCutoffNM / bohrNM
            // No mesh/FFT operator: this is the direct FP64 mode summation.
            // The cell argument deliberately remains in nm; only sources use Bohr.
            let result = try VivoPeriodicElectrostatics.evaluate(sources: sources, cell: box,
                configuration: configuration, primaryPairScales: scales, reciprocalOperator: nil)
            try Task.checkCancellation()
            return .init(energyKJPerMol: result.energyHartree * energyConversion,
                forcesKJPerMolNM: result.forcesHartreePerBohr.map { $0 * forceConversion },
                realEnergyKJPerMol: result.realEnergyHartree * energyConversion,
                reciprocalEnergyKJPerMol: result.reciprocalEnergyHartree * energyConversion,
                selfEnergyKJPerMol: result.selfEnergyHartree * energyConversion,
                backgroundEnergyKJPerMol: result.backgroundEnergyHartree * energyConversion,
                exceptionEnergyKJPerMol: result.exceptionEnergyHartree * energyConversion,
                reciprocalModes: result.reciprocalModeCount, imagePairEvaluations: result.imagePairEvaluations)
        }
        let lowModes = try evaluate(modes: 14, realCutoffNM: expandedOracleCutoffNM)
        let nativeCutoff = try evaluate(modes: 18, realCutoffNM: cutoffNM)
        let reference = try evaluate(modes: 18, realCutoffNM: expandedOracleCutoffNM)
        let reciprocalDifference = OracleDifference(energyAbsoluteKJPerMol: abs(lowModes.energyKJPerMol - reference.energyKJPerMol),
            force: try forceError(lowModes.forcesKJPerMolNM, reference.forcesKJPerMolNM))
        let cutoffDifference = OracleDifference(energyAbsoluteKJPerMol: abs(nativeCutoff.energyKJPerMol - reference.energyKJPerMol),
            force: try forceError(nativeCutoff.forcesKJPerMolNM, reference.forcesKJPerMolNM))
        let background = -Double.pi * coulombPrefactorKJNMPerMol * totalCharge * totalCharge
            / (2 * plan.ewaldBetaPerNM * plan.ewaldBetaPerNM * box.volumeNM3)
        var distances: [Double] = []
        for i in 0..<fixture.chargesE.count { for j in (i + 1)..<fixture.chargesE.count {
            distances.append(try box.minimumImage(fixture.geometry.particlePositionsNM[i] - fixture.geometry.particlePositionsNM[j]).norm)
        } }
        let exception = (fixture.primaryPair01Scale ?? 1) - 1
        let exceptionEnergy = exception * coulombPrefactorKJNMPerMol * fixture.chargesE[0] * fixture.chargesE[1] / distances[0]
        return .init(fixture: fixture, referenceConfiguration: referenceConfiguration,
            referenceBetaPerNM: plan.ewaldBetaPerNM, nativeFP32BetaPerNM: Double(Float(plan.ewaldBetaPerNM)),
            bohrInNM: bohrNM, hartreeInKJPerMol: hartreeKJ, mdCoulombPrefactorKJNMPerMol: coulombPrefactorKJNMPerMol,
            atomicUnitPrefactorNormalization: normalization, energyConversionKJPerMolPerHartree: energyConversion,
            forceConversionKJPerMolNMPerHartreeBohr: forceConversion, modes14ExpandedCutoff: lowModes,
            modes18NativeCutoff: nativeCutoff, referenceModes18ExpandedCutoff: reference,
            reciprocalConvergence: reciprocalDifference, realCutoffConvergence: cutoffDifference,
            analyticBackgroundEnergyKJPerMol: background, analyticPrimaryExceptionEnergyKJPerMol: exceptionEnergy,
            minimumImagePairDistancesNM: distances)
    }
    private func aggregate(_ probes: [ProbeRecord], grid: UInt32) -> Aggregate {
        let selected = probes.filter { $0.gridDimensions[0] == grid }
        return .init(grid: grid, probeCount: selected.count,
            energyRMSAbsoluteKJPerMol: sqrt(selected.reduce(0) { $0 + $1.energyErrorKJPerMol * $1.energyErrorKJPerMol } / Double(selected.count)),
            forceNormalizedRMS: sqrt(selected.reduce(0) { $0 + $1.forceError.normalizedRMS * $1.forceError.normalizedRMS } / Double(selected.count)))
    }

    @Test func plannerDeclaresTheImplementedOrderAndPreservesExplicitCoarseMeshes() throws {
        let coarse = try VivoPMEPlan.make(cell: cell(), cutoffNM: cutoffNM,
            tolerance: pmeTolerance, targetGridSpacingNM: 0.5, fixedGridDimensions: [4, 4, 4])
        #expect(coarse.interpolationOrder == 6)
        #expect(coarse.gridX == 4 && coarse.gridY == 4 && coarse.gridZ == 4 && coarse.gridPointCount == 64)
        // An explicitly requested prior interpolation algorithm cannot silently
        // become a different implementation. Coarse grid accuracy is separate.
        #expect(throws: VivoMDRuntimeError.self) {
            try VivoPMEPlan.make(cell: cell(), cutoffNM: cutoffNM,
                tolerance: pmeTolerance, targetGridSpacingNM: 0.5, interpolationOrder: 4)
        }
    }

    @Test func chargedClassicalPMEConvergesToIndependentDirectEwaldWithoutChangingAcceptedState() async throws {
        let evidence = try PMEEvidenceSink()
        try evidence.record(limits, name: "preregistered-limits")
        let cases = try fixtures()
        try #require(cases.count == 9)
        var references: [OracleRecord] = []
        for fixture in cases {
            let result = try oracle(fixture)
            references.append(result)
            try evidence.record(result, name: "oracle-\(fixture.name)")
        }
        // Validate the independent reference before running or judging mesh probes.
        for reference in references {
            for difference in [reference.reciprocalConvergence, reference.realCutoffConvergence] {
                try #require(difference.energyAbsoluteKJPerMol <= limits.oracleEnergyAbsoluteKJPerMol)
                try #require(difference.force.maximumComponentAbsoluteKJPerMolNM <= limits.oracleForceMaximumComponentKJPerMolNM)
            }
            try #require(abs(reference.referenceModes18ExpandedCutoff.backgroundEnergyKJPerMol - reference.analyticBackgroundEnergyKJPerMol) <= 1e-10)
            try #require(abs(reference.referenceModes18ExpandedCutoff.exceptionEnergyKJPerMol - reference.analyticPrimaryExceptionEnergyKJPerMol) <= 1e-9)
            if reference.fixture.primaryPair01Scale != nil {
                try #require(reference.minimumImagePairDistancesNM[0] < cutoffNM)
            }
        }
        let baseReference = try #require(references.first)
        try #require(baseReference.minimumImagePairDistancesNM[0] < 0.4)
        try #require(baseReference.minimumImagePairDistancesNM[1] > 0.74 && baseReference.minimumImagePairDistancesNM[1] < cutoffNM)
        try #require(baseReference.minimumImagePairDistancesNM.contains { $0 > cutoffNM })

        let device = try VivoMetalDeviceSelector.productionDevice()
        try #require(device.hasUnifiedMemory)
        var probes: [ProbeRecord] = []
        for grid in [UInt32(16), 32, 64] {
            for reference in references {
                let fixture = reference.fixture, model = try system(fixture), cfg = configuration(grid: grid)
                // Every candidate starts from the same accepted geometry
                // and 2nm cell. The affine-cell case exercises a real trial-cell
                // resource change with fixed mesh dimensions and no acceptance.
                let initial = VivoClassicalInitialState(systemFingerprint: try model.fingerprint(),
                    positionsNM: basePositions(), periodicCell: cell(), sourceTimePS: 0.125)
                let runtime = try await VivoMDMetalRuntime.make(system: model, initialState: initial, configuration: cfg,
                    initialVelocitiesNMPerPS: [.init(0.125, 0, 0), .init(0, -0.25, 0), .init(0, 0, 0.375), .init(-0.125, 0.25, -0.375)], device: device)
                let before = try await runtime.checkpoint()
                let evaluation = try await runtime.evaluateHamiltonian(at: fixture.geometry)
                let after = try await runtime.checkpoint()
                let exactGeometry = evaluation.evaluatedGeometry == fixture.geometry
                let oracle = reference.referenceModes18ExpandedCutoff
                let record = ProbeRecord(fixture: fixture, gridDimensions: [grid, grid, grid], configuration: cfg,
                    deviceName: device.name, deviceRegistryID: device.registryID,
                    acceptedCheckpointBefore: before, acceptedCheckpointAfter: after,
                    acceptedCheckpointUnchanged: before == after, evaluatedGeometryMatchesOracle: exactGeometry,
                    native: evaluation, reference: oracle, energyErrorKJPerMol: evaluation.energyKJPerMol - oracle.energyKJPerMol,
                    forceError: try forceError(evaluation.physicalParticleForcesKJPerMolNM, oracle.forcesKJPerMolNM))
                probes.append(record)
                // Persist full component observations before any native-error gate.
                try evidence.record(record, name: "probe-\(fixture.name)-grid-\(grid)")
                try #require(record.acceptedCheckpointUnchanged)
                // There are no virtual sites or constraints here, so the returned
                // physical coordinates must exactly match the oracle geometry.
                try #require(record.evaluatedGeometryMatchesOracle)
            }
        }
        try #require(probes.count == 27)
        let aggregates = [UInt32(16), 32, 64].map { aggregate(probes, grid: $0) }
        try evidence.record(aggregates, name: "grid-error-aggregates")
        for probe in probes where probe.gridDimensions[0] == 64 {
            #expect(abs(probe.energyErrorKJPerMol) <= limits.fineEnergyAbsoluteKJPerMol)
            #expect(probe.forceError.normalizedRMS <= limits.fineForceNormalizedRMS)
            #expect(probe.forceError.maximumComponentNormalized <= limits.fineForceNormalizedMaximumComponent)
        }
        let coarse = aggregates[0], fine = aggregates[2]
        #expect(fine.energyRMSAbsoluteKJPerMol < coarse.energyRMSAbsoluteKJPerMol ||
            max(fine.energyRMSAbsoluteKJPerMol, coarse.energyRMSAbsoluteKJPerMol) <= limits.aggregateEnergyFloorKJPerMol)
        #expect(fine.forceNormalizedRMS < coarse.forceNormalizedRMS ||
            max(fine.forceNormalizedRMS, coarse.forceNormalizedRMS) <= limits.aggregateForceNormalizedFloor)
        // Whole-lattice translation preserves the physical problem and mesh phase.
        // This differs from the subcell shift, whose mesh errors are retained above.
        for grid in [UInt32(16), 32, 64] {
            let base = try #require(probes.first { $0.fixture.name == "neutral-base" && $0.gridDimensions[0] == grid })
            let translated = try #require(probes.first { $0.fixture.name == "neutral-lattice-translation" && $0.gridDimensions[0] == grid })
            #expect(abs(base.native.energyKJPerMol - translated.native.energyKJPerMol) <= limits.aggregateEnergyFloorKJPerMol)
            let error = try forceError(translated.native.physicalParticleForcesKJPerMolNM, base.native.physicalParticleForcesKJPerMolNM)
            #expect(error.normalizedRMS <= limits.aggregateForceNormalizedFloor)
        }
    }
}

/// Raw observations are measurements, not success receipts. Each run owns a new
/// directory and never overwrites earlier evidence. No peak-RSS claim is made.
private struct PMEEvidenceSink {
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
        let directoryURL = parent.appendingPathComponent("classical-pme-conformance-\(UUID().uuidString)")
        guard Darwin.mkdir(directoryURL.path, 0o700) == 0 else {
            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: [NSFilePathErrorKey: directoryURL.path])
        }
        root = directoryURL
        print("NUMIVIVO_PME_EVIDENCE_ROOT=\(directoryURL.path)")
    }
    func record<T: Encodable>(_ value: T, name: String) throws {
        let bytes = try VivoCanonicalJSON.encode(value)
        if let root { try bytes.write(to: root.appendingPathComponent(name + ".json"), options: .withoutOverwriting) }
        print("NUMIVIVO_PME_OBSERVATION[\(name)]=\(String(decoding: bytes, as: UTF8.self))")
    }
}
