import Foundation
import Darwin
@preconcurrency import Metal
import Testing
import NumiVivoKit

/// Two distinct qualifications: actual force-free NVT time evolution, and
/// Gaussian initialization at fixed constrained geometries. No interacting or
/// constrained-trajectory equilibrium claim follows from these tests.
@Suite(.serialized) struct ClassicalThermalConformanceTests {
    private let kB = 0.00831446261815324 // kJ mol^-1 K^-1; mass Da, velocity nm/ps.
    private let temperature = 300.0
    private struct Policy: Encodable {
        let familyFailureProbability = 1e-4
        let comparisonBudget = 512
        // OU: 3*(3 means+6 second moments)=27. Dimers:27+4*4=43.
        // Triatomics:54+4*4=70. Each chi-square mean interval is one
        // two-sided event, plus three CDF events. Deterministic ownership and
        // geometry assertions are not statistical comparisons. Independence
        // across horizons/families is unnecessary for this union bound.
        let declaredStatisticalComparisons = 140
        let relativeFP32MomentAllowance = 1e-5
        let relativeDistanceLimit = 1e-6
        let tangentSpeedLimitNMPerPS = 2e-5
        let maximumConstraintIterations: UInt32 = 32
        var logTail: Double { log(2 * Double(comparisonBudget) / familyFailureProbability) }
    }
    private let policy = Policy()
    private struct MomentReport: Encodable {
        let sampleCount: Int
        let expectedMean: [Double]
        let expectedCovariance: [[Double]]
        let observedMean: [Double]
        /// Centered about the independently known mean, not a fitted sample mean.
        let observedCentralSecondMoment: [[Double]]
        let meanAbsoluteLimits: [Double]
        let secondMomentAbsoluteLimits: [[Double]]
    }
    private struct KineticReport: Encodable {
        let degreesOfFreedom: Int
        let sampleCount: Int
        let samplesTwoKOverKBT: [Double]
        let observedMean: Double
        let observedVariance: Double
        let meanLowerDeviationLimit: Double
        let meanUpperDeviationLimit: Double
        let cdfLocations: [Double]
        let expectedCDF: [Double]
        let observedCDF: [Double]
        let cdfAbsoluteLimit: Double
    }
    private struct Molecule: Encodable {
        let name: String
        let massesDa: [Double]
        let offsetsNM: [VivoVector3D]
        let constraintPairs: [[Int]]
        var degreesOfFreedom: Int { 3 * massesDa.count - constraintPairs.count }
    }
    private struct ProjectionOracle: Encodable {
        let molecule: Molecule
        let temperatureK: Double
        let constraintJacobian: [[Double]]
        let inverseConstraintMassMetric: [[Double]]
        let velocityCovarianceNM2PerPS2: [[Double]]
    }
    private func components(_ v: VivoVector3D) -> [Double] { [v.x, v.y, v.z] }
    private func seeds(_ count: Int, namespace: UInt64) -> [UInt64] {
        (0..<count).map { namespace ^ (UInt64($0 + 1) &* 0x9e3779b97f4a7c15) }
    }
    private func model(name: String, masses: [Double], constraints: [VivoDistanceConstraint] = []) throws -> VivoClassicalSystem {
        .init(identifier: name, structureFingerprint: try VivoCanonicalJSON.fingerprint(Data(name.utf8)),
            particles: masses.enumerated().map { index, mass in
                .init(index: UInt32(index), atomIndex: UInt32(index), typeIdentifier: "noninteracting",
                      massDa: mass, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
            }, constraints: constraints)
    }
    private func momentReport(samples: [[Double]], mean: [Double], covariance: [[Double]]) throws -> MomentReport {
        let n = samples.count, dimension = mean.count
        try #require(n > 1 && dimension > 0 && covariance.count == dimension)
        try #require(samples.allSatisfy { $0.count == dimension && $0.allSatisfy(\.isFinite) })
        var observedMean = [Double](repeating: 0, count: dimension)
        var second = [[Double]](repeating: [Double](repeating: 0, count: dimension), count: dimension)
        for sample in samples {
            for i in 0..<dimension {
                observedMean[i] += sample[i] / Double(n)
                for j in 0..<dimension { second[i][j] += (sample[i] - mean[i]) * (sample[j] - mean[j]) / Double(n) }
            }
        }
        let t = policy.logTail, count = Double(n)
        var meanLimits = [Double](repeating: 0, count: dimension), secondLimits = second
        for i in 0..<dimension {
            try #require(covariance[i].count == dimension && covariance[i][i] > 0)
            meanLimits[i] = sqrt(2 * covariance[i][i] * t / count)
                + policy.relativeFP32MomentAllowance * sqrt(covariance[i][i])
            for j in 0..<dimension {
                let product = covariance[i][i] * covariance[j][j], cross = covariance[i][j]
                // Gaussian quadratic-form concentration, including off-diagonal
                // products. Each two-sided failure probability is <=2*exp(-t).
                secondLimits[i][j] = sqrt(2 * (product + cross * cross) * t / count)
                    + (sqrt(product) + abs(cross)) * t / count
                    + policy.relativeFP32MomentAllowance * sqrt(product)
            }
        }
        return .init(sampleCount: n, expectedMean: mean, expectedCovariance: covariance,
            observedMean: observedMean, observedCentralSecondMoment: second,
            meanAbsoluteLimits: meanLimits, secondMomentAbsoluteLimits: secondLimits)
    }
    private func momentComparisonCount(_ dimension: Int) -> Int { dimension + dimension * (dimension + 1) / 2 }
    private func validateFamilyBudget() throws {
        let planned = 3 * momentComparisonCount(3) + [6, 9].reduce(0) { $0 + momentComparisonCount($1) + 4 * 4 }
        try #require(planned == policy.declaredStatisticalComparisons && planned <= policy.comparisonBudget)
    }
    private func check(_ report: MomentReport, statisticalChecks: inout Int) {
        for i in report.expectedMean.indices {
            statisticalChecks += 1
            #expect(abs(report.observedMean[i] - report.expectedMean[i]) <= report.meanAbsoluteLimits[i])
            for j in i..<report.expectedMean.count {
                statisticalChecks += 1
                #expect(abs(report.observedCentralSecondMoment[i][j] - report.expectedCovariance[i][j]) <= report.secondMomentAbsoluteLimits[i][j])
            }
        }
    }
    private func chiSquareCDF(_ value: Double, degrees: Int) -> Double {
        let x = value / 2, target = Double(degrees) / 2
        var shape = degrees % 2 == 0 ? 1.0 : 0.5
        var gamma = degrees % 2 == 0 ? 1.0 : sqrt(Double.pi)
        var result = degrees % 2 == 0 ? 1 - exp(-x) : erf(sqrt(x))
        while shape < target {
            result -= pow(x, shape) * exp(-x) / (shape * gamma)
            gamma *= shape; shape += 1
        }
        return min(1, max(0, result))
    }
    private func kineticReport(_ values: [Double], degrees: Int) throws -> KineticReport {
        try #require(!values.isEmpty && values.allSatisfy { $0.isFinite && $0 >= -1e-12 })
        let count = Double(values.count), mean = values.reduce(0, +) / count
        let locations = [0.5, 1, 2].map { $0 * Double(degrees) }
        return .init(degreesOfFreedom: degrees, sampleCount: values.count, samplesTwoKOverKBT: values,
            observedMean: mean, observedVariance: values.reduce(0) { $0 + ($1 - mean) * ($1 - mean) } / (count - 1),
            meanLowerDeviationLimit: 2 * sqrt(Double(degrees) * policy.logTail / count)
                + policy.relativeFP32MomentAllowance * Double(degrees),
            meanUpperDeviationLimit: 2 * sqrt(Double(degrees) * policy.logTail / count) + 2 * policy.logTail / count
                + policy.relativeFP32MomentAllowance * Double(degrees),
            cdfLocations: locations, expectedCDF: locations.map { chiSquareCDF($0, degrees: degrees) },
            observedCDF: locations.map { point in Double(values.filter { $0 <= point }.count) / count },
            cdfAbsoluteLimit: sqrt(policy.logTail / (2 * count)) + policy.relativeFP32MomentAllowance)
    }
    private func check(_ report: KineticReport, statisticalChecks: inout Int) {
        statisticalChecks += 1 // A two-sided chi-square mean event.
        #expect(report.observedMean >= Double(report.degreesOfFreedom) - report.meanLowerDeviationLimit)
        #expect(report.observedMean <= Double(report.degreesOfFreedom) + report.meanUpperDeviationLimit)
        for i in report.cdfLocations.indices {
            statisticalChecks += 1
            #expect(abs(report.observedCDF[i] - report.expectedCDF[i]) <= report.cdfAbsoluteLimit)
        }
    }

    @Test func freeNVTMatchesFiniteTimeOrnsteinUhlenbeckMoments() async throws {
        let output = try ThermalEvidenceSink(label: "free-ou")
        try validateFamilyBudget()
        let count = 1_024, mass = 12.0, dt = 1.0 / 256, friction = 32.0
        let horizons = [1, 8, 32], seedSet = seeds(8, namespace: 0x4f55_434f_4e46_3031)
        let initialVelocity = VivoVector3D(0.25, -0.125, 0.0625)
        let positions = (0..<count).map { VivoVector3D(Double($0 % 16) * 2, Double(($0 / 16) % 8) * 2, Double($0 / 128) * 2) }
        let system = try model(name: "free-ou-1024", masses: [Double](repeating: mass, count: count))
        struct Design: Encodable {
            let policy: Policy; let particleCount: Int; let massDa: Double; let timeStepPS: Double
            let frictionPerPS: Double; let temperatureK: Double; let initialVelocityNMPerPS: VivoVector3D
            let horizons: [Int]; let seeds: [UInt64]; let independentSamplesPerHorizon: Int
            let physicalBoltzmannKJPerMolK: Double; let nativeFP32A: Double
        }
        try output.record(Design(policy: policy, particleCount: count, massDa: mass, timeStepPS: dt,
            frictionPerPS: friction, temperatureK: temperature, initialVelocityNMPerPS: initialVelocity,
            horizons: horizons, seeds: seedSet, independentSamplesPerHorizon: count * seedSet.count,
            physicalBoltzmannKJPerMolK: kB, nativeFP32A: Double(Float(exp(-friction * dt)))), name: "design")
        let device = try VivoMetalDeviceSelector.productionDevice()
        try #require(device.hasUnifiedMemory)
        var samples: [Int: [[Double]]] = [:]
        var firstSeedFinal: [VivoVector3D]?
        for seed in seedSet {
            let cfg = VivoMDConfiguration(timeStepPS: dt, cutoffNM: 0.125, neighborSkinNM: 0,
                electrostatics: .cutoff, pmeTolerance: nil, pmeGridSpacingNM: nil, ensemble: .nvt,
                thermostat: .langevinMiddle, targetTemperatureK: temperature, frictionPerPS: friction,
                neighborListEnabled: false, randomSeed: seed)
            let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(), positionsNM: positions)
            let runtime = try await VivoMDMetalRuntime.make(system: system, initialState: initial, configuration: cfg,
                initialVelocitiesNMPerPS: [VivoVector3D](repeating: initialVelocity, count: count), device: device)
            var certificates: [VivoMDStepCertificate] = []
            for step in 1...32 {
                let certificate = try await runtime.step()
                certificates.append(certificate)
                if !certificate.committed { try output.record(certificates, name: "rejected-steps-\(seed)") }
                try #require(certificate.committed && certificate.statusFlags == 0)
                if horizons.contains(step) {
                    let snapshot = try await runtime.snapshot()
                    struct Raw: Encodable {
                        let deviceName: String; let deviceRegistryID: UInt64; let seed: UInt64
                        let horizon: Int; let configuration: VivoMDConfiguration; let state: VivoMDStateSnapshot
                    }
                    try output.record(Raw(deviceName: device.name, deviceRegistryID: device.registryID,
                        seed: seed, horizon: step, configuration: cfg, state: snapshot), name: "sample-\(seed)-\(step)")
                    try #require(snapshot.stepIndex == UInt64(step) && snapshot.timePS == Double(step) * dt)
                    for velocity in snapshot.velocitiesNMPerPS { samples[step, default: []].append(components(velocity)) }
                    if step == 32 {
                        if let firstSeedFinal { try #require(snapshot.velocitiesNMPerPS != firstSeedFinal) }
                        else { firstSeedFinal = snapshot.velocitiesNMPerPS }
                    }
                }
            }
            try output.record(certificates, name: "accepted-steps-\(seed)")
        }
        // Reused time horizons are correlated. They are separate tests with
        // N=8*1024 each, never pooled into a fictitious 3N independent sample.
        var statisticalChecks = 0
        for horizon in horizons {
            let observed = try #require(samples[horizon])
            try #require(observed.count == count * seedSet.count)
            let decay = exp(-friction * dt * Double(horizon))
            let variance = kB * temperature / mass * (1 - decay * decay)
            let report = try momentReport(samples: observed, mean: components(initialVelocity * decay),
                covariance: (0..<3).map { i in (0..<3).map { j in i == j ? variance : 0 } })
            try output.record(report, name: "moments-horizon-\(horizon)")
            check(report, statisticalChecks: &statisticalChecks)
        }
        try output.record(["executedStatisticalComparisons": statisticalChecks, "familyPlanned": policy.declaredStatisticalComparisons,
                           "familyBudget": policy.comparisonBudget], name: "comparison-counts")
        try #require(statisticalChecks == 3 * momentComparisonCount(3) && statisticalChecks <= policy.comparisonBudget)
    }

    private func molecule(_ kind: Int) -> Molecule {
        if kind == 0 {
            // The 4:1 ratio gives power against an incorrect Euclidean projector
            // while both inverse masses remain exactly representable in FP32.
            return .init(name: "unequal-mass-dimers", massesDa: [4, 16],
                offsetsNM: [.zero, .init(0.125, 0.09375, 0)], constraintPairs: [[0, 1]])
        }
        return .init(name: "connected-nonlinear-triatomics", massesDa: [12, 16, 20],
            offsetsNM: [.zero, .init(0.125, 0, 0), .init(0.03125, 0.125, 0.0625)],
            constraintPairs: [[0, 1], [0, 2], [1, 2]])
    }
    private func projectionOracle(_ molecule: Molecule) throws -> ProjectionOracle {
        let dimension = molecule.massesDa.count * 3, rank = molecule.constraintPairs.count
        let inverseMass = (0..<dimension).map { 1 / molecule.massesDa[$0 / 3] }
        var g = [[Double]](repeating: [Double](repeating: 0, count: dimension), count: rank)
        for (row, pair) in molecule.constraintPairs.enumerated() {
            let delta = molecule.offsetsNM[pair[0]] - molecule.offsetsNM[pair[1]]
            let direction = components(delta / delta.norm)
            for axis in 0..<3 { g[row][3 * pair[0] + axis] = direction[axis]; g[row][3 * pair[1] + axis] = -direction[axis] }
        }
        var metric = (0..<rank).map { a in (0..<rank).map { b in
            (0..<dimension).reduce(0) { $0 + g[a][$1] * inverseMass[$1] * g[b][$1] }
        } }
        var inverse = (0..<rank).map { a in (0..<rank).map { b in a == b ? 1.0 : 0.0 } }
        // Tiny independent pivoted elimination; no production constraint solver,
        // thermal initializer, or production matrix inverse supplies this oracle.
        for column in 0..<rank {
            let pivot = try #require((column..<rank).max { abs(metric[$0][column]) < abs(metric[$1][column]) })
            try #require(abs(metric[pivot][column]) > 1e-12)
            metric.swapAt(column, pivot); inverse.swapAt(column, pivot)
            let divisor = metric[column][column]
            for j in 0..<rank { metric[column][j] /= divisor; inverse[column][j] /= divisor }
            let pivotMetric = metric[column], pivotInverse = inverse[column]
            for row in 0..<rank where row != column {
                let factor = metric[row][column]
                for j in 0..<rank { metric[row][j] -= factor * pivotMetric[j]; inverse[row][j] -= factor * pivotInverse[j] }
            }
        }
        var covariance = [[Double]](repeating: [Double](repeating: 0, count: dimension), count: dimension)
        for i in 0..<dimension { for j in 0..<dimension {
            var projection = 0.0
            for a in 0..<rank { for b in 0..<rank { projection += g[a][i] * inverse[a][b] * g[b][j] } }
            covariance[i][j] = kB * temperature * ((i == j ? inverseMass[i] : 0) - inverseMass[i] * projection * inverseMass[j])
        } }
        let thermalRank = (0..<dimension).reduce(0) { $0 + covariance[$1][$1] / inverseMass[$1] / (kB * temperature) }
        try #require(abs(thermalRank - Double(molecule.degreesOfFreedom)) <= 1e-10)
        let projector = (0..<dimension).map { i in (0..<dimension).map { j in
            covariance[i][j] / (kB * temperature * sqrt(inverseMass[i] * inverseMass[j]))
        } }
        for i in 0..<dimension { for j in 0..<dimension {
            try #require(abs(projector[i][j] - projector[j][i]) <= 1e-10)
            let squared = (0..<dimension).reduce(0) { $0 + projector[i][$1] * projector[$1][j] }
            try #require(abs(squared - projector[i][j]) <= 1e-10)
        } }
        for a in 0..<rank { for j in 0..<dimension {
            let tangent = (0..<dimension).reduce(0) { $0 + g[a][$1] * covariance[$1][j] }
            try #require(abs(tangent) <= 1e-11)
        } }
        return .init(molecule: molecule, temperatureK: temperature, constraintJacobian: g,
                     inverseConstraintMassMetric: inverse, velocityCovarianceNM2PerPS2: covariance)
    }

    @Test(arguments: [0, 1])
    func fixedGeometryThermalizationMatchesProjectedGaussian(kind: Int) async throws {
        try validateFamilyBudget()
        let molecule = molecule(kind), copies = 128, particlesPerCopy = molecule.massesDa.count
        let output = try ThermalEvidenceSink(label: molecule.name)
        let oracle = try projectionOracle(molecule), seedSet = seeds(32, namespace: 0x5448_4552_4d41_3031)
        let cfg = VivoMDConfiguration(timeStepPS: 0.002, cutoffNM: 0.0625, neighborSkinNM: 0,
            electrostatics: .cutoff, pmeTolerance: nil, pmeGridSpacingNM: nil, ensemble: .nve,
            thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil,
            constraintTolerance: policy.relativeDistanceLimit, maximumConstraintIterations: policy.maximumConstraintIterations,
            neighborListEnabled: false)
        struct Design: Encodable {
            let policy: Policy; let oracle: ProjectionOracle; let copies: Int; let seeds: [UInt64]
            let configuration: VivoMDConfiguration; let independentMolecularSamples: Int; let boltzmannKJPerMolK: Double
        }
        try output.record(Design(policy: policy, oracle: oracle, copies: copies, seeds: seedSet, configuration: cfg,
            independentMolecularSamples: copies * seedSet.count, boltzmannKJPerMolK: kB), name: "design")
        var positions: [VivoVector3D] = [], masses: [Double] = [], constraints: [VivoDistanceConstraint] = []
        for copy in 0..<copies {
            let origin = VivoVector3D(Double(copy % 8) * 2, Double((copy / 8) % 4) * 2, Double(copy / 32) * 2)
            positions += molecule.offsetsNM.map { $0 + origin }; masses += molecule.massesDa
            for pair in molecule.constraintPairs {
                constraints.append(.init(a: UInt32(copy * particlesPerCopy + pair[0]), b: UInt32(copy * particlesPerCopy + pair[1]),
                    distanceNM: (molecule.offsetsNM[pair[0]] - molecule.offsetsNM[pair[1]]).norm))
            }
        }
        let system = try model(name: molecule.name, masses: masses, constraints: constraints)
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(), positionsNM: positions, sourceTimePS: 0.25)
        let device = try VivoMetalDeviceSelector.productionDevice()
        try #require(device.hasUnifiedMemory)
        let runtime = try await VivoMDMetalRuntime.make(system: system, initialState: initial, configuration: cfg, device: device)
        let initialCheckpoint = try await runtime.checkpoint()
        // No virtual sites or position projection occurs during initialization
        // for this profile. Bind the independent Jacobian to actual readback.
        try #require(initialCheckpoint.positionsNM == positions)
        var samples: [[Double]] = [], totalX: [Double] = [], centerX: [Double] = [], rotationalX: [Double] = []
        var globalCenterX: [Double] = []
        for seed in seedSet {
            let checkpoint: VivoMDCheckpoint
            do { checkpoint = try await runtime.thermalize(temperatureK: temperature, seed: seed) }
            catch {
                struct Failure: Encodable { let seed: UInt64; let error: String; let retained: VivoMDCheckpoint? }
                let retained = try? await runtime.checkpoint()
                try output.record(Failure(seed: seed, error: String(describing: error), retained: retained), name: "thermalize-failed-\(seed)")
                throw error
            }
            let committed = try await runtime.checkpoint()
            var maximumDistanceError = 0.0, maximumTangentSpeed = 0.0
            for constraint in constraints {
                let a = Int(constraint.a), b = Int(constraint.b)
                let delta = checkpoint.positionsNM[a] - checkpoint.positionsNM[b]
                maximumDistanceError = max(maximumDistanceError, abs(delta.norm - constraint.distanceNM) / constraint.distanceNM)
                maximumTangentSpeed = max(maximumTangentSpeed,
                    abs(delta.dot(checkpoint.velocitiesNMPerPS[a] - checkpoint.velocitiesNMPerPS[b])) / delta.norm)
            }
            struct Raw: Encodable {
                let deviceName: String; let deviceRegistryID: UInt64; let seed: UInt64
                let returnedCheckpoint: VivoMDCheckpoint; let committedCheckpoint: VivoMDCheckpoint
                let maximumRelativeDistanceError: Double; let maximumTangentSpeedNMPerPS: Double
            }
            try output.record(Raw(deviceName: device.name, deviceRegistryID: device.registryID, seed: seed,
                returnedCheckpoint: checkpoint, committedCheckpoint: committed,
                maximumRelativeDistanceError: maximumDistanceError, maximumTangentSpeedNMPerPS: maximumTangentSpeed), name: "sample-\(seed)")
            try #require(checkpoint == committed)
            try #require(checkpoint.positionsNM == initialCheckpoint.positionsNM && checkpoint.periodicCell == initialCheckpoint.periodicCell)
            try #require(checkpoint.acceptedStep == initialCheckpoint.acceptedStep && checkpoint.timePS == initialCheckpoint.timePS)
            #expect(maximumDistanceError <= policy.relativeDistanceLimit)
            #expect(maximumTangentSpeed <= policy.tangentSpeedLimitNMPerPS)
            // Local marginals have little power against subtracting the global
            // COM. One independent whole-system momentum sample per seed tests
            // the retained global translational modes without claiming N=4096.
            let globalMomentum = zip(checkpoint.velocitiesNMPerPS, masses).reduce(VivoVector3D.zero) { $0 + $1.0 * $1.1 }
            globalCenterX.append(globalMomentum.squaredNorm / (masses.reduce(0, +) * kB * temperature))
            for copy in 0..<copies {
                let velocities = Array(checkpoint.velocitiesNMPerPS[(copy * particlesPerCopy)..<((copy + 1) * particlesPerCopy)])
                samples.append(velocities.flatMap(components))
                let totalMass = molecule.massesDa.reduce(0, +)
                let momentum = zip(velocities, molecule.massesDa).reduce(VivoVector3D.zero) { $0 + $1.0 * $1.1 }
                let centerVelocity = momentum / totalMass
                let kinetic = zip(velocities, molecule.massesDa).reduce(0.0) { $0 + $1.0.squaredNorm * $1.1 } / (kB * temperature)
                let center = totalMass * centerVelocity.squaredNorm / (kB * temperature)
                let rotation = zip(velocities, molecule.massesDa).reduce(0.0) { $0 + ($1.0 - centerVelocity).squaredNorm * $1.1 } / (kB * temperature)
                totalX.append(kinetic); centerX.append(center); rotationalX.append(rotation)
            }
        }
        try #require(samples.count == copies * seedSet.count)
        let moments = try momentReport(samples: samples, mean: [Double](repeating: 0, count: particlesPerCopy * 3),
                                       covariance: oracle.velocityCovarianceNM2PerPS2)
        let kinetic = try kineticReport(totalX, degrees: molecule.degreesOfFreedom)
        let center = try kineticReport(centerX, degrees: 3)
        let rotation = try kineticReport(rotationalX, degrees: molecule.degreesOfFreedom - 3)
        try #require(globalCenterX.count == seedSet.count)
        let globalCenter = try kineticReport(globalCenterX, degrees: 3)
        try output.record(moments, name: "projected-velocity-moments")
        try output.record(kinetic, name: "total-kinetic-chi-square")
        try output.record(center, name: "center-of-mass-kinetic-chi-square")
        try output.record(rotation, name: "rotational-kinetic-chi-square")
        try output.record(globalCenter, name: "global-center-of-mass-kinetic-chi-square")
        var statisticalChecks = 0
        check(moments, statisticalChecks: &statisticalChecks)
        check(kinetic, statisticalChecks: &statisticalChecks)
        check(center, statisticalChecks: &statisticalChecks)
        check(rotation, statisticalChecks: &statisticalChecks)
        check(globalCenter, statisticalChecks: &statisticalChecks)
        try output.record(["executedStatisticalComparisons": statisticalChecks, "familyPlanned": policy.declaredStatisticalComparisons,
                           "familyBudget": policy.comparisonBudget], name: "comparison-counts")
        try #require(statisticalChecks == momentComparisonCount(particlesPerCopy * 3) + 4 * 4 && statisticalChecks <= policy.comparisonBudget)
    }
}

/// Unique retained raw observations, never a success receipt. With no artifact
/// directory, the complete JSON goes to the test log instead.
private struct ThermalEvidenceSink {
    let root: URL?
    init(label: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { root = nil; return }
        guard !path.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoArtifactValidationError.invalid("NUMIVIVO_TEST_ARTIFACTS must name a nonempty directory")
        }
        let parent = URL(fileURLWithPath: path, isDirectory: true).standardizedFileURL
        try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: true)
        let directory = parent.appendingPathComponent("classical-thermal-\(label)-\(UUID().uuidString)")
        guard Darwin.mkdir(directory.path, 0o700) == 0 else {
            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: [NSFilePathErrorKey: directory.path])
        }
        root = directory
        print("NUMIVIVO_THERMAL_EVIDENCE_ROOT=\(directory.path)")
    }
    func record<T: Encodable>(_ value: T, name: String) throws {
        let data = try VivoCanonicalJSON.encode(value)
        if let root {
            let destination = root.appendingPathComponent(name + ".json")
            try data.write(to: destination, options: .withoutOverwriting)
            print("NUMIVIVO_THERMAL_OBSERVATION=\(destination.path)")
        } else { print("NUMIVIVO_THERMAL_OBSERVATION[\(name)]=\(String(decoding: data, as: UTF8.self))") }
    }
}
