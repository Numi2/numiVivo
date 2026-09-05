import Foundation

public struct VivoExperimentRunOutput: Sendable {
    public let result: VivoResultPack
    public let checkpoints: [VivoExperimentCheckpoint]

    public init(result: VivoResultPack, checkpoints: [VivoExperimentCheckpoint]) {
        self.result = result
        self.checkpoints = checkpoints
    }
}

public actor VivoExperimentRunner {
    public struct Options: Sendable, Codable, Equatable {
        public var runtimeVersion: String
        public var maximumStepsPerReplicate: UInt64
        public var maximumStoredSamples: Int
        public var maximumStoredEvents: Int
        public var maximumStoredCheckpoints: Int
        public var failOnDroppedEvents: Bool
        public var boundaryToleranceSeconds: Double

        public init(
            runtimeVersion: String = "NumiVivo/0.1.0",
            maximumStepsPerReplicate: UInt64 = 10_000_000,
            maximumStoredSamples: Int = 10_000_000,
            maximumStoredEvents: Int = 1_000_000,
            maximumStoredCheckpoints: Int = 10_000,
            failOnDroppedEvents: Bool = true,
            boundaryToleranceSeconds: Double = 1e-7
        ) {
            self.runtimeVersion = runtimeVersion
            self.maximumStepsPerReplicate = maximumStepsPerReplicate
            self.maximumStoredSamples = maximumStoredSamples
            self.maximumStoredEvents = maximumStoredEvents
            self.maximumStoredCheckpoints = maximumStoredCheckpoints
            self.failOnDroppedEvents = failOnDroppedEvents
            self.boundaryToleranceSeconds = boundaryToleranceSeconds
        }

        func validate() throws {
            guard !runtimeVersion.isEmpty,
                  maximumStepsPerReplicate > 0,
                  maximumStoredSamples > 0,
                  maximumStoredEvents > 0,
                  maximumStoredCheckpoints > 0,
                  boundaryToleranceSeconds.isFinite,
                  boundaryToleranceSeconds >= 0 else {
                throw VivoArtifactValidationError.invalid("experiment runner options are invalid")
            }
        }
    }

    private let programPack: VivoProgramPack
    private let hostContext: PreparedVivoHostContext
    private let experiment: PreparedVivoExperiment
    private let coupling: PreparedVivoCouplingPlan?
    private let baseConfiguration: VivoRuntimeConfiguration
    private let options: Options
    private let stepPlanner: VivoAdaptiveStepPlanner

    public init(
        programPack: VivoProgramPack,
        hostContext: PreparedVivoHostContext,
        experiment: PreparedVivoExperiment,
        coupling: PreparedVivoCouplingPlan? = nil,
        baseConfiguration: VivoRuntimeConfiguration,
        options: Options = .init()
    ) throws {
        try options.validate()
        guard programPack.header.contentFingerprint == experiment.programFingerprint,
              hostContext.programFingerprint == experiment.programFingerprint,
              hostContext.contextFingerprint == experiment.hostContextFingerprint,
              baseConfiguration.environmentCount == experiment.environmentCount,
              baseConfiguration.fidelity == experiment.fidelity else {
            throw VivoArtifactValidationError.incompatible("runner inputs disagree on program, context, environment count, or fidelity")
        }
        if let expected = experiment.couplingFingerprint {
            guard coupling?.fingerprint == expected else {
                throw VivoArtifactValidationError.incompatible("runner coupling plan does not match the experiment")
            }
        }
        self.programPack = programPack
        self.hostContext = hostContext
        self.experiment = experiment
        self.coupling = coupling
        self.baseConfiguration = baseConfiguration
        self.options = options
        self.stepPlanner = try VivoAdaptiveStepPlanner()
    }

    public func run() async throws -> VivoExperimentRunOutput {
        let startedAt = Date()
        let monitorMetadata = try programPack.monitorMetadata()
        let stopConditions = Dictionary(
            uniqueKeysWithValues: experiment.stopConditions.map { ($0.monitorIdentifier, $0) }
        )
        var allSamples: [VivoMeasurementSample] = []
        var allSummaries: [VivoMeasurementSummary] = []
        var allEvents: [VivoRecordedEvent] = []
        var allCheckpoints: [VivoExperimentCheckpoint] = []
        var allLedgers: [VivoRunLedger] = []
        var replicateResults: [VivoReplicateResult] = []
        var limitations = Set<String>()
        var deviceName = "uninitialized"
        var deviceRegistryID: UInt64 = 0

        if experiment.fidelity == .tissue, coupling == nil {
            limitations.insert("F4 was executed without an external coupling plan; only the local reaction–transport state advanced.")
        }

        for replicateIndex in 0..<experiment.replicateCount {
            var configuration = baseConfiguration
            configuration.seed = experiment.replicateSeeds[Int(replicateIndex)]
            configuration.timeStep = try checkedFloat(experiment.preferredTimeStepSeconds, label: "preferred time step")
            configuration.minimumTimeStep = try checkedFloat(experiment.minimumTimeStepSeconds, label: "minimum time step")
            configuration.maximumTimeStep = try checkedFloat(experiment.maximumTimeStepSeconds, label: "maximum time step")
            let runtime = try await VivoRuntime.make(pack: programPack, configuration: configuration)
            if replicateIndex == 0 {
                deviceName = runtime.capabilities.deviceName
                deviceRegistryID = runtime.capabilities.registryID
            }
            _ = try await runtime.apply(hostContext)

            let ledgerGenesis = try makeLedgerGenesis(replicateIndex: replicateIndex, seed: configuration.seed)
            var ledger = VivoRunLedger(genesis: ledgerGenesis)
            var schedules = experiment.measurements.map(MeasurementSchedule.init)
            var summaries = Dictionary(
                uniqueKeysWithValues: experiment.measurements.map {
                    ($0.identifier, MeasurementAccumulator(measurement: $0, replicateIndex: replicateIndex))
                }
            )
            var interventionCursor = 0
            var currentTransport = hostContext.transport
            var nextCheckpointTime = experiment.checkpointCadenceSeconds
            var checkpointFingerprints: [VivoFingerprint] = []
            var committedSteps: UInt64 = 0
            var rejectedSteps: UInt64 = 0
            var stepAttempts: UInt64 = 0
            var replicateSampleCount: UInt64 = 0
            var termination: VivoReplicateResult.Termination = .completedDuration
            var terminationReason = "configured duration reached"

            try await applyDueInterventions(
                runtime: runtime,
                currentTime: 0,
                cursor: &interventionCursor,
                transport: &currentTransport
            )
            _ = try await runtime.primePendingInputs()
            let initialSnapshot = try await runtime.snapshot()
            try collectSnapshotMeasurements(
                snapshot: initialSnapshot,
                replicateIndex: replicateIndex,
                schedules: &schedules,
                summaries: &summaries,
                samples: &allSamples,
                replicateSampleCount: &replicateSampleCount,
                limitations: &limitations
            )

            while Double(await runtime.time()) + options.boundaryToleranceSeconds < experiment.durationSeconds {
                guard stepAttempts < options.maximumStepsPerReplicate else {
                    termination = .runtimeFailure
                    terminationReason = "maximumStepsPerReplicate was reached"
                    break
                }
                stepAttempts &+= 1
                let currentTime = Double(await runtime.time())

                try await applyDueInterventions(
                    runtime: runtime,
                    currentTime: currentTime,
                    cursor: &interventionCursor,
                    transport: &currentTransport
                )
                let lifecycle = await runtime.lifecycle()
                switch lifecycle {
                case .ready:
                    break
                case .reversiblyStopped(let reason):
                    termination = .reversibleShutdown
                    terminationReason = reason
                    break
                case .permanentlyStopped(let reason):
                    termination = .permanentShutdown
                    terminationReason = reason
                    break
                case .failed(let reason):
                    termination = .runtimeFailure
                    terminationReason = reason
                    break
                }
                if termination != .completedDuration { break }

                let boundaries = exactBoundaries(
                    currentTime: currentTime,
                    interventionCursor: interventionCursor,
                    schedules: schedules,
                    nextCheckpointTime: nextCheckpointTime
                )
                let plan = try stepPlanner.plan(
                    requestedSeconds: experiment.preferredTimeStepSeconds,
                    pack: programPack,
                    configuration: configuration,
                    transport: currentTransport,
                    exactBoundaries: boundaries
                )
                let plannedAfter = min(currentTime + plan.selectedSeconds, experiment.durationSeconds)
                let publicationPlan = makePublicationPlan(
                    schedules: schedules,
                    plannedAfter: plannedAfter
                )
                let request = VivoStepRequest(
                    timeStep: try checkedFloat(plannedAfter - currentTime, label: "planned step"),
                    coupling: [],
                    publications: publicationPlan.requests,
                    permitAdaptiveReduction: true
                )
                let stepResult = try await runtime.step(request)
                try ledger.append(.init(replicateIndex: replicateIndex, certificate: stepResult.certificate))
                if ledger.entries.count > Int(options.maximumStepsPerReplicate) {
                    throw VivoArtifactValidationError.invalid("step ledger exceeded configured limit")
                }

                if options.failOnDroppedEvents, stepResult.certificate.status.eventDropped > 0 {
                    termination = .runtimeFailure
                    terminationReason = "runtime event buffer overflowed"
                }
                guard allEvents.count <= options.maximumStoredEvents - stepResult.events.count else {
                    throw VivoArtifactValidationError.invalid("stored event limit exceeded")
                }
                allEvents.append(contentsOf: stepResult.events.map {
                    VivoRecordedEvent(replicateIndex: replicateIndex, event: $0)
                })

                switch stepResult.certificate.disposition {
                case .committed, .committedWithReducedStep:
                    committedSteps &+= 1
                    let actualAfter = Double(stepResult.certificate.timeAfter)
                    try collectPublishedMeasurements(
                        publicationPlan: publicationPlan,
                        result: stepResult,
                        actualTime: actualAfter,
                        replicateIndex: replicateIndex,
                        schedules: &schedules,
                        summaries: &summaries,
                        samples: &allSamples,
                        replicateSampleCount: &replicateSampleCount,
                        limitations: &limitations
                    )
                    if let monitorIndex = stepResult.certificate.status.firstMonitor,
                       Int(monitorIndex) < monitorMetadata.count,
                       let stop = stopConditions[monitorMetadata[Int(monitorIndex)].identifier] {
                        switch stop.response {
                        case .finishReplicate:
                            termination = .stopCondition
                            terminationReason = "stop condition \(stop.identifier) fired"
                        case .rejectReplicate:
                            termination = .rejectedStep
                            terminationReason = "stop condition \(stop.identifier) rejected the replicate"
                        case .reversibleShutdown:
                            termination = .reversibleShutdown
                            terminationReason = "stop condition \(stop.identifier) requested reversible shutdown"
                        case .permanentShutdown:
                            termination = .permanentShutdown
                            terminationReason = "stop condition \(stop.identifier) requested permanent shutdown"
                        }
                    }
                    while let checkpointTime = nextCheckpointTime,
                          actualAfter + options.boundaryToleranceSeconds >= checkpointTime {
                        guard allCheckpoints.count < options.maximumStoredCheckpoints else {
                            throw VivoArtifactValidationError.invalid("stored checkpoint limit exceeded")
                        }
                        let checkpoint = try await VivoExperimentCheckpoint(
                            state: runtime.resumeCheckpoint(), replicateIndex: replicateIndex,
                            experimentFingerprint: experiment.fingerprint,
                            hostContextFingerprint: hostContext.contextFingerprint,
                            couplingFingerprint: coupling?.fingerprint,
                            ledgerHead: ledger.head
                        )
                        let fingerprint = try checkpoint.fingerprint()
                        checkpointFingerprints.append(fingerprint)
                        allCheckpoints.append(checkpoint)
                        nextCheckpointTime = checkpointTime + (experiment.checkpointCadenceSeconds ?? .infinity)
                    }
                case .rejected:
                    rejectedSteps &+= 1
                    termination = .rejectedStep
                    terminationReason = "runtime rejected step \(stepResult.certificate.stepIndex)"
                case .reversibleShutdown:
                    termination = .reversibleShutdown
                    terminationReason = "runtime safety system requested reversible shutdown"
                case .permanentShutdown:
                    termination = .permanentShutdown
                    terminationReason = "runtime safety system requested permanent shutdown"
                }
                if termination != .completedDuration { break }
            }

            for measurement in experiment.measurements {
                if let accumulator = summaries[measurement.identifier], let summary = accumulator.finalize() {
                    allSummaries.append(summary)
                }
            }
            allLedgers.append(ledger)
            let finalTime = Double(await runtime.time())
            replicateResults.append(VivoReplicateResult(
                replicateIndex: replicateIndex,
                seed: configuration.seed,
                termination: termination,
                terminationReason: terminationReason,
                finalStepIndex: await runtime.stepIndex(),
                finalTimeSeconds: finalTime,
                ledgerHead: ledger.head,
                committedSteps: committedSteps,
                rejectedSteps: rejectedSteps,
                measurementCount: replicateSampleCount,
                checkpointFingerprints: checkpointFingerprints
            ))
        }

        let result = VivoResultPack(
            experimentFingerprint: experiment.fingerprint,
            programFingerprint: programPack.header.contentFingerprint,
            hostContextFingerprint: hostContext.contextFingerprint,
            couplingFingerprint: coupling?.fingerprint,
            sourceProgramFingerprint: programPack.header.sourceFingerprint,
            fidelity: experiment.fidelity,
            runtimeVersion: options.runtimeVersion,
            deviceName: deviceName,
            deviceRegistryID: deviceRegistryID,
            startedAt: startedAt,
            finishedAt: Date(),
            replicates: replicateResults,
            measurements: allSamples,
            measurementSummaries: allSummaries,
            events: allEvents.map(\.event),
            recordedEvents: allEvents,
            ledgers: allLedgers,
            limitations: limitations.sorted(),
            annotations: [:]
        )
        try result.validate()
        return VivoExperimentRunOutput(result: result, checkpoints: allCheckpoints)
    }

    private func applyDueInterventions(
        runtime: VivoRuntime,
        currentTime: Double,
        cursor: inout Int,
        transport: inout [VivoSpeciesTransportABI]
    ) async throws {
        while cursor < experiment.interventions.count,
              experiment.interventions[cursor].timeSeconds <= currentTime + options.boundaryToleranceSeconds {
            let intervention = experiment.interventions[cursor]
            if case .transport(let index, let value) = intervention.operation {
                guard Int(index) < transport.count else {
                    throw VivoArtifactValidationError.invalid("prepared transport intervention index is out of bounds")
                }
                transport[Int(index)] = value
            }
            try await runtime.apply(intervention: intervention.operation)
            cursor += 1
        }
    }

    private func exactBoundaries(
        currentTime: Double,
        interventionCursor: Int,
        schedules: [MeasurementSchedule],
        nextCheckpointTime: Double?
    ) -> [(source: VivoStepLimit.Source, timeRemaining: Double, subject: String)] {
        var result: [(VivoStepLimit.Source, Double, String)] = [
            (.experimentBoundary, experiment.durationSeconds - currentTime, "experiment.end")
        ]
        if interventionCursor < experiment.interventions.count {
            let intervention = experiment.interventions[interventionCursor]
            result.append((.interventionBoundary, intervention.timeSeconds - currentTime, intervention.identifier))
        }
        if let nextMeasurement = schedules.compactMap(\.nextTime).filter({ $0 > currentTime + options.boundaryToleranceSeconds }).min() {
            result.append((.measurementBoundary, nextMeasurement - currentTime, "measurement"))
        }
        if let nextCheckpointTime,
           nextCheckpointTime > currentTime + options.boundaryToleranceSeconds {
            result.append((.checkpointBoundary, nextCheckpointTime - currentTime, "checkpoint"))
        }
        if let coupling {
            for channel in coupling.inbound + coupling.outbound where channel.cadenceSeconds > 0 {
                let periods = floor(currentTime / channel.cadenceSeconds) + 1
                let next = periods * channel.cadenceSeconds
                result.append((.couplingCadence, next - currentTime, channel.identifier))
            }
        }
        return result
    }

    private func makePublicationPlan(
        schedules: [MeasurementSchedule],
        plannedAfter: Double
    ) -> PublicationPlan {
        var requests: [VivoPublicationRequest] = []
        var slices: [String: Range<Int>] = [:]
        for schedule in schedules {
            let measurement = schedule.measurement
            let scheduled = schedule.nextTime.map {
                $0 <= plannedAfter + options.boundaryToleranceSeconds
            } ?? false
            guard scheduled || measurement.storage == .eventOnly else { continue }
            let start = requests.count
            requests.append(contentsOf: measurement.lanes.map {
                VivoPublicationRequest(speciesIndex: measurement.speciesIndex, laneIndex: $0)
            })
            slices[measurement.identifier] = start..<requests.count
        }
        return PublicationPlan(requests: requests, slices: slices)
    }

    private func collectSnapshotMeasurements(
        snapshot: VivoStateSnapshot,
        replicateIndex: UInt32,
        schedules: inout [MeasurementSchedule],
        summaries: inout [String: MeasurementAccumulator],
        samples: inout [VivoMeasurementSample],
        replicateSampleCount: inout UInt64,
        limitations: inout Set<String>
    ) throws {
        for index in schedules.indices {
            guard let nextTime = schedules[index].nextTime,
                  nextTime <= options.boundaryToleranceSeconds,
                  schedules[index].measurement.storage != .eventOnly else { continue }
            let measurement = schedules[index].measurement
            let raw = measurement.lanes.compactMap {
                snapshot.value(species: measurement.speciesIndex, lane: $0)
            }.map(Double.init)
            guard raw.count == measurement.lanes.count else {
                throw VivoArtifactValidationError.invalid("snapshot measurement lane resolution failed")
            }
            let values = try transformAndAggregate(raw, measurement: measurement, limitations: &limitations)
            try storeMeasurement(
                measurement: measurement,
                values: values,
                time: 0,
                stepIndex: snapshot.stepIndex,
                replicateIndex: replicateIndex,
                summaries: &summaries,
                samples: &samples,
                replicateSampleCount: &replicateSampleCount
            )
            schedules[index].advance(past: 0, tolerance: options.boundaryToleranceSeconds)
        }
    }

    private func collectPublishedMeasurements(
        publicationPlan: PublicationPlan,
        result: VivoStepResult,
        actualTime: Double,
        replicateIndex: UInt32,
        schedules: inout [MeasurementSchedule],
        summaries: inout [String: MeasurementAccumulator],
        samples: inout [VivoMeasurementSample],
        replicateSampleCount: inout UInt64,
        limitations: inout Set<String>
    ) throws {
        for index in schedules.indices {
            let measurement = schedules[index].measurement
            guard let slice = publicationPlan.slices[measurement.identifier],
                  slice.upperBound <= result.publications.count else { continue }
            let due = schedules[index].nextTime.map {
                $0 <= actualTime + options.boundaryToleranceSeconds
            } ?? false
            let eventTriggered = measurement.storage == .eventOnly && !result.events.isEmpty && due
            guard due && (measurement.storage != .eventOnly || eventTriggered) else { continue }
            let raw = result.publications[slice].map(Double.init)
            let values = try transformAndAggregate(raw, measurement: measurement, limitations: &limitations)
            try storeMeasurement(
                measurement: measurement,
                values: values,
                time: actualTime,
                stepIndex: result.certificate.stepIndex,
                replicateIndex: replicateIndex,
                summaries: &summaries,
                samples: &samples,
                replicateSampleCount: &replicateSampleCount
            )
            schedules[index].advance(past: actualTime, tolerance: options.boundaryToleranceSeconds)
        }
    }

    private func transformAndAggregate(
        _ raw: [Double],
        measurement: PreparedVivoMeasurement,
        limitations: inout Set<String>
    ) throws -> [Double] {
        let transformed = try raw.map(measurement.transform.apply)
        guard !transformed.isEmpty else {
            throw VivoArtifactValidationError.invalid("measurement \(measurement.identifier) produced no values")
        }
        switch measurement.aggregation {
        case .direct:
            return transformed
        case .mean:
            return [transformed.reduce(0, +) / Double(transformed.count)]
        case .sum:
            return [transformed.reduce(0, +)]
        case .minimum:
            return [transformed.min()!]
        case .maximum:
            return [transformed.max()!]
        case .volumeWeightedMean:
            limitations.insert("volumeWeightedMean currently uses equal voxel volume; dynamic cut-cell volume weights must be supplied by the coupled spatial participant.")
            return [transformed.reduce(0, +) / Double(transformed.count)]
        }
    }

    private func storeMeasurement(
        measurement: PreparedVivoMeasurement,
        values: [Double],
        time: Double,
        stepIndex: UInt32,
        replicateIndex: UInt32,
        summaries: inout [String: MeasurementAccumulator],
        samples: inout [VivoMeasurementSample],
        replicateSampleCount: inout UInt64
    ) throws {
        guard var accumulator = summaries[measurement.identifier] else {
            throw VivoArtifactValidationError.unresolved("measurement accumulator is missing")
        }
        try accumulator.append(values: values, time: time)
        summaries[measurement.identifier] = accumulator
        replicateSampleCount &+= 1

        guard measurement.storage != .summaryOnly else { return }
        guard samples.count < options.maximumStoredSamples else {
            throw VivoArtifactValidationError.invalid("stored measurement sample limit exceeded")
        }
        samples.append(VivoMeasurementSample(
            measurementIdentifier: measurement.identifier,
            replicateIndex: replicateIndex,
            stepIndex: stepIndex,
            timeSeconds: time,
            values: values,
            unit: measurement.outputUnit,
            laneIndices: measurement.aggregation == .direct ? measurement.lanes : nil
        ))
    }

    private func makeLedgerGenesis(replicateIndex: UInt32, seed: VivoRuntimeSeed) throws -> VivoFingerprint {
        struct Genesis: Codable {
            let experiment: VivoFingerprint
            let program: VivoFingerprint
            let host: VivoFingerprint
            let replicate: UInt32
            let seed: VivoRuntimeSeed
        }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Genesis(
            experiment: experiment.fingerprint,
            program: programPack.header.contentFingerprint,
            host: hostContext.contextFingerprint,
            replicate: replicateIndex,
            seed: seed
        )))
    }

    private func checkedFloat(_ value: Double, label: String) throws -> Float {
        guard value.isFinite,
              value >= 0,
              value <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid("\(label) cannot be represented as non-negative finite FP32")
        }
        return Float(value)
    }
}

private struct PublicationPlan {
    let requests: [VivoPublicationRequest]
    let slices: [String: Range<Int>]
}

private struct MeasurementSchedule {
    let measurement: PreparedVivoMeasurement
    private(set) var nextTime: Double?

    init(measurement: PreparedVivoMeasurement) {
        self.measurement = measurement
        self.nextTime = measurement.startSeconds <= measurement.endSeconds
            ? measurement.startSeconds
            : nil
    }

    mutating func advance(past time: Double, tolerance: Double) {
        guard var candidate = nextTime else { return }
        repeat {
            candidate += measurement.cadenceSeconds
        } while candidate <= time + tolerance
        nextTime = candidate <= measurement.endSeconds + tolerance ? candidate : nil
    }
}

private struct MeasurementAccumulator {
    let measurement: PreparedVivoMeasurement
    let replicateIndex: UInt32
    private var count: UInt64 = 0
    private var firstTime: Double = 0
    private var lastTime: Double = 0
    private var sum: [Double] = []
    private var minimum: [Double] = []
    private var maximum: [Double] = []
    private var last: [Double] = []

    init(measurement: PreparedVivoMeasurement, replicateIndex: UInt32) {
        self.measurement = measurement
        self.replicateIndex = replicateIndex
    }

    mutating func append(values: [Double], time: Double) throws {
        guard time.isFinite, values.allSatisfy(\.isFinite), !values.isEmpty else {
            throw VivoArtifactValidationError.invalid("measurement accumulator received invalid values")
        }
        if count == 0 {
            firstTime = time
            sum = values
            minimum = values
            maximum = values
        } else {
            guard values.count == sum.count else {
                throw VivoArtifactValidationError.invalid("measurement vector width changed during a replicate")
            }
            for index in values.indices {
                sum[index] += values[index]
                minimum[index] = Swift.min(minimum[index], values[index])
                maximum[index] = Swift.max(maximum[index], values[index])
            }
        }
        last = values
        lastTime = time
        count &+= 1
    }

    func finalize() -> VivoMeasurementSummary? {
        guard count > 0 else { return nil }
        return VivoMeasurementSummary(
            measurementIdentifier: measurement.identifier,
            replicateIndex: replicateIndex,
            sampleCount: count,
            firstTimeSeconds: firstTime,
            lastTimeSeconds: lastTime,
            mean: sum.map { $0 / Double(count) },
            minimum: minimum,
            maximum: maximum,
            last: last,
            unit: measurement.outputUnit
        )
    }
}
