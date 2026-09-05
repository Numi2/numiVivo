import Foundation

public enum VivoOccupancyObservable: String, Codable, Sendable {
    case drugOccupancy, covalentOccupancy, freeTargetFraction, competitorOccupancy
    public func value(_ sample: VivoTargetEngagementSample) -> Double {
        switch self {
        case .drugOccupancy: return sample.drugOccupancy
        case .covalentOccupancy: return sample.covalentOccupancy
        case .freeTargetFraction: return sample.fractions.free / sample.fractions.total
        case .competitorOccupancy: return sample.fractions.competitor / sample.fractions.total
        }
    }
}

public enum VivoStudyPartition: String, Codable, Sendable, CaseIterable { case calibration, validation, test }

/// Observations are fractions, not percentages. Values may fall outside [0,1]
/// through measurement noise; neither observations nor predictions are clipped.
public struct VivoOccupancyObservation: Codable, Equatable, Sendable {
    public let identifier: String
    public let timeSeconds: Double
    public let observable: VivoOccupancyObservable
    public let value: Double
    public let standardDeviation: Double?
    public let origin: VivoKineticOrigin
    public let evidence: VivoKineticEvidence
    public init(identifier: String, timeSeconds: Double, observable: VivoOccupancyObservable,
                value: Double, standardDeviation: Double? = nil, origin: VivoKineticOrigin,
                evidence: VivoKineticEvidence) {
        self.identifier = identifier; self.timeSeconds = timeSeconds; self.observable = observable
        self.value = value; self.standardDeviation = standardDeviation; self.origin = origin; self.evidence = evidence
    }
}

public struct VivoTargetEngagementStudyCase: Codable, Equatable, Sendable {
    public let identifier: String
    /// All measurements from one donor/condition/assay grouping must remain in
    /// one partition. Defining this grouping is an explicit dataset responsibility.
    public let leakageGroup: String
    public let partition: VivoStudyPartition
    public let experiment: VivoTargetEngagementExperiment
    public let observations: [VivoOccupancyObservation]
    public init(identifier: String, leakageGroup: String, partition: VivoStudyPartition,
                experiment: VivoTargetEngagementExperiment, observations: [VivoOccupancyObservation]) {
        self.identifier = identifier; self.leakageGroup = leakageGroup; self.partition = partition
        self.experiment = experiment; self.observations = observations
    }
}

public struct VivoTargetEngagementStudy: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let identifier: String
    public let cases: [VivoTargetEngagementStudyCase]
    public init(identifier: String, cases: [VivoTargetEngagementStudyCase]) {
        schemaVersion = 1; self.identifier = identifier; self.cases = cases
    }
    public func validate() throws {
        guard schemaVersion == 1, !identifier.isEmpty, identifier.utf8.count <= 1024,
              (1...1024).contains(cases.count) else { throw VivoKineticsError.invalid("study identity or capacity") }
        var ids = Set<String>(), groups: [String: VivoStudyPartition] = [:]
        var totalSamples = 0, totalObservations = 0
        for item in cases {
            try item.experiment.validate()
            guard !item.identifier.isEmpty, item.identifier.utf8.count <= 1024,
                  !item.leakageGroup.isEmpty, item.leakageGroup.utf8.count <= 1024,
                  ids.insert(item.identifier).inserted, !item.observations.isEmpty else {
                throw VivoKineticsError.invalid("study case identity, grouping or observations")
            }
            if let previous = groups[item.leakageGroup], previous != item.partition {
                throw VivoKineticsError.invalid("one leakage group cannot occur in multiple partitions")
            }
            groups[item.leakageGroup] = item.partition
            totalSamples += item.experiment.sampleTimesSeconds.count
            totalObservations += item.observations.count
            guard totalSamples <= 262_144, totalObservations <= 262_144 else { throw VivoKineticsError.capacity("study samples") }
            let times = Set(item.experiment.sampleTimesSeconds)
            var observationIDs = Set<String>()
            for observation in item.observations {
                try observation.evidence.validate(origin: observation.origin)
                guard !observation.identifier.isEmpty, observation.identifier.utf8.count <= 1024,
                      observationIDs.insert(observation.identifier).inserted,
                      observation.value.isFinite, times.contains(observation.timeSeconds) else {
                    throw VivoKineticsError.invalid("observation identity or exact sampling time")
                }
                if let sd = observation.standardDeviation {
                    guard sd.isFinite, sd > 0 else { throw VivoKineticsError.invalid("observation SD must be positive or explicitly unknown") }
                }
            }
        }
    }
}

public struct VivoOccupancyResidual: Codable, Equatable, Sendable {
    public let observationIdentifier: String
    public let observable: VivoOccupancyObservable
    public let timeSeconds: Double
    public let predicted: Double
    public let observed: Double
    public let error: Double
    public let standardizedError: Double?
}

public struct VivoTargetStudyCaseResult: Codable, Equatable, Sendable {
    public let identifier: String
    public let partition: VivoStudyPartition
    public let prediction: VivoTargetEngagementResult?
    public let residuals: [VivoOccupancyResidual]
    public let failure: String?
}

public struct VivoTargetStudyPartitionResult: Codable, Equatable, Sendable {
    public let partition: VivoStudyPartition
    public let requestedCases: Int
    public let failedCases: Int
    public let observationCount: Int
    public let rootMeanSquaredError: Double?
    public let meanAbsoluteError: Double?
    public let measuredObservationCount: Int
}

public struct VivoTargetEngagementStudyResult: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let study: VivoTargetEngagementStudy
    public let numerics: VivoTargetEngagementNumerics
    public let cases: [VivoTargetStudyCaseResult]
    public let partitions: [VivoTargetStudyPartitionResult]
    public let limitations: [String]
}

/// Evaluation only: no fitting, no posterior claims, no held-out observations
/// used to select parameters. All requested cases remain in the output. A failed
/// case makes that partition's aggregate metrics unavailable, not more favorable.
public enum VivoTargetEngagementStudyEvaluator {
    public static func evaluate(_ study: VivoTargetEngagementStudy,
                                numerics: VivoTargetEngagementNumerics = .init()) throws -> VivoTargetEngagementStudyResult {
        try study.validate(); try numerics.validate()
        var results: [VivoTargetStudyCaseResult] = []
        for item in study.cases {
            try Task.checkCancellation()
            do {
                let prediction = try VivoTargetEngagementReference.run(item.experiment, numerics: numerics)
                let byTime = Dictionary(uniqueKeysWithValues: prediction.samples.map { ($0.timeSeconds, $0) })
                var residuals: [VivoOccupancyResidual] = []
                for observation in item.observations {
                    guard let sample = byTime[observation.timeSeconds] else { throw VivoKineticsError.numerical("missing requested observation") }
                    let value = observation.observable.value(sample)
                    let error = value - observation.value
                    let standardized = observation.standardDeviation.map { error / $0 }
                    guard value.isFinite, error.isFinite, (error * error).isFinite,
                          standardized?.isFinite != false else { throw VivoKineticsError.numerical("residual overflow") }
                    residuals.append(.init(observationIdentifier: observation.identifier,
                        observable: observation.observable, timeSeconds: observation.timeSeconds,
                        predicted: value, observed: observation.value, error: error, standardizedError: standardized))
                }
                results.append(.init(identifier: item.identifier, partition: item.partition,
                                     prediction: prediction, residuals: residuals, failure: nil))
            } catch is CancellationError { throw CancellationError() }
            catch {
                results.append(.init(identifier: item.identifier, partition: item.partition,
                                     prediction: nil, residuals: [], failure: error.localizedDescription))
            }
        }
        let partitions = VivoStudyPartition.allCases.map { role -> VivoTargetStudyPartitionResult in
            let input = study.cases.filter { $0.partition == role }
            let output = results.filter { $0.partition == role }
            let failed = output.filter { $0.failure != nil }.count
            let residuals = output.flatMap(\.residuals)
            let n = Double(residuals.count)
            // Scale each contribution before adding to avoid overflow in a large
            // sum of individually representable squared errors.
            let mse = n > 0 ? residuals.reduce(0) { $0 + ($1.error * $1.error) / n } : Double.nan
            let mae = n > 0 ? residuals.reduce(0) { $0 + abs($1.error) / n } : Double.nan
            let eligible = failed == 0 && !input.isEmpty && mse.isFinite && mae.isFinite
            return .init(partition: role, requestedCases: input.count, failedCases: failed,
                observationCount: input.reduce(0) { $0 + $1.observations.count },
                rootMeanSquaredError: eligible ? sqrt(mse) : nil, meanAbsoluteError: eligible ? mae : nil,
                measuredObservationCount: input.flatMap(\.observations).filter { $0.origin == .measured }.count)
        }
        return .init(schemaVersion: 1, study: study, numerics: numerics, cases: results, partitions: partitions,
            limitations: ["Nominal occupancy validation only; no cellular efficacy or whole-body safety endpoint is inferred.",
                "Partitions and leakage groups are declared by the dataset author; grouping cannot discover undisclosed shared provenance.",
                "Standardized residuals condition on supplied measurement SD; they are not predictive coverage estimates.",
                "Synthetic/calculated observations are not experimental biological validation.",
                "Numerical failures remain visible and suppress aggregate metrics for the affected partition."])
    }
}
