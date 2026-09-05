import Foundation

/// A published cohort summary, not independent patient observations or an assay
/// confidence interval. The range is descriptive. No implicit time conversion
/// from a clinical visit label to simulation seconds is made.
public struct VivoAggregateOccupancyObservation: Codable, Equatable, Sendable {
    public let identifier: String
    public let compound: String
    public let target: String
    public let tissue: String
    public let visit: String
    public let regimenLabel: String
    public let patientCount: Int
    public let medianFraction: Double
    public let minimumFraction: Double
    public let maximumFraction: Double
    public let strictThreshold: Double
    public let reportedFractionAboveThreshold: Double
    public let sourceLocator: String
    public init(identifier: String, compound: String, target: String, tissue: String, visit: String,
                regimenLabel: String, patientCount: Int, medianFraction: Double, minimumFraction: Double,
                maximumFraction: Double, strictThreshold: Double, reportedFractionAboveThreshold: Double,
                sourceLocator: String) {
        self.identifier = identifier; self.compound = compound; self.target = target; self.tissue = tissue
        self.visit = visit; self.regimenLabel = regimenLabel; self.patientCount = patientCount
        self.medianFraction = medianFraction; self.minimumFraction = minimumFraction; self.maximumFraction = maximumFraction
        self.strictThreshold = strictThreshold; self.reportedFractionAboveThreshold = reportedFractionAboveThreshold
        self.sourceLocator = sourceLocator
    }
    public func validate() throws {
        guard [identifier, compound, target, tissue, visit, regimenLabel, sourceLocator].allSatisfy({ !$0.isEmpty && $0.utf8.count <= 4096 }),
              (1...10000).contains(patientCount), [medianFraction, minimumFraction, maximumFraction, strictThreshold, reportedFractionAboveThreshold]
                .allSatisfy({ $0.isFinite && (0...1).contains($0) }), minimumFraction <= medianFraction, medianFraction <= maximumFraction else {
            throw VivoPosteriorError.invalid("aggregate occupancy observation")
        }
    }
}

public struct VivoAggregateOccupancyBenchmark: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let identifier: String
    public let sourceDOI: String
    public let sourceURI: String
    public let transcriptionMethod: String
    public let sourceAvailability: String
    public let observations: [VivoAggregateOccupancyObservation]
    public let missingInputs: [String]
    public func validate() throws {
        guard schemaVersion == 1, !identifier.isEmpty, !sourceDOI.isEmpty,
              sourceURI.hasPrefix("https://"), !transcriptionMethod.isEmpty, !sourceAvailability.isEmpty,
              (1...256).contains(observations.count), Set(observations.map(\.identifier)).count == observations.count else {
            throw VivoPosteriorError.invalid("aggregate dataset provenance or shape")
        }
        for item in observations { try item.validate() }
    }
}

public struct VivoPredictedOccupancyCohort: Codable, Equatable, Sendable {
    public let observationIdentifier: String
    public let compound: String
    public let target: String
    public let tissue: String
    public let visit: String
    public let regimenLabel: String
    /// Each row is a complete virtual cohort of the observed sample size. Each
    /// column within that row is an individual. Posterior particles for ONE
    /// individual must never be relabelled as a population of individuals.
    public let cohortDraws: [[Double]]
    public let populationModelEvidence: VivoKineticEvidence
    public let timingMappingEvidence: VivoKineticEvidence
    public let predictionEvidence: VivoKineticEvidence
}

public struct VivoAggregateBenchmarkRequest: Codable, Equatable, Sendable {
    public let benchmark: VivoAggregateOccupancyBenchmark
    public let predictions: [VivoPredictedOccupancyCohort]
}

public struct VivoAggregateBenchmarkComparison: Codable, Equatable, Sendable {
    public let observationIdentifier: String
    public let observedPatientCount: Int
    public let virtualCohortDrawCount: Int
    public let observedMedian: Double
    public let expectedCohortMedian: Double
    public let cohortMedian90PercentInterval: [Double]
    public let medianDifference: Double
    public let observedFractionAboveThreshold: Double
    public let expectedCohortFractionAboveThreshold: Double
    public let cohortFraction90PercentInterval: [Double]
    public let thresholdFractionDifference: Double
    public let publishedRange: [Double]
}

public struct VivoAggregateBenchmarkReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let requestFingerprint: VivoFingerprint
    public let comparisons: [VivoAggregateBenchmarkComparison]
    public let limitations: [String]
}

public enum VivoAggregateOccupancyEvaluator {
    public static func compare(_ request: VivoAggregateBenchmarkRequest) throws -> VivoAggregateBenchmarkReport {
        try request.benchmark.validate()
        let expectedIDs = Set(request.benchmark.observations.map(\.identifier))
        guard request.predictions.count == expectedIDs.count,
              Set(request.predictions.map(\.observationIdentifier)) == expectedIDs else {
            throw VivoPosteriorError.invalid("all benchmark cohorts must be supplied exactly once")
        }
        let lookup = Dictionary(uniqueKeysWithValues: request.predictions.map { ($0.observationIdentifier, $0) })
        var comparisons: [VivoAggregateBenchmarkComparison] = [], totalValues = 0
        for observation in request.benchmark.observations {
            try Task.checkCancellation()
            let prediction = lookup[observation.identifier]!
            for evidence in [prediction.populationModelEvidence, prediction.timingMappingEvidence, prediction.predictionEvidence] {
                try evidence.validate(origin: .calculated)
            }
            guard prediction.compound == observation.compound, prediction.target == observation.target,
                  prediction.tissue == observation.tissue, prediction.visit == observation.visit,
                  prediction.regimenLabel == observation.regimenLabel, (1...4096).contains(prediction.cohortDraws.count) else {
                throw VivoPosteriorError.invalid("prediction context or virtual cohort shape")
            }
            var medians: [Double] = [], fractions: [Double] = []
            for row in prediction.cohortDraws {
                totalValues += row.count
                guard totalValues <= 2_097_152, row.count == observation.patientCount,
                      row.allSatisfy({ $0.isFinite && (0...1).contains($0) }) else {
                    throw VivoPosteriorError.invalid("virtual individuals, sample size or cohort work capacity")
                }
                medians.append(VivoPosteriorNumerics.quantile(sorted: row.sorted(), probability: 0.5))
                fractions.append(Double(row.filter { $0 > observation.strictThreshold }.count) / Double(row.count))
            }
            let expectedMedian = medians.reduce(0) { $0 + $1 / Double(medians.count) }
            let expectedFraction = fractions.reduce(0) { $0 + $1 / Double(fractions.count) }
            medians.sort(); fractions.sort()
            comparisons.append(.init(observationIdentifier: observation.identifier, observedPatientCount: observation.patientCount,
                virtualCohortDrawCount: medians.count, observedMedian: observation.medianFraction,
                expectedCohortMedian: expectedMedian,
                cohortMedian90PercentInterval: [VivoPosteriorNumerics.quantile(sorted: medians, probability: 0.05), VivoPosteriorNumerics.quantile(sorted: medians, probability: 0.95)],
                medianDifference: expectedMedian - observation.medianFraction,
                observedFractionAboveThreshold: observation.reportedFractionAboveThreshold,
                expectedCohortFractionAboveThreshold: expectedFraction,
                cohortFraction90PercentInterval: [VivoPosteriorNumerics.quantile(sorted: fractions, probability: 0.05), VivoPosteriorNumerics.quantile(sorted: fractions, probability: 0.95)],
                thresholdFractionDifference: expectedFraction - observation.reportedFractionAboveThreshold,
                publishedRange: [observation.minimumFraction, observation.maximumFraction]))
        }
        return try .init(schemaVersion: 1, requestFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),
            comparisons: comparisons, limitations: [
                "Descriptive comparison with published aggregates, not an individual-observation likelihood or a validation certificate.",
                "Published minima/maxima are ranges, not confidence limits or standard deviations. Rounded percentages remain rounded reports.",
                "Virtual cohort draws must reflect a declared population model and visit/exposure mapping. One person's posterior uncertainty is not population variability.",
                "No raw patient values, missing exposure curves, biochemical rates or measurement SD are reconstructed from these summaries.",
                "Model intervals describe its declared virtual cohort distribution; clinical efficacy, dosing choice and treatment safety are not assessed."
            ])
    }
}
