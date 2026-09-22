import Foundation
import MLX
import MLXNN
import MLXOptimizers
import NumiVivoKit

/// Errors raised by the MLX response learner before an artifact is published.
public enum VivoCellResponseLearningError: Error, LocalizedError, Sendable {
    case invalid(String)
    case incompatible(String)
    case limit(String)

    public var errorDescription: String? {
        switch self {
        case .invalid(let value): return "cell-response learning: \(value)"
        case .incompatible(let value): return "cell-response learning incompatible: \(value)"
        case .limit(let value): return "cell-response learning resource limit: \(value)"
        }
    }
}

/// The bounded neural architecture used for one source-bound response model.
/// Context cells, perturbation identity, and supplied target descriptors are
/// fused before separate mean and variance heads produce log-CPM outputs.
public struct VivoCellResponseArchitecture: Codable, Sendable, Equatable {
    public let hiddenWidth: Int
    public let contextCells: Int
    public let maximumParameters: Int

    public init(hiddenWidth: Int = 256, contextCells: Int = 32, maximumParameters: Int = 30_000_000) {
        self.hiddenWidth = hiddenWidth
        self.contextCells = contextCells
        self.maximumParameters = maximumParameters
    }

    private enum CodingKeys: String, CodingKey { case hiddenWidth, contextCells, maximumParameters }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["hiddenWidth", "contextCells", "maximumParameters"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        hiddenWidth = try values.decode(Int.self, forKey: .hiddenWidth)
        contextCells = try values.decode(Int.self, forKey: .contextCells)
        maximumParameters = try values.decode(Int.self, forKey: .maximumParameters)
    }

    private static func product(_ left: Int, _ right: Int) throws -> Int {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow, result.partialValue >= 0 else {
            throw VivoCellResponseLearningError.limit("parameter arithmetic")
        }
        return result.partialValue
    }

    private static func sum(_ values: [Int]) throws -> Int {
        var total = 0
        for value in values {
            let result = total.addingReportingOverflow(value)
            guard !result.overflow else { throw VivoCellResponseLearningError.limit("parameter arithmetic") }
            total = result.partialValue
        }
        return total
    }

    func parameterCount(featureCount: Int, targetCount: Int, descriptorCount: Int) throws -> Int {
        let projection = try Self.product(featureCount, hiddenWidth)
        let embedding = try Self.product(targetCount, hiddenWidth)
        // A target-specific response prior is fitted exclusively from the
        // training partition. It gives known perturbations a measured starting
        // response instead of asking a randomly initialized full-gene decoder
        // to discover a 33k-dimensional effect from a few updates.
        let responsePrior = try Self.product(targetCount, featureCount)
        let descriptor = try Self.product(descriptorCount, hiddenWidth)
        let fusion = try Self.product(try Self.product(hiddenWidth, 3), hiddenWidth)
        let output = try Self.product(try Self.product(featureCount, hiddenWidth), 2)
        return try Self.sum([
            projection, hiddenWidth,
            embedding, responsePrior,
            descriptor, hiddenWidth,
            fusion, hiddenWidth,
            output, try Self.product(featureCount, 2)
        ])
    }

    func validate(featureCount: Int, targetCount: Int, descriptorCount: Int) throws {
        guard (16...1_024).contains(hiddenWidth), (1...512).contains(contextCells),
              (100_000...30_000_000).contains(maximumParameters),
              featureCount > 0, targetCount > 0, descriptorCount > 0 else {
            throw VivoCellResponseLearningError.invalid("architecture dimensions or parameter budget")
        }
        let count = try parameterCount(featureCount: featureCount, targetCount: targetCount, descriptorCount: descriptorCount)
        guard count <= maximumParameters else {
            throw VivoCellResponseLearningError.limit("architecture has \(count) parameters, above \(maximumParameters)")
        }
    }
}

/// A deterministic training run. The learner uses momentum-free SGD so a
/// checkpoint can resume exactly from the stored model, absolute step, and
/// stateless sampler seed.
public struct VivoCellResponseTrainingPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let architecture: VivoCellResponseArchitecture
    public let steps: Int
    public let batchSize: Int
    public let learningRate: Double
    public let weightDecay: Double
    public let seed: UInt64
    public let validationEvery: Int
    public let validationExamples: Int

    public init(id: String, architecture: VivoCellResponseArchitecture = .init(), steps: Int = 1_000,
                batchSize: Int = 4, learningRate: Double = 1e-3, weightDecay: Double = 1e-5,
                seed: UInt64 = 0, validationEvery: Int = 100, validationExamples: Int = 32) {
        schemaVersion = 1
        self.id = id
        self.architecture = architecture
        self.steps = steps
        self.batchSize = batchSize
        self.learningRate = learningRate
        self.weightDecay = weightDecay
        self.seed = seed
        self.validationEvery = validationEvery
        self.validationExamples = validationExamples
    }

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, id, architecture, steps, batchSize, learningRate, weightDecay, seed, validationEvery, validationExamples
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "id", "architecture", "steps", "batchSize", "learningRate", "weightDecay", "seed", "validationEvery", "validationExamples"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        architecture = try values.decode(VivoCellResponseArchitecture.self, forKey: .architecture)
        steps = try values.decode(Int.self, forKey: .steps)
        batchSize = try values.decode(Int.self, forKey: .batchSize)
        learningRate = try values.decode(Double.self, forKey: .learningRate)
        weightDecay = try values.decode(Double.self, forKey: .weightDecay)
        seed = try values.decode(UInt64.self, forKey: .seed)
        validationEvery = try values.decode(Int.self, forKey: .validationEvery)
        validationExamples = try values.decode(Int.self, forKey: .validationExamples)
    }

    func validateStatic() throws {
        guard schemaVersion == 1, !id.isEmpty, id.utf8.count <= 1_024,
              (1...1_000_000).contains(steps), (1...128).contains(batchSize),
              learningRate.isFinite, (1e-8...1).contains(learningRate),
              weightDecay.isFinite, (0...1).contains(weightDecay),
              (1...steps).contains(validationEvery), (1...4_096).contains(validationExamples) else {
            throw VivoCellResponseLearningError.invalid("training plan fields")
        }
        try architecture.validate(featureCount: 1, targetCount: 1, descriptorCount: 1)
    }

    func validate(for corpus: VivoCellResponseCompositeCorpusReader) throws {
        try validateStatic()
        try architecture.validate(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                  descriptorCount: corpus.descriptorCount)
        let denseValues = Double(batchSize) * Double(architecture.contextCells) * Double(corpus.featureCount)
        guard denseValues.isFinite, denseValues <= 134_217_728 else {
            throw VivoCellResponseLearningError.limit("dense context batch exceeds 512 MiB")
        }
        let trainingExamples = try corpus.examples(in: .training)
        let grouped = Dictionary(grouping: trainingExamples, by: \.stratumIndex)
        let strata = try grouped.keys.sorted().map { index -> [VivoCellResponseTrainingExample] in
            guard let rows = grouped[index], !rows.isEmpty else {
                throw VivoCellResponseLearningError.invalid("empty response stratum")
            }
            return rows.sorted { $0.targetRow < $1.targetRow }
        }
        let sourceStrata = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: strata)
        guard let largestSourceStrata = sourceStrata.map(\.count).max() else {
            throw VivoCellResponseLearningError.invalid("training response strata")
        }
        let coverageSlots = largestSourceStrata.multipliedReportingOverflow(by: sourceStrata.count)
        guard !coverageSlots.overflow else {
            throw VivoCellResponseLearningError.limit("training source schedule")
        }
        let scheduled = steps.multipliedReportingOverflow(by: batchSize)
        guard !scheduled.overflow, scheduled.partialValue >= coverageSlots.partialValue else {
            throw VivoCellResponseLearningError.invalid("training plan cannot cover every source response stratum")
        }
        // A training-only development corpus has no validation partition. The
        // run loop already treats that as an intentional absence; only impose
        // per-source validation coverage when a declared treated validation
        // sample exists. Once it does, propagate reader errors rather than
        // mistaking a malformed held-out partition for an intentional absence.
        let hasValidationTreatment = corpus.sourceReaders.contains { reader in
            reader.plan.assignments.contains {
                $0.role == .perturbed && $0.partition == .validation
            }
        }
        if hasValidationTreatment {
            let validationRows = try corpus.examples(in: .validation)
            let validationGroups = Dictionary(grouping: validationRows, by: \.stratumIndex)
            let validationStrata = try validationGroups.keys.sorted().map { index -> [VivoCellResponseTrainingExample] in
                guard let rows = validationGroups[index], !rows.isEmpty else {
                    throw VivoCellResponseLearningError.invalid("empty validation response stratum")
                }
                return rows.sorted { $0.targetRow < $1.targetRow }
            }
            _ = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: validationStrata)
            guard validationExamples >= corpus.sourceCount else {
                throw VivoCellResponseLearningError.invalid(
                    "training validation example budget cannot represent every source")
            }
        }
        _ = try corpus.trainingTargetBindings()
    }

    func validate(for corpus: VivoCellResponseCorpusReader) throws {
        try validate(for: VivoCellResponseCompositeCorpusReader(readers: [corpus]))
    }
}

/// A resume request can only add deterministic steps. Architecture, optimizer
/// settings, and sampler seed remain bound to the parent checkpoint.
public struct VivoCellResponseResumePlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let additionalSteps: Int

    public init(id: String, additionalSteps: Int) {
        schemaVersion = 1
        self.id = id
        self.additionalSteps = additionalSteps
    }

    private enum CodingKeys: String, CodingKey { case schemaVersion, id, additionalSteps }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "additionalSteps"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        additionalSteps = try values.decode(Int.self, forKey: .additionalSteps)
    }

    func validate() throws {
        guard schemaVersion == 1, !id.isEmpty, id.utf8.count <= 1_024,
              (1...1_000_000).contains(additionalSteps) else {
            throw VivoCellResponseLearningError.invalid("resume plan fields")
        }
    }
}

/// A prediction asks for a response distribution from declared control cells.
/// A target can use a trained target embedding only when its exact ID was in
/// the fitted corpus; novel targets use their supplied descriptor vector alone.
public struct VivoCellResponsePredictionPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let target: VivoCellResponseTarget
    public let contextSampleIDs: [String]
    /// Required for a multi-source composite so equal sample labels in two
    /// original studies cannot be silently combined into one query context.
    public let contextCorpus: VivoFingerprint?
    public let useTrainedTargetEmbedding: Bool

    public init(id: String, target: VivoCellResponseTarget, contextSampleIDs: [String],
                contextCorpus: VivoFingerprint? = nil,
                useTrainedTargetEmbedding: Bool = false) {
        schemaVersion = contextCorpus == nil ? 1 : 2
        self.id = id
        self.target = target
        self.contextSampleIDs = contextSampleIDs
        self.contextCorpus = contextCorpus
        self.useTrainedTargetEmbedding = useTrainedTargetEmbedding
    }

    private enum CodingKeys: String, CodingKey { case schemaVersion, id, target, contextSampleIDs, contextCorpus, useTrainedTargetEmbedding }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "target", "contextSampleIDs", "contextCorpus", "useTrainedTargetEmbedding"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        target = try values.decode(VivoCellResponseTarget.self, forKey: .target)
        contextSampleIDs = try values.decode([String].self, forKey: .contextSampleIDs)
        contextCorpus = try values.decodeIfPresent(VivoFingerprint.self, forKey: .contextCorpus)
        useTrainedTargetEmbedding = try values.decode(Bool.self, forKey: .useTrainedTargetEmbedding)
    }

    public func encode(to encoder: Encoder) throws {
        var values = encoder.container(keyedBy: CodingKeys.self)
        try values.encode(schemaVersion, forKey: .schemaVersion)
        try values.encode(id, forKey: .id)
        try values.encode(target, forKey: .target)
        try values.encode(contextSampleIDs, forKey: .contextSampleIDs)
        if let contextCorpus {
            try values.encode(contextCorpus, forKey: .contextCorpus)
        }
        try values.encode(useTrainedTargetEmbedding, forKey: .useTrainedTargetEmbedding)
    }

    public func validate(expectedDescriptorCount: Int? = nil) throws {
        try target.validate(expectedDescriptorCount: expectedDescriptorCount)
        guard (schemaVersion == 1 && contextCorpus == nil) || (schemaVersion == 2 && contextCorpus != nil) else {
            throw VivoCellResponseLearningError.invalid("prediction query")
        }
        guard vivoOmicsID(id), !contextSampleIDs.isEmpty,
              contextSampleIDs.count <= 4_096, Set(contextSampleIDs).count == contextSampleIDs.count,
              contextSampleIDs.allSatisfy(vivoOmicsID), target.descriptors.allSatisfy({ Float($0).isFinite }) else {
            throw VivoCellResponseLearningError.invalid("prediction query")
        }
    }
}

public struct VivoCellResponseTrainingMetrics: Codable, Sendable, Equatable {
    public let step: UInt64
    public let trainNegativeLogLikelihood: Double
    public let validationNegativeLogLikelihood: Double?
    public let validationRMSE: Double?
    /// RMSE of the exact matched-control mean on the same frozen validation
    /// rows as `validationRMSE`. This is a gate, not a fitted metric.
    public let validationMatchedControlRMSE: Double?
}

/// State carried in the native checkpoint. All fields that affect sampling or
/// numerical continuation are explicit rather than inferred from a filename.
public struct VivoCellResponseModelState: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let architecture: VivoCellResponseArchitecture
    public let featureAxis: VivoCellResponseFeatureAxis
    public let targetIDs: [String]
    public let trainedTargetBindings: [VivoCellResponseTrainedTargetBinding]
    public let descriptorCount: Int
    public let corpus: VivoFingerprint
    public let trainingPlan: VivoFingerprint
    /// Canonical source-replayable cohort contract carried by current models.
    /// Legacy v5 artifacts omit it and remain inspectable only.
    public let cohort: VivoFingerprint?
    public let step: UInt64
    public let seed: UInt64
    public let samplerVersion: UInt32
    /// Algorithm used to compute the frozen training-only response prior.
    public let targetResponsePriorVersion: UInt32
    /// Deterministic optimizer routing used for the learned response decoder.
    public let optimizerVersion: UInt32
    public let learningRate: Double
    public let weightDecay: Double
    public let latestMetrics: VivoCellResponseTrainingMetrics

    public init(schemaVersion: Int, architecture: VivoCellResponseArchitecture,
                featureAxis: VivoCellResponseFeatureAxis, targetIDs: [String],
                trainedTargetBindings: [VivoCellResponseTrainedTargetBinding], descriptorCount: Int,
                corpus: VivoFingerprint, trainingPlan: VivoFingerprint, cohort: VivoFingerprint? = nil,
                step: UInt64, seed: UInt64,
                samplerVersion: UInt32, targetResponsePriorVersion: UInt32,
                optimizerVersion: UInt32,
                learningRate: Double, weightDecay: Double,
                latestMetrics: VivoCellResponseTrainingMetrics) {
        self.schemaVersion = schemaVersion
        self.architecture = architecture
        self.featureAxis = featureAxis
        self.targetIDs = targetIDs
        self.trainedTargetBindings = trainedTargetBindings
        self.descriptorCount = descriptorCount
        self.corpus = corpus
        self.trainingPlan = trainingPlan
        self.cohort = cohort
        self.step = step
        self.seed = seed
        self.samplerVersion = samplerVersion
        self.targetResponsePriorVersion = targetResponsePriorVersion
        self.optimizerVersion = optimizerVersion
        self.learningRate = learningRate
        self.weightDecay = weightDecay
        self.latestMetrics = latestMetrics
    }

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, architecture, featureAxis, targetIDs, trainedTargetBindings, descriptorCount,
             corpus, trainingPlan, cohort, step, seed, samplerVersion, targetResponsePriorVersion,
             optimizerVersion,
             learningRate, weightDecay, latestMetrics
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "architecture", "featureAxis", "targetIDs", "trainedTargetBindings", "descriptorCount",
            "corpus", "trainingPlan", "cohort", "step", "seed", "samplerVersion", "targetResponsePriorVersion",
            "optimizerVersion",
            "learningRate", "weightDecay", "latestMetrics"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        architecture = try values.decode(VivoCellResponseArchitecture.self, forKey: .architecture)
        featureAxis = try values.decode(VivoCellResponseFeatureAxis.self, forKey: .featureAxis)
        targetIDs = try values.decode([String].self, forKey: .targetIDs)
        trainedTargetBindings = try values.decode([VivoCellResponseTrainedTargetBinding].self, forKey: .trainedTargetBindings)
        descriptorCount = try values.decode(Int.self, forKey: .descriptorCount)
        corpus = try values.decode(VivoFingerprint.self, forKey: .corpus)
        trainingPlan = try values.decode(VivoFingerprint.self, forKey: .trainingPlan)
        cohort = try values.decodeIfPresent(VivoFingerprint.self, forKey: .cohort)
        step = try values.decode(UInt64.self, forKey: .step)
        seed = try values.decode(UInt64.self, forKey: .seed)
        samplerVersion = try values.decode(UInt32.self, forKey: .samplerVersion)
        targetResponsePriorVersion = try values.decode(UInt32.self, forKey: .targetResponsePriorVersion)
        optimizerVersion = try values.decode(UInt32.self, forKey: .optimizerVersion)
        learningRate = try values.decode(Double.self, forKey: .learningRate)
        weightDecay = try values.decode(Double.self, forKey: .weightDecay)
        latestMetrics = try values.decode(VivoCellResponseTrainingMetrics.self, forKey: .latestMetrics)
    }

    public func encode(to encoder: Encoder) throws {
        var values = encoder.container(keyedBy: CodingKeys.self)
        try values.encode(schemaVersion, forKey: .schemaVersion)
        try values.encode(architecture, forKey: .architecture)
        try values.encode(featureAxis, forKey: .featureAxis)
        try values.encode(targetIDs, forKey: .targetIDs)
        try values.encode(trainedTargetBindings, forKey: .trainedTargetBindings)
        try values.encode(descriptorCount, forKey: .descriptorCount)
        try values.encode(corpus, forKey: .corpus)
        try values.encode(trainingPlan, forKey: .trainingPlan)
        try values.encodeIfPresent(cohort, forKey: .cohort)
        try values.encode(step, forKey: .step)
        try values.encode(seed, forKey: .seed)
        try values.encode(samplerVersion, forKey: .samplerVersion)
        try values.encode(targetResponsePriorVersion, forKey: .targetResponsePriorVersion)
        try values.encode(optimizerVersion, forKey: .optimizerVersion)
        try values.encode(learningRate, forKey: .learningRate)
        try values.encode(weightDecay, forKey: .weightDecay)
        try values.encode(latestMetrics, forKey: .latestMetrics)
    }
}

public struct VivoCellResponseModelReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let checkpoint: VivoFingerprint
    public let trainingPlan: VivoFingerprint
    public let corpus: VivoFingerprint
    /// Present only for v6 cohort-bound models.
    public let cohort: VivoFingerprint?
    public let implementation: VivoFingerprint

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, format, checkpoint, trainingPlan, corpus, cohort, implementation
    }

    public init(schemaVersion: Int, format: String, checkpoint: VivoFingerprint,
                trainingPlan: VivoFingerprint, corpus: VivoFingerprint, cohort: VivoFingerprint? = nil,
                implementation: VivoFingerprint) {
        self.schemaVersion = schemaVersion
        self.format = format
        self.checkpoint = checkpoint
        self.trainingPlan = trainingPlan
        self.corpus = corpus
        self.cohort = cohort
        self.implementation = implementation
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "format", "checkpoint", "trainingPlan", "corpus", "cohort", "implementation"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        format = try values.decode(String.self, forKey: .format)
        checkpoint = try values.decode(VivoFingerprint.self, forKey: .checkpoint)
        trainingPlan = try values.decode(VivoFingerprint.self, forKey: .trainingPlan)
        corpus = try values.decode(VivoFingerprint.self, forKey: .corpus)
        cohort = try values.decodeIfPresent(VivoFingerprint.self, forKey: .cohort)
        implementation = try values.decode(VivoFingerprint.self, forKey: .implementation)
    }

    public func encode(to encoder: Encoder) throws {
        var values = encoder.container(keyedBy: CodingKeys.self)
        try values.encode(schemaVersion, forKey: .schemaVersion)
        try values.encode(format, forKey: .format)
        try values.encode(checkpoint, forKey: .checkpoint)
        try values.encode(trainingPlan, forKey: .trainingPlan)
        try values.encode(corpus, forKey: .corpus)
        try values.encodeIfPresent(cohort, forKey: .cohort)
        try values.encode(implementation, forKey: .implementation)
    }
}

/// Distributional prediction on the corpus feature axis. `featureMask` marks
/// axis features measured in the query source; values where it is zero are
/// structural absence rather than a measured zero-expression prediction.
public struct VivoCellResponsePrediction: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let id: String
    public let featureAxis: VivoCellResponseFeatureAxis
    public let featureMask: [Float]
    /// The deterministic aggregate of the declared matched control cells.
    /// Keeping it separate makes the learned perturbation response auditable.
    public let contextBaselineLogCPM: [Float]
    /// The training-only target response prior plus learned residual over
    /// `contextBaselineLogCPM` (zero for descriptor-only queries).
    public let meanDeltaLogCPM: [Float]
    /// Exact source rows used to form the published control baseline.
    public let contextSourceRows: [Int]
    public let meanLogCPM: [Float]
    public let varianceLogCPM: [Float]
}

public struct VivoCellResponsePredictionReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let model: VivoFingerprint
    public let corpus: VivoFingerprint
    public let query: VivoFingerprint
    public let prediction: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// One immutable held-out target and the exact control rows used to score it.
public struct VivoCellResponseEvaluationExample: Codable, Sendable, Equatable {
    public let targetRow: Int
    public let contextRows: [Int]
}

/// A source-specific held-out score. Current composite evaluations macro-average
/// these values so a high-cell-count or wider-panel source cannot hide a weak
/// response fit in another bound study.
public struct VivoCellResponseSourceEvaluation: Codable, Sendable, Equatable {
    public let corpus: VivoFingerprint
    public let examples: Int
    public let observedFeatures: Int
    public let negativeLogLikelihood: Double
    public let rmse: Double
    public let matchedControlRMSE: Double

    public init(corpus: VivoFingerprint, examples: Int, observedFeatures: Int,
                negativeLogLikelihood: Double, rmse: Double, matchedControlRMSE: Double) {
        self.corpus = corpus
        self.examples = examples
        self.observedFeatures = observedFeatures
        self.negativeLogLikelihood = negativeLogLikelihood
        self.rmse = rmse
        self.matchedControlRMSE = matchedControlRMSE
    }
}

/// A raw-source and perturbation-specific held-out score. Prepared corpus
/// receipts can be separate views of one raw source, so `source` deliberately
/// names the immutable raw count source rather than a receipt-bound corpus.
/// Qualification requires every such observed held-out pair to beat its own
/// exact matched-control baseline.
public struct VivoCellResponseSourceTargetEvaluation: Codable, Sendable, Equatable {
    public let source: VivoFingerprint
    public let targetID: String
    public let examples: Int
    public let observedFeatures: Int
    public let negativeLogLikelihood: Double
    public let rmse: Double
    public let matchedControlRMSE: Double

    public init(source: VivoFingerprint, targetID: String, examples: Int, observedFeatures: Int,
                negativeLogLikelihood: Double, rmse: Double, matchedControlRMSE: Double) {
        self.source = source
        self.targetID = targetID
        self.examples = examples
        self.observedFeatures = observedFeatures
        self.negativeLogLikelihood = negativeLogLikelihood
        self.rmse = rmse
        self.matchedControlRMSE = matchedControlRMSE
    }
}

/// A bounded, replayable software evaluation. It records the chosen rows so
/// reported metrics cannot drift with source ordering or a later sampler.
public struct VivoCellResponseEvaluation: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let partition: VivoCellResponsePartition
    public let maximumExamples: Int
    public let samplerVersion: UInt32
    public let seed: UInt64
    public let examples: Int
    public let observedFeatures: Int
    public let negativeLogLikelihood: Double
    public let rmse: Double
    /// RMSE of the exact declared control mean used as the prediction
    /// baseline, scored on exactly the rows and features in `selection`.
    public let matchedControlRMSE: Double
    public let selection: [VivoCellResponseEvaluationExample]
    /// One entry for every source in the bound composite, receipt-sorted.
    /// Aggregate metrics are the unweighted mean of these source scores.
    public let sourceMetrics: [VivoCellResponseSourceEvaluation]
    /// One entry for every raw-source and target pair selected for this
    /// evaluation. The entries are source then target sorted. Unlike
    /// `sourceMetrics`, these retain target-level regressions that a source
    /// aggregate could conceal.
    public let sourceTargetMetrics: [VivoCellResponseSourceTargetEvaluation]

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, format, partition, maximumExamples, samplerVersion, seed, examples,
             observedFeatures, negativeLogLikelihood, rmse, matchedControlRMSE, selection, sourceMetrics,
             sourceTargetMetrics
    }

    public init(schemaVersion: Int, format: String, partition: VivoCellResponsePartition,
                maximumExamples: Int, samplerVersion: UInt32, seed: UInt64, examples: Int,
                observedFeatures: Int, negativeLogLikelihood: Double, rmse: Double,
                matchedControlRMSE: Double, selection: [VivoCellResponseEvaluationExample],
                sourceMetrics: [VivoCellResponseSourceEvaluation] = [],
                sourceTargetMetrics: [VivoCellResponseSourceTargetEvaluation] = []) {
        self.schemaVersion = schemaVersion
        self.format = format
        self.partition = partition
        self.maximumExamples = maximumExamples
        self.samplerVersion = samplerVersion
        self.seed = seed
        self.examples = examples
        self.observedFeatures = observedFeatures
        self.negativeLogLikelihood = negativeLogLikelihood
        self.rmse = rmse
        self.matchedControlRMSE = matchedControlRMSE
        self.selection = selection
        self.sourceMetrics = sourceMetrics
        self.sourceTargetMetrics = sourceTargetMetrics
    }

    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "format", "partition", "maximumExamples", "samplerVersion", "seed", "examples",
            "observedFeatures", "negativeLogLikelihood", "rmse", "matchedControlRMSE", "selection", "sourceMetrics",
            "sourceTargetMetrics"
        ])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        format = try values.decode(String.self, forKey: .format)
        partition = try values.decode(VivoCellResponsePartition.self, forKey: .partition)
        maximumExamples = try values.decode(Int.self, forKey: .maximumExamples)
        samplerVersion = try values.decode(UInt32.self, forKey: .samplerVersion)
        seed = try values.decode(UInt64.self, forKey: .seed)
        examples = try values.decode(Int.self, forKey: .examples)
        observedFeatures = try values.decode(Int.self, forKey: .observedFeatures)
        negativeLogLikelihood = try values.decode(Double.self, forKey: .negativeLogLikelihood)
        rmse = try values.decode(Double.self, forKey: .rmse)
        matchedControlRMSE = try values.decode(Double.self, forKey: .matchedControlRMSE)
        selection = try values.decode([VivoCellResponseEvaluationExample].self, forKey: .selection)
        sourceMetrics = try values.decodeIfPresent([VivoCellResponseSourceEvaluation].self, forKey: .sourceMetrics) ?? []
        sourceTargetMetrics = try values.decodeIfPresent([VivoCellResponseSourceTargetEvaluation].self,
                                                          forKey: .sourceTargetMetrics) ?? []
    }

    public func encode(to encoder: Encoder) throws {
        var values = encoder.container(keyedBy: CodingKeys.self)
        try values.encode(schemaVersion, forKey: .schemaVersion)
        try values.encode(format, forKey: .format)
        try values.encode(partition, forKey: .partition)
        try values.encode(maximumExamples, forKey: .maximumExamples)
        try values.encode(samplerVersion, forKey: .samplerVersion)
        try values.encode(seed, forKey: .seed)
        try values.encode(examples, forKey: .examples)
        try values.encode(observedFeatures, forKey: .observedFeatures)
        try values.encode(negativeLogLikelihood, forKey: .negativeLogLikelihood)
        try values.encode(rmse, forKey: .rmse)
        try values.encode(matchedControlRMSE, forKey: .matchedControlRMSE)
        try values.encode(selection, forKey: .selection)
        if schemaVersion == 4 || schemaVersion == 5 {
            try values.encode(sourceMetrics, forKey: .sourceMetrics)
        }
        if schemaVersion == 5 {
            try values.encode(sourceTargetMetrics, forKey: .sourceTargetMetrics)
        }
    }
}

public struct VivoCellResponseEvaluationReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let format: String
    public let model: VivoFingerprint
    public let corpus: VivoFingerprint
    public let evaluation: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// Apple-native MLX implementation. The context branch is order-invariant over
/// declared control cells; it predicts a treatment residual over their mean,
/// so the model does not have to relearn the unperturbed expression profile.
final class VivoCellResponseMLXModel: Module {
    @ModuleInfo var cellProjection: Linear
    @ModuleInfo var targetEmbedding: Embedding
    @ModuleInfo var targetResponsePrior: Embedding
    @ModuleInfo var descriptorProjection: Linear
    @ModuleInfo var fusion: Linear
    @ModuleInfo var meanHead: Linear
    @ModuleInfo var varianceHead: Linear

    let featureCount: Int
    let hiddenWidth: Int

    init(featureCount: Int, targetCount: Int, descriptorCount: Int, architecture: VivoCellResponseArchitecture,
         targetResponsePriorValues: MLXArray? = nil) {
        self.featureCount = featureCount
        self.hiddenWidth = architecture.hiddenWidth
        cellProjection = Linear(featureCount, architecture.hiddenWidth)
        targetEmbedding = Embedding(embeddingCount: targetCount, dimensions: architecture.hiddenWidth)
        let prior = targetResponsePriorValues ?? MLXArray(
            Array(repeating: Float(0), count: targetCount * featureCount), [targetCount, featureCount])
        // Preserve the empirical training-only prior in the normal checkpoint
        // parameter tree, while excluding it from `valueAndGrad` and the
        // optimizer. Learned layers can only add a residual correction.
        let frozenPrior = Embedding(weight: prior)
        frozenPrior.freeze()
        targetResponsePrior = frozenPrior
        descriptorProjection = Linear(descriptorCount, architecture.hiddenWidth)
        fusion = Linear(architecture.hiddenWidth * 3, architecture.hiddenWidth)
        // Start from the training-only empirical response prior. The residual
        // decoder learns only corrections, so it cannot begin by degrading a
        // matched-control-plus-response prediction with random full-axis noise.
        meanHead = Linear(
            weight: MLXArray(Array(repeating: Float(0), count: featureCount * architecture.hiddenWidth),
                             [featureCount, architecture.hiddenWidth]),
            bias: MLXArray(Array(repeating: Float(0), count: featureCount)))
        varianceHead = Linear(architecture.hiddenWidth, featureCount)
    }

    func outputs(context: MLXArray, contextFeatureMask: MLXArray? = nil,
                 targetIDs: MLXArray, knownTargetMask: MLXArray,
                 descriptors: MLXArray) -> (mean: MLXArray, variance: MLXArray,
                                            baseline: MLXArray, delta: MLXArray) {
        let batch = context.shape[0]
        let contexts = context.shape[1]
        // A missing source-axis gene is structural absence, not a large
        // negative expression value. Mask before the shared cell encoder as
        // well as at the target likelihood, so heterogeneous source panels
        // cannot teach the context branch a missingness signature.
        let observedContext = contextFeatureMask.map { context * $0 } ?? context
        let projected = relu(cellProjection(observedContext.reshaped([batch * contexts, featureCount])))
            .reshaped([batch, contexts, hiddenWidth])
        let pooled = projected.mean(axis: 1)
        let baseline = observedContext.mean(axis: 1)
        let target = targetEmbedding(targetIDs) * knownTargetMask
        let prior = targetResponsePrior(targetIDs) * knownTargetMask
        let descriptor = relu(descriptorProjection(descriptors))
        let state = relu(fusion(concatenated([pooled, target, descriptor], axis: 1)))
        let delta = prior + meanHead(state)
        return (baseline + delta, softplus(varianceHead(state)) + 1e-4, baseline, delta)
    }
}

/// Versioned MLX response-learning surface. Artifacts remain source-bound by
/// the corpus receipt and use the shared native checkpoint codec for weights.
public enum VivoCellResponseLearning {
    public static let format = "numivivo-cell-response-learning/v2"
    /// Read-only compatibility identifier for pre-cohort technical artifacts.
    public static let legacyModelFormat = "numivivo-cell-response-model/v5"
    /// Read-only compatibility identifier for the one-corpus cohort format.
    public static let cohortV6ModelFormat = "numivivo-cell-response-model/v6"
    /// Current models bind a replayable, potentially multi-source cohort.
    public static let modelFormat = "numivivo-cell-response-model/v7"
    public static let legacyPredictionFormat = "numivivo-cell-response-prediction/v5"
    public static let predictionFormat = "numivivo-cell-response-prediction/v6"
    public static let legacyEvaluationFormat = "numivivo-cell-response-evaluation/v2"
    public static let cohortV3EvaluationFormat = "numivivo-cell-response-evaluation/v3"
    /// Read-only compatibility identifier for receipt-source macro metrics.
    public static let cohortV4EvaluationFormat = "numivivo-cell-response-evaluation/v4"
    /// Current evaluations retain raw-source and target-level held-out scores.
    public static let evaluationFormat = "numivivo-cell-response-evaluation/v5"
    static let samplerVersion: UInt32 = 3
    static let targetResponsePriorVersion: UInt32 = 1
    static let optimizerVersion: UInt32 = 1
    static let optimizerFormat = "sgd-momentum-0-mean-head-output-axis-v1"
}

private struct VivoCellResponseBatch {
    let context: MLXArray
    let contextFeatureMask: MLXArray
    let targetIDs: MLXArray
    let knownTargetMask: MLXArray
    let descriptors: MLXArray
    let targets: MLXArray
    let targetFeatureMask: MLXArray
    let targetValues: [Float]
    let contextRows: [Int]
}

private struct VivoCellResponseEvaluationTotals {
    let corpus: VivoFingerprint
    var examples = 0
    var observedFeatures = 0
    var negativeLogLikelihood = 0.0
    var squaredError = 0.0
    var baselineSquaredError = 0.0
}

private struct VivoCellResponseSourceTargetEvaluationTotals {
    let source: VivoFingerprint
    let targetID: String
    var examples = 0
    var observedFeatures = 0
    var negativeLogLikelihood = 0.0
    var squaredError = 0.0
    var baselineSquaredError = 0.0
}

private struct VivoCellResponseRawSourceTargetKey: Hashable {
    let source: String
    let targetID: String
}

/// Bind every receipt-qualified composite member to the raw source identity
/// retained by its cohort admission. A raw source can legitimately have more
/// than one prepared corpus, so callers must not use a corpus receipt as a
/// source-target qualification key.
private func vivoCellResponseRawSourceByCorpus(
    cohort: VivoCellResponseCohortAdmission,
    corpus: VivoCellResponseCompositeCorpusReader
) throws -> [String: VivoFingerprint] {
    let corpusReceipts = corpus.identity.sources.map(\.corpus.hex)
    guard cohort.sources.map(\.corpus.hex) == corpusReceipts else {
        throw VivoCellResponseLearningError.incompatible("evaluation cohort sources differ from corpus")
    }
    let mapping = Dictionary(uniqueKeysWithValues: cohort.sources.map { ($0.corpus.hex, $0.source) })
    guard mapping.count == corpus.sourceCount else {
        throw VivoCellResponseLearningError.invalid("evaluation raw source mapping")
    }
    return mapping
}

/// A source-qualified target vocabulary entry used only while assembling
/// training statistics. Global row offsets make the source component stable
/// across caller input order, while retaining each source's own biological
/// contexts and panel mask.
private struct VivoCellResponseSourceTargetKey: Hashable {
    let source: Int
    let target: Int
}

/// Partition response strata by their composite source. A matched-control
/// response must never cross source boundaries: otherwise an identical local
/// sample ID or a panel-specific missing feature could be mistaken for a
/// shared biological context.
private func vivoCellResponseSourceStratumIndices(
    corpus: VivoCellResponseCompositeCorpusReader,
    strata: [[VivoCellResponseTrainingExample]],
    requireEverySource: Bool = true
) throws -> [[Int]] {
    guard !strata.isEmpty else {
        throw VivoCellResponseLearningError.invalid("empty source response strata")
    }
    var result = Array(repeating: [Int](), count: corpus.sourceCount)
    for (stratumIndex, candidates) in strata.enumerated() {
        guard let first = candidates.first, !first.contextRows.isEmpty else {
            throw VivoCellResponseLearningError.invalid("source response stratum")
        }
        let source = try corpus.sourceIndex(forGlobalRow: first.targetRow)
        for candidate in candidates {
            guard try corpus.sourceIndex(forGlobalRow: candidate.targetRow) == source,
                  !candidate.contextRows.isEmpty else {
                throw VivoCellResponseLearningError.invalid("source response target row")
            }
            for contextRow in candidate.contextRows {
                guard try corpus.sourceIndex(forGlobalRow: contextRow) == source else {
                    throw VivoCellResponseLearningError.invalid("source response control row")
                }
            }
        }
        result[source].append(stratumIndex)
    }
    if requireEverySource, !result.allSatisfy({ !$0.isEmpty }) {
        throw VivoCellResponseLearningError.invalid("every bound source needs a response stratum")
    }
    return result
}

/// Ephemeral, corpus-scoped normalized control rows used only while fitting.
/// It is rebuilt from the verified source store for each train/resume call and
/// is never published as model state or evidence.
private struct VivoCellResponseTrainingControlCache {
    let featureCount: Int
    let values: [Float]
    let offsets: [Int: Int]
    let rankedContexts: [[Int]]

    init(corpus: VivoCellResponseCompositeCorpusReader, rows: [Int], rankedContexts: [[Int]]) throws {
        featureCount = corpus.featureCount
        self.rankedContexts = rankedContexts
        guard !rows.isEmpty else {
            throw VivoCellResponseLearningError.invalid("training control cache is empty")
        }
        let valueCount = rows.count.multipliedReportingOverflow(by: featureCount)
        guard !valueCount.overflow, valueCount.partialValue <= 268_435_456 else {
            throw VivoCellResponseLearningError.limit("training control cache exceeds 1 GiB")
        }
        var storage: [Float] = []
        storage.reserveCapacity(valueCount.partialValue)
        var rowOffsets: [Int: Int] = [:]
        rowOffsets.reserveCapacity(rows.count)
        for (offset, row) in rows.enumerated() {
            if offset % 8 == 0 { try Task.checkCancellation() }
            let normalized = try corpus.normalizedRow(row)
            guard normalized.count == featureCount else {
                throw VivoCellResponseLearningError.invalid("training control cache row width")
            }
            rowOffsets[row] = storage.count
            storage.append(contentsOf: normalized)
        }
        guard storage.count == valueCount.partialValue else {
            throw VivoCellResponseLearningError.invalid("training control cache size")
        }
        values = storage
        offsets = rowOffsets
    }

    func rows(count: Int, seed: UInt64, step: UInt64, lane: UInt64, stratum: Int) throws -> [Int] {
        guard rankedContexts.indices.contains(stratum) else {
            throw VivoCellResponseLearningError.invalid("training control cache stratum")
        }
        let candidates = rankedContexts[stratum]
        guard count > 0, !candidates.isEmpty else {
            throw VivoCellResponseLearningError.invalid("training control cache context")
        }
        let start = try VivoCellResponseLearning.sampledIndex(
            count: candidates.count, seed: seed, step: step, lane: lane,
            salt: 0x721a5be9 ^ UInt64(stratum))
        return (0..<count).map { candidates[(start + $0) % candidates.count] }
    }

    func appendNormalizedRow(_ row: Int, to destination: inout [Float]) throws {
        guard let offset = offsets[row], offset >= 0,
              offset <= values.count - featureCount else {
            throw VivoCellResponseLearningError.invalid("training control cache row")
        }
        destination.append(contentsOf: values[offset..<(offset + featureCount)])
    }
}

private struct VivoCellResponseLoadedArtifact {
    let model: VivoCellResponseMLXModel
    let state: VivoCellResponseModelState
    let receipt: VivoCellResponseModelReceipt
    let trainingPlan: VivoCellResponseTrainingPlan
    let trainingPlanBytes: Data
    let cohortAdmission: VivoCellResponseCohortAdmission?
    let cohortAdmissionBytes: Data?
}

private enum VivoCellResponseArtifactIO {
    static let checkpointName = "checkpoint.nvckpt"
    static let planName = "training-plan.json"
    static let cohortAdmissionName = "cohort-admission.json"
    static let receiptName = "receipt.json"
    static let weightsSection = "weights.safetensors"
    static let stateSection = "state.json"

    static func requireNew(_ destination: URL) throws {
        let manager = FileManager.default
        guard destination.isFileURL, destination.path.utf8.count <= 8_192,
              !destination.path.contains("\0"), !manager.fileExists(atPath: destination.path) else {
            throw VivoCellResponseLearningError.invalid("output directory already exists or path is invalid")
        }
        var directory = ObjCBool(false)
        guard manager.fileExists(atPath: destination.deletingLastPathComponent().path, isDirectory: &directory),
              directory.boolValue else {
            throw VivoCellResponseLearningError.invalid("output parent directory")
        }
    }

    static func staging(_ parent: URL, prefix: String) throws -> URL {
        let directory = parent.appendingPathComponent("." + prefix + "-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false,
                                                attributes: [.posixPermissions: 0o700])
        return directory
    }

    static func read(_ root: VivoRootedDocumentStore, _ name: String, maximumBytes: Int) throws -> Data {
        do {
            return try root.readDocument(name, maximumBytes: maximumBytes)
        } catch {
            throw VivoCellResponseLearningError.invalid("artifact document \(name)")
        }
    }

    static func write(_ data: Data, to root: URL, name: String, maximumBytes: Int) throws {
        guard data.count <= maximumBytes else { throw VivoCellResponseLearningError.limit("artifact document \(name)") }
        try data.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
    }

    static func corpusFingerprint(_ corpus: VivoCellResponseCompositeCorpusReader) throws -> VivoFingerprint {
        corpus.fingerprint
    }

    static func modelFingerprint(_ receipt: VivoCellResponseModelReceipt) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
    }

    static func format(for state: VivoCellResponseModelState) -> String {
        guard state.cohort != nil else { return VivoCellResponseLearning.legacyModelFormat }
        return state.schemaVersion == 6 ? VivoCellResponseLearning.modelFormat : VivoCellResponseLearning.cohortV6ModelFormat
    }

    static func checkpoint(model: VivoCellResponseMLXModel, state: VivoCellResponseModelState,
                           parentCheckpoint: VivoFingerprint?) throws -> Data {
        let format = format(for: state)
        let arrays = Dictionary(uniqueKeysWithValues: model.parameters().flattened())
        let metadata = [
            "format": format,
            "optimizer": VivoCellResponseLearning.optimizerFormat,
            "sampler": "splitmix64-v1",
            "targetResponsePrior": "balanced-training-mean-v1",
            "featureCount": String(state.featureAxis.featureIDs.count),
            "step": String(state.step)
        ]
        let weights = try MLX.saveToData(arrays: arrays, metadata: metadata)
        let request = VivoCheckpointBuildRequest(
            runtime: .cellResponseMLX,
            artifactFingerprint: state.corpus.hex,
            sourceFingerprint: state.corpus.hex,
            experimentFingerprint: state.trainingPlan.hex,
            parentCheckpointFingerprint: parentCheckpoint?.hex,
            stepIndex: state.step,
            logicalTime: Double(state.step),
            randomStreamVersion: state.samplerVersion,
            sections: [
                .init(id: weightsSection, encoding: .rawBytes, elementCount: UInt64(weights.count),
                      elementStride: 1, data: weights),
                try .canonicalJSON(id: stateSection, value: state)
            ],
            metadata: [
                "format": format,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1"
            ])
        return try VivoCheckpointCodec.encode(request)
    }

    static func validateState(_ state: VivoCellResponseModelState,
                              trainingPlan: VivoCellResponseTrainingPlan) throws {
        let qualified = state.cohort != nil
        guard (qualified ? (state.schemaVersion == 5 || state.schemaVersion == 6) : state.schemaVersion == 4),
              state.samplerVersion == VivoCellResponseLearning.samplerVersion,
              state.targetResponsePriorVersion == VivoCellResponseLearning.targetResponsePriorVersion,
              state.optimizerVersion == VivoCellResponseLearning.optimizerVersion,
              state.featureAxis.featureIDs.count > 0, state.featureAxis.featureIDs.count <= 200_000,
              state.targetIDs.count > 0, Set(state.targetIDs).count == state.targetIDs.count,
              state.targetIDs.allSatisfy(vivoOmicsID), !state.trainedTargetBindings.isEmpty,
              state.trainedTargetBindings.count <= state.targetIDs.count,
              Set(state.trainedTargetBindings.map(\.id)).count == state.trainedTargetBindings.count,
              state.trainedTargetBindings.allSatisfy({ state.targetIDs.contains($0.id) }),
              state.descriptorCount > 0, state.learningRate.isFinite, state.learningRate > 0,
              state.weightDecay.isFinite, state.weightDecay >= 0,
              state.latestMetrics.step == state.step,
              state.latestMetrics.trainNegativeLogLikelihood.isFinite,
              state.latestMetrics.validationNegativeLogLikelihood?.isFinite ?? true,
              state.latestMetrics.validationRMSE?.isFinite ?? true,
              state.latestMetrics.validationRMSE.map({ $0 >= 0 }) ?? true,
              state.latestMetrics.validationMatchedControlRMSE?.isFinite ?? true,
              state.latestMetrics.validationMatchedControlRMSE.map({ $0 >= 0 }) ?? true,
              state.architecture == trainingPlan.architecture,
              state.seed == trainingPlan.seed,
              state.learningRate == trainingPlan.learningRate,
              state.weightDecay == trainingPlan.weightDecay else {
            throw VivoCellResponseLearningError.invalid("checkpoint state")
        }
        try state.featureAxis.validate()
        try state.architecture.validate(featureCount: state.featureAxis.featureIDs.count,
                                        targetCount: state.targetIDs.count,
                                        descriptorCount: state.descriptorCount)
    }

    static func load(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseLoadedArtifact {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("model artifact directory")
        }
        let receiptBytes = try read(root, receiptName, maximumBytes: 1_048_576)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponseModelReceipt.self, from: receiptBytes)
        let qualified = receipt.format == VivoCellResponseLearning.modelFormat ||
            receipt.format == VivoCellResponseLearning.cohortV6ModelFormat
        guard ((receipt.format == VivoCellResponseLearning.modelFormat && receipt.schemaVersion == 3 && receipt.cohort != nil) ||
               (receipt.format == VivoCellResponseLearning.cohortV6ModelFormat && receipt.schemaVersion == 2 && receipt.cohort != nil) ||
               (!qualified && receipt.schemaVersion == 1 &&
                receipt.format == VivoCellResponseLearning.legacyModelFormat && receipt.cohort == nil)),
              receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes else {
            throw VivoCellResponseLearningError.invalid("model receipt")
        }
        let trainingPlanBytes = try read(root, planName, maximumBytes: 4_194_304)
        let trainingPlan = try VivoCanonicalJSON.decode(VivoCellResponseTrainingPlan.self, from: trainingPlanBytes)
        guard try VivoCanonicalJSON.encode(trainingPlan) == trainingPlanBytes else {
            throw VivoCellResponseLearningError.invalid("noncanonical training plan")
        }
        try trainingPlan.validateStatic()
        var cohortAdmission: VivoCellResponseCohortAdmission?
        var cohortAdmissionBytes: Data?
        if qualified {
            guard let cohortFingerprint = receipt.cohort else {
            throw VivoCellResponseLearningError.invalid("cohort-bound model receipt")
            }
            let bytes = try read(root, cohortAdmissionName, maximumBytes: 268_435_456)
            let admission = try VivoCanonicalJSON.decode(VivoCellResponseCohortAdmission.self, from: bytes)
            guard try VivoCanonicalJSON.encode(admission) == bytes,
                  try VivoCanonicalJSON.fingerprint(bytes) == cohortFingerprint else {
            throw VivoCellResponseLearningError.invalid("cohort-bound model artifact")
            }
            try VivoCellResponseCohort.validate(admission)
            cohortAdmission = admission
            cohortAdmissionBytes = bytes
        }
        let checkpointBytes = try read(root, checkpointName, maximumBytes: 1_073_741_824)
        guard try VivoCanonicalJSON.fingerprint(trainingPlanBytes) == receipt.trainingPlan,
              try VivoCanonicalJSON.fingerprint(checkpointBytes) == receipt.checkpoint else {
            throw VivoCellResponseLearningError.invalid("model artifact fingerprints")
        }
        let decoded = try VivoCheckpointCodec.decode(checkpointBytes)
        let stateBytes = try decoded.section(id: stateSection)
        let weightsBytes = try decoded.section(id: weightsSection)
        let state = try VivoCanonicalJSON.decode(VivoCellResponseModelState.self, from: stateBytes)
        let sections = decoded.manifest.sections
        let weightsInfo = sections.first(where: { $0.id == weightsSection })
        let stateInfo = sections.first(where: { $0.id == stateSection })
        guard decoded.manifest.runtime == .cellResponseMLX,
              decoded.manifest.artifactFingerprint == receipt.corpus.hex,
              decoded.manifest.sourceFingerprint == receipt.corpus.hex,
              decoded.manifest.experimentFingerprint == receipt.trainingPlan.hex,
              decoded.manifest.logicalTime == Double(state.step),
              decoded.manifest.randomStreamVersion == VivoCellResponseLearning.samplerVersion,
              sections.count == 2, Set(sections.map(\.id)) == Set([weightsSection, stateSection]),
              let weightsInfo, weightsInfo.encoding == .rawBytes,
              weightsInfo.length == UInt64(weightsBytes.count), weightsInfo.elementCount == UInt64(weightsBytes.count),
              weightsInfo.elementStride == 1,
              let stateInfo, stateInfo.encoding == .canonicalJSON,
              stateInfo.length == UInt64(stateBytes.count), stateInfo.elementCount == 1,
              stateInfo.elementStride == UInt32(clamping: stateBytes.count),
              decoded.manifest.metadata == [
                "format": receipt.format,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1"
              ] else {
            throw VivoCellResponseLearningError.invalid("cell-response checkpoint manifest")
        }
        guard try VivoCanonicalJSON.encode(state) == stateBytes,
              state.corpus == receipt.corpus, state.trainingPlan == receipt.trainingPlan,
              state.cohort == receipt.cohort,
              format(for: state) == receipt.format,
              state.step == decoded.manifest.stepIndex else {
            throw VivoCellResponseLearningError.invalid("cell-response checkpoint state binding")
        }
        try validateState(state, trainingPlan: trainingPlan)
        let (arrays, metadata) = try MLX.loadArraysAndMetadata(data: weightsBytes)
        guard metadata["format"] == receipt.format,
              metadata["optimizer"] == VivoCellResponseLearning.optimizerFormat,
              metadata["sampler"] == "splitmix64-v1",
              metadata["targetResponsePrior"] == "balanced-training-mean-v1",
              metadata["featureCount"] == String(state.featureAxis.featureIDs.count),
              metadata["step"] == String(state.step) else {
            throw VivoCellResponseLearningError.invalid("cell-response weight metadata")
        }
        let model = VivoCellResponseMLXModel(featureCount: state.featureAxis.featureIDs.count,
                                             targetCount: state.targetIDs.count,
                                             descriptorCount: state.descriptorCount,
                                             architecture: state.architecture)
        try model.update(parameters: ModuleParameters.unflattened(arrays), verify: .all)
        eval(model)
        for (_, parameter) in model.parameters().flattened() {
            guard parameter.asArray(Float.self).allSatisfy(\.isFinite) else {
                throw VivoCellResponseLearningError.invalid("nonfinite model parameter")
            }
        }
        return .init(model: model, state: state, receipt: receipt, trainingPlan: trainingPlan,
                     trainingPlanBytes: trainingPlanBytes, cohortAdmission: cohortAdmission,
                     cohortAdmissionBytes: cohortAdmissionBytes)
    }

    static func publish(model: VivoCellResponseMLXModel, state: VivoCellResponseModelState,
                        trainingPlanBytes: Data, implementation: VivoFingerprint,
                        parentCheckpoint: VivoFingerprint?,
                        cohortAdmission: VivoCellResponseCohortAdmission?,
                        cohortAdmissionBytes: Data?,
                        to destination: URL) throws -> VivoCellResponseModelReceipt {
        try requireNew(destination)
        let trainingPlan = try VivoCanonicalJSON.decode(VivoCellResponseTrainingPlan.self, from: trainingPlanBytes)
        guard try VivoCanonicalJSON.encode(trainingPlan) == trainingPlanBytes else {
            throw VivoCellResponseLearningError.invalid("noncanonical training plan")
        }
        try trainingPlan.validateStatic()
        let qualified = state.cohort != nil
        guard qualified == (cohortAdmission != nil), qualified == (cohortAdmissionBytes != nil) else {
            throw VivoCellResponseLearningError.invalid("model cohort publication binding")
        }
        if let cohortAdmission, let cohortAdmissionBytes, let cohort = state.cohort {
            try VivoCellResponseCohort.validate(cohortAdmission)
            guard try VivoCanonicalJSON.encode(cohortAdmission) == cohortAdmissionBytes,
                  try VivoCanonicalJSON.fingerprint(cohortAdmissionBytes) == cohort else {
                throw VivoCellResponseLearningError.invalid("model cohort publication artifact")
            }
        }
        let temporary = try staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-model")
        defer { try? FileManager.default.removeItem(at: temporary) }
        let checkpoint = try checkpoint(model: model, state: state, parentCheckpoint: parentCheckpoint)
        let receipt = try VivoCellResponseModelReceipt(
            schemaVersion: qualified ? (state.schemaVersion == 6 ? 3 : 2) : 1, format: format(for: state),
            checkpoint: VivoCanonicalJSON.fingerprint(checkpoint),
            trainingPlan: VivoCanonicalJSON.fingerprint(trainingPlanBytes), corpus: state.corpus,
            cohort: state.cohort,
            implementation: implementation)
        try write(trainingPlanBytes, to: temporary, name: planName, maximumBytes: 4_194_304)
        if let cohortAdmissionBytes {
            try write(cohortAdmissionBytes, to: temporary, name: cohortAdmissionName, maximumBytes: 268_435_456)
        }
        try write(checkpoint, to: temporary, name: checkpointName, maximumBytes: 1_073_741_824)
        try write(VivoCanonicalJSON.encode(receipt), to: temporary, name: receiptName, maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }
}

extension VivoCellResponseLearning {
    private static func withGPUExecution<R>(_ body: () throws -> R) rethrows -> R {
        return try Device.withDefaultDevice(.gpu) {
            try Stream.withNewDefaultStream(device: .gpu) {
                try body()
            }
        }
    }

    private static func mix(seed: UInt64, step: UInt64, lane: UInt64, salt: UInt64) -> UInt64 {
        var value = seed &+ 0x9e3779b97f4a7c15 &* (step &+ 1)
        value = value &+ 0xbf58476d1ce4e5b9 &* (lane &+ 1)
        value = value &+ salt
        value = (value ^ (value >> 30)) &* 0xbf58476d1ce4e5b9
        value = (value ^ (value >> 27)) &* 0x94d049bb133111eb
        return value ^ (value >> 31)
    }

    fileprivate static func sampledIndex(count: Int, seed: UInt64, step: UInt64, lane: UInt64, salt: UInt64) throws -> Int {
        guard count > 0 else { throw VivoCellResponseLearningError.invalid("empty sample population") }
        return Int(mix(seed: seed, step: step, lane: lane, salt: salt) % UInt64(count))
    }

    /// Select a reproducible subset without replacement. Sorting the selected
    /// indices makes traversal independent of hash-table iteration order.
    private static func sampledIndicesWithoutReplacement(count: Int, take: Int, seed: UInt64,
                                                         salt: UInt64) throws -> [Int] {
        guard count > 0, take > 0, take <= count else {
            throw VivoCellResponseLearningError.invalid("evaluation sample population")
        }
        var selected = Set<Int>()
        selected.reserveCapacity(take)
        for position in 0..<take {
            let upper = count - take + position
            let candidate = Int(mix(seed: seed, step: UInt64(position), lane: 0, salt: salt) % UInt64(upper + 1))
            if !selected.insert(candidate).inserted {
                selected.insert(upper)
            }
        }
        guard selected.count == take else { throw VivoCellResponseLearningError.invalid("evaluation subset") }
        return selected.sorted()
    }

    /// Group each perturbed source sample with its declared matched-control
    /// context. Sampling a group before a cell avoids letting high-cell-count
    /// guides dominate training or held-out metrics.
    private static func stratifiedExamples(_ examples: [VivoCellResponseTrainingExample]) throws -> [[VivoCellResponseTrainingExample]] {
        let groups = Dictionary(grouping: examples, by: \.stratumIndex)
        guard !groups.isEmpty else { throw VivoCellResponseLearningError.invalid("empty response strata") }
        return try groups.keys.sorted().map { index in
            guard let rows = groups[index], !rows.isEmpty else {
                throw VivoCellResponseLearningError.invalid("empty response stratum")
            }
            return rows.sorted { $0.targetRow < $1.targetRow }
        }
    }

    /// Choose a source-order-independent control context shared by train,
    /// evaluation, and inference. A fixed hash ranks the immutable row IDs;
    /// each row is consumed once before cycling only when the model width
    /// exceeds the declared control population.
    private static func rankedContextRows(_ rows: [Int]) throws -> [Int] {
        guard !rows.isEmpty else {
            throw VivoCellResponseLearningError.invalid("empty deterministic control context")
        }
        return rows.sorted { left, right in
            let leftRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(left), lane: 0, salt: 0x6841f29d)
            let rightRank = mix(seed: 0x0d4e3c2b1a987654, step: UInt64(right), lane: 0, salt: 0x6841f29d)
            return leftRank == rightRank ? left < right : leftRank < rightRank
        }
    }

    private static func deterministicContextRows(_ rows: [Int], count: Int) throws -> [Int] {
        let candidates = try rankedContextRows(rows)
        guard count > 0 else {
            throw VivoCellResponseLearningError.invalid("empty deterministic control context")
        }
        return (0..<count).map { candidates[$0 % candidates.count] }
    }

    /// A replay-bound cyclic window through a deterministic control ordering.
    /// Unlike the inference/evaluation context, the starting point changes for
    /// each optimizer lane so training eventually sees the declared control
    /// population instead of repeatedly fitting one fixed subset.
    private static func trainingContextRows(_ rows: [Int], count: Int, seed: UInt64,
                                            step: UInt64, lane: UInt64, stratum: Int) throws -> [Int] {
        let candidates = try rankedContextRows(rows)
        guard count > 0, stratum >= 0 else {
            throw VivoCellResponseLearningError.invalid("empty training control context")
        }
        let start = try sampledIndex(count: candidates.count, seed: seed, step: step, lane: lane,
                                     salt: 0x721a5be9 ^ UInt64(stratum))
        return (0..<count).map { candidates[(start + $0) % candidates.count] }
    }

    private static func trainingControlCache(corpus: VivoCellResponseCompositeCorpusReader,
                                             strata: [[VivoCellResponseTrainingExample]]) throws -> VivoCellResponseTrainingControlCache {
        let maximumContextReferences = 1_048_576
        guard !strata.isEmpty, strata.count <= maximumContextReferences else {
            throw VivoCellResponseLearningError.limit("training control cache strata")
        }
        var rankedContexts: [[Int]] = []
        rankedContexts.reserveCapacity(strata.count)
        var totalContextReferences = 0
        var uniqueRows = Set<Int>()
        for candidates in strata {
            guard let first = candidates.first,
                  candidates.allSatisfy({ $0.contextRows == first.contextRows }) else {
                throw VivoCellResponseLearningError.invalid("training control cache strata")
            }
            let nextCount = totalContextReferences.addingReportingOverflow(first.contextRows.count)
            guard !nextCount.overflow, nextCount.partialValue <= maximumContextReferences else {
                throw VivoCellResponseLearningError.limit("training control cache context references")
            }
            totalContextReferences = nextCount.partialValue
            let ranked = try rankedContextRows(first.contextRows)
            rankedContexts.append(ranked)
            uniqueRows.formUnion(ranked)
        }
        return try .init(corpus: corpus, rows: uniqueRows.sorted(), rankedContexts: rankedContexts)
    }

    /// The empirical perturbation response used to initialize known targets.
    /// For every training source×target×context stratum, it computes the
    /// difference between the mean treated log-CPM vector and the mean of its
    /// declared matched controls. Contexts receive equal weight within a
    /// source, then sources receive equal weight for each observed feature.
    /// Validation and test rows are never opened here.
    private static func trainingOnlyTargetResponsePrior(corpus: VivoCellResponseCompositeCorpusReader,
                                                        trainedTargetIndices: Set<Int>) throws -> MLXArray {
        let strata = try stratifiedExamples(try corpus.examples(in: .training))
        let features = corpus.featureCount
        var prior = [Double](repeating: 0, count: corpus.targetCount * features)
        var sourceObservations = [UInt32](repeating: 0, count: corpus.targetCount * features)
        var controlMeans: [String: [Double]] = [:]
        var grouped: [VivoCellResponseSourceTargetKey: [[VivoCellResponseTrainingExample]]] = [:]

        func mean(_ rows: [Int]) throws -> [Double] {
            guard !rows.isEmpty else { throw VivoCellResponseLearningError.invalid("empty response-prior rows") }
            var result = [Double](repeating: 0, count: features)
            for (rowOffset, row) in rows.enumerated() {
                if rowOffset % 32 == 0 { try Task.checkCancellation() }
                let values = try corpus.normalizedRow(row)
                guard values.count == features else { throw VivoCellResponseLearningError.invalid("response-prior row width") }
                for feature in 0..<features {
                    result[feature] += Double(values[feature])
                }
            }
            let divisor = Double(rows.count)
            return result.map { $0 / divisor }
        }

        for candidates in strata {
            try Task.checkCancellation()
            guard let first = candidates.first,
                  trainedTargetIndices.contains(first.targetIndex),
                  candidates.allSatisfy({ $0.targetIndex == first.targetIndex && $0.contextRows == first.contextRows }) else {
                throw VivoCellResponseLearningError.invalid("response-prior stratum")
            }
            let source = try corpus.sourceIndex(forGlobalRow: first.targetRow)
            let targetMask = try corpus.featureMask(forGlobalRow: first.targetRow)
            guard targetMask.count == features, targetMask.contains(where: { $0 > 0 }) else {
                throw VivoCellResponseLearningError.invalid("response-prior target feature mask")
            }
            for candidate in candidates {
                guard try corpus.sourceIndex(forGlobalRow: candidate.targetRow) == source,
                      try corpus.featureMask(forGlobalRow: candidate.targetRow) == targetMask else {
                    throw VivoCellResponseLearningError.invalid("response-prior target source or feature mask")
                }
            }
            for row in first.contextRows {
                guard try corpus.sourceIndex(forGlobalRow: row) == source,
                      try corpus.featureMask(forGlobalRow: row) == targetMask else {
                    throw VivoCellResponseLearningError.invalid("response-prior control source or feature mask")
                }
            }
            grouped[.init(source: source, target: first.targetIndex), default: []].append(candidates)
        }

        for key in grouped.keys.sorted(by: { left, right in
            left.source == right.source ? left.target < right.target : left.source < right.source
        }) {
            try Task.checkCancellation()
            guard let sourceStrata = grouped[key], !sourceStrata.isEmpty else {
                throw VivoCellResponseLearningError.invalid("response-prior source target")
            }
            var sourceSum = [Double](repeating: 0, count: features)
            var sourceContexts = [UInt32](repeating: 0, count: features)
            for candidates in sourceStrata {
                guard let first = candidates.first else {
                    throw VivoCellResponseLearningError.invalid("response-prior source stratum")
                }
                let targetMask = try corpus.featureMask(forGlobalRow: first.targetRow)
                let contextKey = first.contextRows.map(String.init).joined(separator: ",")
                let controlMean: [Double]
                if let existing = controlMeans[contextKey] {
                    controlMean = existing
                } else {
                    let computed = try mean(first.contextRows)
                    controlMeans[contextKey] = computed
                    controlMean = computed
                }
                let treatedMean = try mean(candidates.map(\.targetRow))
                for feature in 0..<features where targetMask[feature] > 0 {
                    guard sourceContexts[feature] < UInt32.max else {
                        throw VivoCellResponseLearningError.limit("response-prior source contexts")
                    }
                    sourceSum[feature] += treatedMean[feature] - controlMean[feature]
                    sourceContexts[feature] += 1
                }
            }
            let offset = key.target * features
            for feature in 0..<features where sourceContexts[feature] > 0 {
                guard sourceObservations[offset + feature] < UInt32.max else {
                    throw VivoCellResponseLearningError.limit("response-prior source observations")
                }
                prior[offset + feature] += sourceSum[feature] / Double(sourceContexts[feature])
                sourceObservations[offset + feature] += 1
            }
        }
        for target in trainedTargetIndices {
            let offset = target * features
            guard sourceObservations[offset..<(offset + features)].contains(where: { $0 > 0 }) else {
                throw VivoCellResponseLearningError.invalid("response prior has an unobserved target")
            }
            for feature in 0..<features where sourceObservations[offset + feature] > 0 {
                prior[offset + feature] /= Double(sourceObservations[offset + feature])
            }
        }
        guard prior.allSatisfy(\.isFinite) else { throw VivoCellResponseLearningError.invalid("nonfinite response prior") }
        return MLXArray(prior.map(Float.init), [corpus.targetCount, features])
    }

    private static func trainedTargetIndices(corpus: VivoCellResponseCompositeCorpusReader,
                                             bindings: [VivoCellResponseTrainedTargetBinding]) throws -> Set<Int> {
        guard !bindings.isEmpty, bindings.count <= corpus.targetCount,
              Set(bindings.map(\.id)).count == bindings.count else {
            throw VivoCellResponseLearningError.incompatible("trained target bindings")
        }
        let vocabulary = Dictionary(uniqueKeysWithValues: corpus.targets.enumerated().map { ($0.element.id, $0.offset) })
        var indices = Set<Int>()
        for binding in bindings {
            guard let index = vocabulary[binding.id],
                  try corpus.targets[index].fingerprint() == binding.descriptorFingerprint else {
                throw VivoCellResponseLearningError.incompatible("trained target binding differs from corpus")
            }
            indices.insert(index)
        }
        return indices
    }

    private static func validateCheckpointCorpus(_ state: VivoCellResponseModelState,
                                                 corpus: VivoCellResponseCompositeCorpusReader) throws -> Set<Int> {
        guard state.featureAxis == corpus.featureAxis,
              state.targetIDs == corpus.targets.map(\.id),
              state.descriptorCount == corpus.descriptorCount else {
            throw VivoCellResponseLearningError.incompatible("checkpoint corpus schema differs")
        }
        let actualBindings = try corpus.trainingTargetBindings()
        guard state.trainedTargetBindings == actualBindings else {
            throw VivoCellResponseLearningError.incompatible("checkpoint trained targets differ from corpus")
        }
        return try trainedTargetIndices(corpus: corpus, bindings: state.trainedTargetBindings)
    }

    private static func loss(model: VivoCellResponseMLXModel, context: MLXArray,
                             contextFeatureMask: MLXArray, targetIDs: MLXArray,
                             knownTargetMask: MLXArray, descriptors: MLXArray, targets: MLXArray,
                             targetFeatureMask: MLXArray) -> MLXArray {
        let output = model.outputs(context: context, contextFeatureMask: contextFeatureMask,
                                   targetIDs: targetIDs, knownTargetMask: knownTargetMask,
                                   descriptors: descriptors)
        // Learn the perturbation response after subtracting the explicit
        // matched-control baseline. Algebraically this preserves the absolute
        // Gaussian likelihood while making the learned head a delta decoder.
        let response = targets - output.baseline
        let residual = output.delta - response
        let gaussianConstant = Float(log(2.0 * Double.pi))
        let terms = (residual * residual / output.variance + log(output.variance) + gaussianConstant) * 0.5
        return sourceBalancedMaskedMean(terms, mask: targetFeatureMask)
    }

    /// Give every scheduled source lane equal likelihood weight even when its
    /// source panel has fewer observed genes. `trainingBatch` rejects an empty
    /// feature mask for every lane, so a zero denominator is never admitted.
    static func sourceBalancedMaskedMean(_ terms: MLXArray, mask: MLXArray) -> MLXArray {
        let perLaneTerms = (terms * mask).sum(axis: 1)
        let observedPerLane = mask.sum(axis: 1)
        return (perLaneTerms / observedPerLane).mean()
    }

    private static func trainingBatch(corpus: VivoCellResponseCompositeCorpusReader,
                                      strata: [[VivoCellResponseTrainingExample]],
                                      sourceStrata: [[Int]],
                                      plan: VivoCellResponseTrainingPlan, trainedTargetIndices: Set<Int>,
                                      step: UInt64,
                                      controlCache: VivoCellResponseTrainingControlCache?) throws -> VivoCellResponseBatch {
        let features = corpus.featureCount
        let contexts = plan.architecture.contextCells
        let batches = plan.batchSize
        guard sourceStrata.count == corpus.sourceCount,
              sourceStrata.allSatisfy({ !$0.isEmpty }),
              corpus.featureMask.count == features, corpus.featureMask.contains(where: { $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("corpus has no measured prediction features")
        }
        var contextValues: [Float] = []
        var contextMasks: [Float] = []
        var targetValues: [Float] = []
        var targetMasks: [Float] = []
        var targetIDs: [Int32] = []
        var known: [Float] = []
        var descriptors: [Float] = []
        contextValues.reserveCapacity(batches * contexts * features)
        contextMasks.reserveCapacity(batches * contexts * features)
        targetValues.reserveCapacity(batches * features)
        targetMasks.reserveCapacity(batches * features)
        targetIDs.reserveCapacity(batches); known.reserveCapacity(batches); descriptors.reserveCapacity(batches * corpus.descriptorCount)
        for lane in 0..<batches {
            // Cycle sources before response strata. Each source therefore gets
            // equal optimizer exposure even when one study contributes many
            // more guide/context strata than another.
            let completed = step.subtractingReportingOverflow(1)
            let scaled = completed.partialValue.multipliedReportingOverflow(by: UInt64(batches))
            let scheduled = scaled.partialValue.addingReportingOverflow(UInt64(lane))
            guard !completed.overflow, !scaled.overflow, !scheduled.overflow else {
                throw VivoCellResponseLearningError.limit("training schedule")
            }
            let sourceIndex = Int(scheduled.partialValue % UInt64(sourceStrata.count))
            let sourceSchedule = scheduled.partialValue / UInt64(sourceStrata.count)
            let sourceCandidates = sourceStrata[sourceIndex]
            let stratumIndex = sourceCandidates[Int(sourceSchedule % UInt64(sourceCandidates.count))]
            let candidates = strata[stratumIndex]
            let exampleIndex = try sampledIndex(count: candidates.count, seed: plan.seed, step: step,
                                                lane: UInt64(lane), salt: 0x4c7912ae)
            let example = candidates[exampleIndex]
            guard trainedTargetIndices.contains(example.targetIndex),
                  let targetID = Int32(exactly: example.targetIndex), example.contextRows.count > 0 else {
                throw VivoCellResponseLearningError.invalid("training target or context rows")
            }
            let targetMask = try corpus.featureMask(forGlobalRow: example.targetRow)
            guard targetMask.count == features, targetMask.contains(where: { $0 > 0 }) else {
                throw VivoCellResponseLearningError.invalid("training target feature mask")
            }
            let contextRows: [Int]
            if let controlCache {
                contextRows = try controlCache.rows(count: contexts, seed: plan.seed, step: step,
                                                     lane: UInt64(lane), stratum: stratumIndex)
            } else {
                contextRows = try trainingContextRows(example.contextRows, count: contexts, seed: plan.seed,
                                                       step: step, lane: UInt64(lane), stratum: stratumIndex)
            }
            for row in contextRows {
                guard try corpus.featureMask(forGlobalRow: row) == targetMask else {
                    throw VivoCellResponseLearningError.invalid("training context source feature mask")
                }
                if let controlCache {
                    try controlCache.appendNormalizedRow(row, to: &contextValues)
                } else {
                    contextValues.append(contentsOf: try corpus.normalizedRow(row))
                }
                contextMasks.append(contentsOf: targetMask)
            }
            targetValues.append(contentsOf: try corpus.normalizedRow(example.targetRow))
            targetMasks.append(contentsOf: targetMask)
            targetIDs.append(targetID); known.append(1)
            descriptors.append(contentsOf: corpus.targets[example.targetIndex].descriptors.map(Float.init))
        }
        return .init(
            context: MLXArray(contextValues, [batches, contexts, features]),
            contextFeatureMask: MLXArray(contextMasks, [batches, contexts, features]),
            targetIDs: MLXArray(targetIDs, [batches]),
            knownTargetMask: MLXArray(known, [batches, 1]),
            descriptors: MLXArray(descriptors, [batches, corpus.descriptorCount]),
            targets: MLXArray(targetValues, [batches, features]),
            targetFeatureMask: MLXArray(targetMasks, [batches, features]),
            targetValues: targetValues, contextRows: [])
    }

    /// Test-only equivalence probe for the ephemeral training-control cache.
    /// It compares the exact dense MLX batch built from the source store with
    /// the batch assembled from cached, source-derived control rows.
    static func trainingControlCacheMatchesUncachedBatchForTesting(corpus: VivoCellResponseCompositeCorpusReader,
                                                                    plan: VivoCellResponseTrainingPlan,
                                                                    step: UInt64) throws -> Bool {
        try plan.validate(for: corpus)
        let strata = try stratifiedExamples(try corpus.examples(in: .training))
        let sourceStrata = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: strata)
        let bindings = try corpus.trainingTargetBindings()
        let targetIndices = try trainedTargetIndices(corpus: corpus, bindings: bindings)
        let uncached = try trainingBatch(corpus: corpus, strata: strata, sourceStrata: sourceStrata, plan: plan,
                                         trainedTargetIndices: targetIndices, step: step, controlCache: nil)
        let cache = try trainingControlCache(corpus: corpus, strata: strata)
        let cached = try trainingBatch(corpus: corpus, strata: strata, sourceStrata: sourceStrata, plan: plan,
                                       trainedTargetIndices: targetIndices, step: step, controlCache: cache)
        return uncached.context.shape == cached.context.shape &&
            uncached.context.asArray(Float.self) == cached.context.asArray(Float.self) &&
            uncached.contextFeatureMask.asArray(Float.self) == cached.contextFeatureMask.asArray(Float.self) &&
            uncached.targets.asArray(Float.self) == cached.targets.asArray(Float.self) &&
            uncached.targetFeatureMask.asArray(Float.self) == cached.targetFeatureMask.asArray(Float.self) &&
            uncached.targetIDs.asArray(Int32.self) == cached.targetIDs.asArray(Int32.self) &&
            uncached.descriptors.asArray(Float.self) == cached.descriptors.asArray(Float.self)
    }

    static func trainingControlCacheMatchesUncachedBatchForTesting(corpus: VivoCellResponseCorpusReader,
                                                                    plan: VivoCellResponseTrainingPlan,
                                                                    step: UInt64) throws -> Bool {
        try trainingControlCacheMatchesUncachedBatchForTesting(
            corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]), plan: plan, step: step)
    }

    /// Test-only in-memory comparison between the production full-axis decoder
    /// route and the former unscaled SGD route. Neither model is serialized,
    /// so no test optimizer can acquire production artifact provenance.
    static func outputAxisRouteDifferenceForTesting(corpus: VivoCellResponseCompositeCorpusReader,
                                                     plan: VivoCellResponseTrainingPlan) throws -> Float {
        try withGPUExecution {
            try plan.validate(for: corpus)
            let bindings = try corpus.trainingTargetBindings()
            let targetIndices = try trainedTargetIndices(corpus: corpus, bindings: bindings)
            let responsePrior = try trainingOnlyTargetResponsePrior(corpus: corpus,
                                                                     trainedTargetIndices: targetIndices)
            MLXRandom.seed(plan.seed)
            let scaled = VivoCellResponseMLXModel(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                                   descriptorCount: corpus.descriptorCount, architecture: plan.architecture,
                                                   targetResponsePriorValues: responsePrior)
            _ = try runSteps(model: scaled, corpus: corpus, plan: plan, startingStep: 0,
                             additionalSteps: plan.steps, trainedTargetIndices: targetIndices)
            MLXRandom.seed(plan.seed)
            let unscaled = VivoCellResponseMLXModel(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                                     descriptorCount: corpus.descriptorCount, architecture: plan.architecture,
                                                     targetResponsePriorValues: responsePrior)
            _ = try runSteps(model: unscaled, corpus: corpus, plan: plan, startingStep: 0,
                             additionalSteps: plan.steps, trainedTargetIndices: targetIndices,
                             meanHeadScaleOverrideForTesting: 1)
            let strata = try stratifiedExamples(try corpus.examples(in: .training))
            let sourceStrata = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: strata)
            let cache = try trainingControlCache(corpus: corpus, strata: strata)
            let batch = try trainingBatch(corpus: corpus, strata: strata, sourceStrata: sourceStrata, plan: plan,
                                          trainedTargetIndices: targetIndices, step: 1, controlCache: cache)
            let scaledDelta = scaled.outputs(context: batch.context, contextFeatureMask: batch.contextFeatureMask,
                                              targetIDs: batch.targetIDs,
                                              knownTargetMask: batch.knownTargetMask,
                                              descriptors: batch.descriptors).delta.asArray(Float.self)
            let unscaledDelta = unscaled.outputs(context: batch.context, contextFeatureMask: batch.contextFeatureMask,
                                                  targetIDs: batch.targetIDs,
                                                  knownTargetMask: batch.knownTargetMask,
                                                  descriptors: batch.descriptors).delta.asArray(Float.self)
            guard scaledDelta.count == unscaledDelta.count, !scaledDelta.isEmpty else {
                throw VivoCellResponseLearningError.invalid("output-axis test output shape")
            }
            let difference = zip(scaledDelta, unscaledDelta).map { abs($0.0 - $0.1) }.max() ?? 0
            guard difference.isFinite else {
                throw VivoCellResponseLearningError.invalid("output-axis test difference")
            }
            return difference
        }
    }

    static func outputAxisRouteDifferenceForTesting(corpus: VivoCellResponseCorpusReader,
                                                     plan: VivoCellResponseTrainingPlan) throws -> Float {
        try outputAxisRouteDifferenceForTesting(
            corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]), plan: plan)
    }

    private static func evaluationBatch(corpus: VivoCellResponseCompositeCorpusReader,
                                        example: VivoCellResponseTrainingExample,
                                        architecture: VivoCellResponseArchitecture,
                                        trainedTargetIndices: Set<Int>) throws -> VivoCellResponseBatch {
        guard example.contextRows.count > 0 else {
            throw VivoCellResponseLearningError.invalid("evaluation target or context rows")
        }
        let known = trainedTargetIndices.contains(example.targetIndex)
        let targetID: Int32
        if known {
            guard let value = Int32(exactly: example.targetIndex) else {
                throw VivoCellResponseLearningError.invalid("evaluation target index")
            }
            targetID = value
        } else {
            targetID = 0
        }
        let selectedRows = try deterministicContextRows(example.contextRows, count: architecture.contextCells)
        let targetMask = try corpus.featureMask(forGlobalRow: example.targetRow)
        guard targetMask.count == corpus.featureCount, targetMask.contains(where: { $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("evaluation target feature mask")
        }
        var contextValues: [Float] = []
        var contextMasks: [Float] = []
        for row in selectedRows {
            guard try corpus.featureMask(forGlobalRow: row) == targetMask else {
                throw VivoCellResponseLearningError.invalid("evaluation context source feature mask")
            }
            contextValues.append(contentsOf: try corpus.normalizedRow(row))
            contextMasks.append(contentsOf: targetMask)
        }
        let targetValues = try corpus.normalizedRow(example.targetRow)
        return .init(
            context: MLXArray(contextValues, [1, architecture.contextCells, corpus.featureCount]),
            contextFeatureMask: MLXArray(contextMasks, [1, architecture.contextCells, corpus.featureCount]),
            targetIDs: MLXArray([targetID], [1]), knownTargetMask: MLXArray([known ? Float(1) : Float(0)], [1, 1]),
            descriptors: MLXArray(corpus.targets[example.targetIndex].descriptors.map(Float.init), [1, corpus.descriptorCount]),
            targets: MLXArray(targetValues, [1, corpus.featureCount]),
            targetFeatureMask: MLXArray(targetMask, [1, corpus.featureCount]), targetValues: targetValues,
            contextRows: selectedRows)
    }

    private static func evaluate(model: VivoCellResponseMLXModel, corpus: VivoCellResponseCompositeCorpusReader,
                                 examples: [VivoCellResponseTrainingExample], architecture: VivoCellResponseArchitecture,
                                 trainedTargetIndices: Set<Int>, maximumExamples: Int, seed: UInt64,
                                 salt: UInt64,
                                 rawSourceByCorpus: [String: VivoFingerprint]? = nil) throws -> (Double, Double, Double, Int,
                                                                                                     [VivoCellResponseEvaluationExample],
                                                                                                     [VivoCellResponseSourceEvaluation],
                                                                                                     [VivoCellResponseSourceTargetEvaluation]) {
        guard (1...4_096).contains(maximumExamples) else { throw VivoCellResponseLearningError.invalid("evaluation example count") }
        if let rawSourceByCorpus {
            guard Set(rawSourceByCorpus.keys) == Set(corpus.identity.sources.map(\.corpus.hex)) else {
                throw VivoCellResponseLearningError.invalid("evaluation raw source mapping")
            }
        }
        let strata = try stratifiedExamples(examples)
        let sourceStrata = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: strata)
        let take = min(maximumExamples, strata.count)
        guard take >= sourceStrata.count else {
            throw VivoCellResponseLearningError.invalid("evaluation maximum must represent every source")
        }
        var allocations = Array(repeating: 1, count: sourceStrata.count)
        var remaining = take - sourceStrata.count
        var cursor = 0
        while remaining > 0 {
            let source = cursor % sourceStrata.count
            if allocations[source] < sourceStrata[source].count {
                allocations[source] += 1
                remaining -= 1
            }
            cursor += 1
            guard cursor <= strata.count * max(1, sourceStrata.count) else {
                throw VivoCellResponseLearningError.invalid("evaluation source allocation")
            }
        }
        var selectedStrata: [(source: Int, index: Int)] = []
        selectedStrata.reserveCapacity(take)
        for source in sourceStrata.indices {
            let local = try sampledIndicesWithoutReplacement(
                count: sourceStrata[source].count, take: allocations[source], seed: seed,
                salt: salt ^ (UInt64(source) &* 0x9e3779b97f4a7c15))
            selectedStrata.append(contentsOf: local.map { (source, sourceStrata[source][$0]) })
        }
        guard selectedStrata.count == take else {
            throw VivoCellResponseLearningError.invalid("evaluation source selection")
        }
        var selection: [VivoCellResponseEvaluationExample] = []
        selection.reserveCapacity(take)
        var totals: [Int: VivoCellResponseEvaluationTotals] = [:]
        var sourceTargetTotals: [VivoCellResponseRawSourceTargetKey: VivoCellResponseSourceTargetEvaluationTotals] = [:]
        for selected in selectedStrata {
            let candidates = strata[selected.index]
            let row = try sampledIndex(count: candidates.count, seed: seed, step: UInt64(selected.index),
                                       lane: UInt64(selected.source), salt: salt ^ 0x38b64ea1)
            let example = candidates[row]
            let source = try corpus.sourceIndex(forGlobalRow: example.targetRow)
            guard source == selected.source else {
                throw VivoCellResponseLearningError.invalid("evaluation source selection binding")
            }
            let sourceCorpus = try corpus.sourceCorpus(forGlobalRow: example.targetRow)
            var total = totals[source] ?? .init(corpus: sourceCorpus)
            guard total.corpus == sourceCorpus else {
                throw VivoCellResponseLearningError.invalid("evaluation source receipt")
            }
            total.examples += 1
            let sourceTargetKey: VivoCellResponseRawSourceTargetKey?
            if let rawSourceByCorpus {
                guard let rawSource = rawSourceByCorpus[sourceCorpus.hex] else {
                    throw VivoCellResponseLearningError.invalid("evaluation raw source")
                }
                let targetID = corpus.targets[example.targetIndex].id
                let key = VivoCellResponseRawSourceTargetKey(source: rawSource.hex, targetID: targetID)
                var sourceTargetTotal = sourceTargetTotals[key] ?? .init(source: rawSource, targetID: targetID)
                guard sourceTargetTotal.source == rawSource, sourceTargetTotal.targetID == targetID else {
                    throw VivoCellResponseLearningError.invalid("evaluation raw source target")
                }
                sourceTargetTotal.examples += 1
                sourceTargetTotals[key] = sourceTargetTotal
                sourceTargetKey = key
            } else {
                sourceTargetKey = nil
            }
            let batch = try evaluationBatch(corpus: corpus, example: example, architecture: architecture,
                                            trainedTargetIndices: trainedTargetIndices)
            selection.append(.init(targetRow: example.targetRow, contextRows: batch.contextRows))
            let output = model.outputs(context: batch.context, contextFeatureMask: batch.contextFeatureMask,
                                       targetIDs: batch.targetIDs,
                                       knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            let means = output.mean.asArray(Float.self), variances = output.variance.asArray(Float.self),
                baselines = output.baseline.asArray(Float.self)
            guard means.count == corpus.featureCount, variances.count == corpus.featureCount,
                  baselines.count == corpus.featureCount else {
                throw VivoCellResponseLearningError.invalid("evaluation output shape")
            }
            let targetMask = batch.targetFeatureMask.asArray(Float.self)
            guard targetMask.count == corpus.featureCount else {
                throw VivoCellResponseLearningError.invalid("evaluation target mask shape")
            }
            for feature in 0..<corpus.featureCount where targetMask[feature] > 0 {
                let mean = Double(means[feature]), variance = Double(variances[feature]),
                    baseline = Double(baselines[feature]), target = Double(batch.targetValues[feature])
                guard mean.isFinite, variance.isFinite, variance > 0, baseline.isFinite, target.isFinite else {
                    throw VivoCellResponseLearningError.invalid("nonfinite evaluation output")
                }
                let residual = mean - target
                total.negativeLogLikelihood += 0.5 * (residual * residual / variance + log(variance) + log(2.0 * Double.pi))
                total.squaredError += residual * residual
                let baselineResidual = baseline - target
                total.baselineSquaredError += baselineResidual * baselineResidual
                total.observedFeatures += 1
                if let sourceTargetKey {
                    guard var sourceTargetTotal = sourceTargetTotals[sourceTargetKey] else {
                        throw VivoCellResponseLearningError.invalid("evaluation raw source target")
                    }
                    sourceTargetTotal.negativeLogLikelihood +=
                        0.5 * (residual * residual / variance + log(variance) + log(2.0 * Double.pi))
                    sourceTargetTotal.squaredError += residual * residual
                    sourceTargetTotal.baselineSquaredError += baselineResidual * baselineResidual
                    sourceTargetTotal.observedFeatures += 1
                    sourceTargetTotals[sourceTargetKey] = sourceTargetTotal
                }
            }
            totals[source] = total
        }
        let sourceMetrics = try totals.keys.sorted().map { source -> VivoCellResponseSourceEvaluation in
            guard let total = totals[source], total.examples > 0, total.observedFeatures > 0 else {
                throw VivoCellResponseLearningError.invalid("evaluation source has no measured features")
            }
            return .init(corpus: total.corpus, examples: total.examples, observedFeatures: total.observedFeatures,
                         negativeLogLikelihood: total.negativeLogLikelihood / Double(total.observedFeatures),
                         rmse: sqrt(total.squaredError / Double(total.observedFeatures)),
                         matchedControlRMSE: sqrt(total.baselineSquaredError / Double(total.observedFeatures)))
        }
        guard sourceMetrics.count == corpus.sourceCount else {
            throw VivoCellResponseLearningError.invalid("evaluation does not cover every source")
        }
        let sourceTargetMetrics: [VivoCellResponseSourceTargetEvaluation]
        if rawSourceByCorpus != nil {
            sourceTargetMetrics = try sourceTargetTotals.values.sorted { left, right in
                left.source.hex == right.source.hex ? left.targetID < right.targetID : left.source.hex < right.source.hex
            }.map { total in
                guard total.examples > 0, total.observedFeatures > 0 else {
                    throw VivoCellResponseLearningError.invalid("evaluation source target has no measured features")
                }
                return .init(source: total.source, targetID: total.targetID,
                             examples: total.examples, observedFeatures: total.observedFeatures,
                             negativeLogLikelihood: total.negativeLogLikelihood / Double(total.observedFeatures),
                             rmse: sqrt(total.squaredError / Double(total.observedFeatures)),
                             matchedControlRMSE: sqrt(total.baselineSquaredError / Double(total.observedFeatures)))
            }
            guard !sourceTargetMetrics.isEmpty else {
                throw VivoCellResponseLearningError.invalid("evaluation source target metrics")
            }
        } else {
            sourceTargetMetrics = []
        }
        let observed = sourceMetrics.reduce(0) { $0 + $1.observedFeatures }
        guard observed > 0 else { throw VivoCellResponseLearningError.invalid("evaluation has no measured features") }
        let divisor = Double(sourceMetrics.count)
        return (sourceMetrics.reduce(0) { $0 + $1.negativeLogLikelihood } / divisor,
                sourceMetrics.reduce(0) { $0 + $1.rmse } / divisor,
                sourceMetrics.reduce(0) { $0 + $1.matchedControlRMSE } / divisor,
                observed, selection, sourceMetrics, sourceTargetMetrics)
    }

    private static func runSteps(model: VivoCellResponseMLXModel, corpus: VivoCellResponseCompositeCorpusReader,
                                 plan: VivoCellResponseTrainingPlan, startingStep: UInt64,
                                 additionalSteps: Int,
                                 trainedTargetIndices: Set<Int>,
                                 meanHeadScaleOverrideForTesting: Float? = nil) throws -> VivoCellResponseTrainingMetrics {
        let examples = try corpus.examples(in: .training)
        let strata = try stratifiedExamples(examples)
        let sourceStrata = try vivoCellResponseSourceStratumIndices(corpus: corpus, strata: strata)
        let validation = (try? corpus.examples(in: .validation)) ?? []
        let controlCache = try trainingControlCache(corpus: corpus, strata: strata)
        let observedFeatureCount = corpus.featureMask.reduce(into: 0) { total, value in
            if value > 0 { total += 1 }
        }
        guard observedFeatureCount > 0 else {
            throw VivoCellResponseLearningError.invalid("training corpus has no measured features")
        }
        let outputAxisScale = meanHeadScaleOverrideForTesting ?? Float(observedFeatureCount)
        let baseLearningRate = Float(plan.learningRate)
        let meanHeadLearningRate = baseLearningRate * outputAxisScale
        let meanHeadWeightDecay = Float(plan.weightDecay) / outputAxisScale
        guard outputAxisScale.isFinite, baseLearningRate.isFinite,
              meanHeadLearningRate.isFinite, meanHeadLearningRate > 0,
              meanHeadWeightDecay.isFinite, meanHeadWeightDecay >= 0 else {
            throw VivoCellResponseLearningError.limit("mean-head optimizer scale")
        }
        // The likelihood remains averaged over measured genes. Scale only the
        // independent full-axis mean decoder back to a per-feature update;
        // scaling shared or variance parameters would multiply their update by
        // the whole transcriptome width.
        let optimizer = MultiOptimizer(
            optimizers: [
                SGD(learningRate: meanHeadLearningRate, momentum: 0, weightDecay: meanHeadWeightDecay),
                SGD(learningRate: baseLearningRate, momentum: 0, weightDecay: Float(plan.weightDecay))
            ],
            filters: [{ key, _ in key == "meanHead.weight" || key == "meanHead.bias" }])
        let lossAndGrad = valueAndGrad(model: model) { model, arrays in
            [loss(model: model, context: arrays[0], contextFeatureMask: arrays[1], targetIDs: arrays[2],
                  knownTargetMask: arrays[3], descriptors: arrays[4], targets: arrays[5],
                  targetFeatureMask: arrays[6])]
        }
        var latest: VivoCellResponseTrainingMetrics?
        for offset in 0..<additionalSteps {
            try Task.checkCancellation()
            let increment = UInt64(offset + 1)
            let stepResult = startingStep.addingReportingOverflow(increment)
            guard !stepResult.overflow else { throw VivoCellResponseLearningError.limit("checkpoint step") }
            let step = stepResult.partialValue
            let batch = try trainingBatch(corpus: corpus, strata: strata, sourceStrata: sourceStrata, plan: plan,
                                          trainedTargetIndices: trainedTargetIndices, step: step,
                                          controlCache: controlCache)
            let (values, gradients) = lossAndGrad(model, [batch.context, batch.contextFeatureMask,
                                                           batch.targetIDs, batch.knownTargetMask,
                                                           batch.descriptors, batch.targets, batch.targetFeatureMask])
            optimizer.update(model: model, gradients: gradients)
            eval(model, optimizer)
            guard let scalar = values[0].asArray(Float.self).first, scalar.isFinite else {
                throw VivoCellResponseLearningError.invalid("nonfinite training loss")
            }
            var validationNLL: Double?, validationRMSE: Double?, validationMatchedControlRMSE: Double?
            if step % UInt64(plan.validationEvery) == 0 || offset + 1 == additionalSteps, !validation.isEmpty {
                let result = try evaluate(model: model, corpus: corpus, examples: validation, architecture: plan.architecture,
                                          trainedTargetIndices: trainedTargetIndices, maximumExamples: plan.validationExamples,
                                          seed: plan.seed, salt: 0x6f3a5c29)
                validationNLL = result.0; validationRMSE = result.1; validationMatchedControlRMSE = result.2
            }
            latest = .init(step: step, trainNegativeLogLikelihood: Double(scalar),
                           validationNegativeLogLikelihood: validationNLL, validationRMSE: validationRMSE,
                           validationMatchedControlRMSE: validationMatchedControlRMSE)
        }
        guard let latest else { throw VivoCellResponseLearningError.invalid("training had no steps") }
        return latest
    }

    private static func inferenceBatch(corpus: VivoCellResponseCompositeCorpusReader, state: VivoCellResponseModelState,
                                       plan: VivoCellResponsePredictionPlan) throws -> VivoCellResponseBatch {
        try plan.validate(expectedDescriptorCount: state.descriptorCount)
        guard corpus.featureAxis == state.featureAxis,
              corpus.featureCount == state.featureAxis.featureIDs.count else {
            throw VivoCellResponseLearningError.invalid("prediction query")
        }
        let contextCorpus: VivoFingerprint
        if corpus.sourceCount == 1 {
            contextCorpus = plan.contextCorpus ?? corpus.identity.sources[0].corpus
        } else {
            guard let declared = plan.contextCorpus else {
                throw VivoCellResponseLearningError.invalid("multi-source prediction requires a context corpus")
            }
            contextCorpus = declared
        }
        var rows: [Int] = []
        for sampleID in plan.contextSampleIDs.sorted() {
            let assignment = try corpus.assignment(forSampleID: sampleID, corpus: contextCorpus)
            guard assignment.role == .control else {
                throw VivoCellResponseLearningError.invalid("prediction context sample is not a declared control")
            }
            rows.append(contentsOf: try corpus.rows(forSampleID: sampleID, corpus: contextCorpus))
        }
        guard !rows.isEmpty else { throw VivoCellResponseLearningError.invalid("prediction context has no cells") }
        let targetIndex: Int32
        let knownTarget: Float
        if plan.useTrainedTargetEmbedding {
            guard let binding = state.trainedTargetBindings.first(where: { $0.id == plan.target.id }),
                  try plan.target.fingerprint() == binding.descriptorFingerprint,
                  let index = state.targetIDs.firstIndex(of: plan.target.id), let value = Int32(exactly: index) else {
                throw VivoCellResponseLearningError.incompatible("prediction target has no trained embedding")
            }
            targetIndex = value; knownTarget = 1
        } else {
            targetIndex = 0; knownTarget = 0
        }
        let selectedRows = try deterministicContextRows(rows, count: state.architecture.contextCells)
        guard let firstRow = selectedRows.first else { throw VivoCellResponseLearningError.invalid("prediction context") }
        let featureMask = try corpus.featureMask(forGlobalRow: firstRow)
        guard featureMask.count == corpus.featureCount, featureMask.contains(where: { $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("prediction feature mask")
        }
        var contexts: [Float] = []
        var contextMasks: [Float] = []
        contexts.reserveCapacity(state.architecture.contextCells * corpus.featureCount)
        for row in selectedRows {
            guard try corpus.featureMask(forGlobalRow: row) == featureMask else {
                throw VivoCellResponseLearningError.invalid("prediction context feature mask")
            }
            contexts.append(contentsOf: try corpus.normalizedRow(row))
            contextMasks.append(contentsOf: featureMask)
        }
        let descriptors = plan.target.descriptors.map(Float.init)
        return .init(
            context: MLXArray(contexts, [1, state.architecture.contextCells, corpus.featureCount]),
            contextFeatureMask: MLXArray(contextMasks, [1, state.architecture.contextCells, corpus.featureCount]),
            targetIDs: MLXArray([targetIndex], [1]), knownTargetMask: MLXArray([knownTarget], [1, 1]),
            descriptors: MLXArray(descriptors, [1, state.descriptorCount]),
            targets: MLXArray([Float](repeating: 0, count: corpus.featureCount), [1, corpus.featureCount]),
            targetFeatureMask: MLXArray(featureMask, [1, corpus.featureCount]), targetValues: [],
            contextRows: selectedRows)
    }

    private static func publishPrediction(_ prediction: VivoCellResponsePrediction,
                                          query: VivoCellResponsePredictionPlan,
                                          modelReceipt: VivoCellResponseModelReceipt,
                                          corpus: VivoCellResponseCompositeCorpusReader,
                                          implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponsePredictionReceipt {
        try VivoCellResponseArtifactIO.requireNew(destination)
        let queryBytes = try VivoCanonicalJSON.encode(query)
        let predictionBytes = try VivoCanonicalJSON.encode(prediction)
        let receipt = try VivoCellResponsePredictionReceipt(
            schemaVersion: 2, format: predictionFormat,
            model: VivoCellResponseArtifactIO.modelFingerprint(modelReceipt),
            corpus: VivoCellResponseArtifactIO.corpusFingerprint(corpus),
            query: VivoCanonicalJSON.fingerprint(queryBytes), prediction: VivoCanonicalJSON.fingerprint(predictionBytes),
            implementation: implementation)
        let temporary = try VivoCellResponseArtifactIO.staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-prediction")
        defer { try? FileManager.default.removeItem(at: temporary) }
        try VivoCellResponseArtifactIO.write(queryBytes, to: temporary, name: "query.json", maximumBytes: 4_194_304)
        try VivoCellResponseArtifactIO.write(predictionBytes, to: temporary, name: "prediction.json", maximumBytes: 268_435_456)
        try VivoCellResponseArtifactIO.write(VivoCanonicalJSON.encode(receipt), to: temporary, name: "receipt.json", maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    private static func publishEvaluation(_ evaluation: VivoCellResponseEvaluation,
                                          modelReceipt: VivoCellResponseModelReceipt,
                                          corpus: VivoCellResponseCompositeCorpusReader,
                                          implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseEvaluationReceipt {
        try validateEvaluation(evaluation)
        try VivoCellResponseArtifactIO.requireNew(destination)
        let evaluationBytes = try VivoCanonicalJSON.encode(evaluation)
        let receipt = VivoCellResponseEvaluationReceipt(
            schemaVersion: 4, format: evaluationFormat,
            model: try VivoCellResponseArtifactIO.modelFingerprint(modelReceipt),
            corpus: try VivoCellResponseArtifactIO.corpusFingerprint(corpus),
            evaluation: try VivoCanonicalJSON.fingerprint(evaluationBytes), implementation: implementation)
        let temporary = try VivoCellResponseArtifactIO.staging(destination.deletingLastPathComponent(), prefix: "numivivo-cell-response-evaluation")
        defer { try? FileManager.default.removeItem(at: temporary) }
        try VivoCellResponseArtifactIO.write(evaluationBytes, to: temporary, name: "evaluation.json", maximumBytes: 134_217_728)
        try VivoCellResponseArtifactIO.write(try VivoCanonicalJSON.encode(receipt), to: temporary, name: "receipt.json", maximumBytes: 1_048_576)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    private static func validateEvaluation(_ evaluation: VivoCellResponseEvaluation) throws {
        let current = evaluation.schemaVersion == 5 && evaluation.format == evaluationFormat
        let previousSourceMetrics = evaluation.schemaVersion == 4 && evaluation.format == cohortV4EvaluationFormat
        let previousCohort = evaluation.schemaVersion == 3 && evaluation.format == cohortV3EvaluationFormat
        let legacy = evaluation.schemaVersion == 2 && evaluation.format == legacyEvaluationFormat
        let hasSourceMetrics = current || previousSourceMetrics
        guard hasSourceMetrics || previousCohort || legacy else {
            throw VivoCellResponseLearningError.invalid("evaluation artifact")
        }
        guard
              (1...4_096).contains(evaluation.maximumExamples),
              evaluation.samplerVersion == samplerVersion,
              evaluation.examples == evaluation.selection.count, !evaluation.selection.isEmpty,
              evaluation.observedFeatures > 0, evaluation.negativeLogLikelihood.isFinite,
              evaluation.rmse.isFinite, evaluation.rmse >= 0,
              evaluation.matchedControlRMSE.isFinite, evaluation.matchedControlRMSE >= 0,
              Set(evaluation.selection.map(\.targetRow)).count == evaluation.selection.count,
              evaluation.selection.allSatisfy({ $0.targetRow >= 0 && (1...512).contains($0.contextRows.count) && $0.contextRows.allSatisfy({ $0 >= 0 }) }) else {
            throw VivoCellResponseLearningError.invalid("evaluation artifact")
        }
        guard hasSourceMetrics || evaluation.sourceMetrics.isEmpty else {
            throw VivoCellResponseLearningError.invalid("historical evaluation source metrics")
        }
        guard !hasSourceMetrics ||
                ((1...64).contains(evaluation.sourceMetrics.count) &&
                 evaluation.sourceMetrics.map(\.corpus.hex) == evaluation.sourceMetrics.map(\.corpus.hex).sorted() &&
                 Set(evaluation.sourceMetrics.map(\.corpus.hex)).count == evaluation.sourceMetrics.count &&
                 evaluation.sourceMetrics.allSatisfy { metric in
                     metric.examples > 0 && metric.observedFeatures > 0 && metric.negativeLogLikelihood.isFinite &&
                     metric.rmse.isFinite && metric.rmse >= 0 && metric.matchedControlRMSE.isFinite &&
                     metric.matchedControlRMSE >= 0
                 } &&
                 evaluation.sourceMetrics.reduce(0, { $0 + $1.examples }) == evaluation.examples &&
                 evaluation.sourceMetrics.reduce(0, { $0 + $1.observedFeatures }) == evaluation.observedFeatures) else {
            throw VivoCellResponseLearningError.invalid("evaluation source metrics")
        }
        guard current || evaluation.sourceTargetMetrics.isEmpty else {
            throw VivoCellResponseLearningError.invalid("historical evaluation source target metrics")
        }
        guard !current ||
                ((1...evaluation.examples).contains(evaluation.sourceTargetMetrics.count) &&
                 evaluation.sourceTargetMetrics.allSatisfy { metric in
                     vivoOmicsID(metric.targetID) && metric.examples > 0 && metric.observedFeatures > 0 &&
                         metric.negativeLogLikelihood.isFinite && metric.rmse.isFinite && metric.rmse >= 0 &&
                         metric.matchedControlRMSE.isFinite && metric.matchedControlRMSE >= 0
                 } &&
                 evaluation.sourceTargetMetrics.indices.dropFirst().allSatisfy { index in
                     let previous = evaluation.sourceTargetMetrics[index - 1]
                     let current = evaluation.sourceTargetMetrics[index]
                     return previous.source.hex < current.source.hex ||
                         (previous.source == current.source && previous.targetID < current.targetID)
                 } &&
                 Set(evaluation.sourceTargetMetrics.map {
                     VivoCellResponseRawSourceTargetKey(source: $0.source.hex, targetID: $0.targetID)
                 }).count ==
                    evaluation.sourceTargetMetrics.count &&
                 evaluation.sourceTargetMetrics.reduce(0, { $0 + $1.examples }) == evaluation.examples &&
                 evaluation.sourceTargetMetrics.reduce(0, { $0 + $1.observedFeatures }) == evaluation.observedFeatures) else {
            throw VivoCellResponseLearningError.invalid("evaluation source target metrics")
        }
        if hasSourceMetrics {
            let divisor = Double(evaluation.sourceMetrics.count)
            let sourceNLL = evaluation.sourceMetrics.reduce(0) { $0 + $1.negativeLogLikelihood } / divisor
            let sourceRMSE = evaluation.sourceMetrics.reduce(0) { $0 + $1.rmse } / divisor
            let sourceBaseline = evaluation.sourceMetrics.reduce(0) { $0 + $1.matchedControlRMSE } / divisor
            func agrees(_ actual: Double, _ expected: Double) -> Bool {
                abs(actual - expected) <= 1e-10 * max(1, abs(expected))
            }
            guard agrees(evaluation.negativeLogLikelihood, sourceNLL), agrees(evaluation.rmse, sourceRMSE),
                  agrees(evaluation.matchedControlRMSE, sourceBaseline) else {
                throw VivoCellResponseLearningError.invalid("evaluation aggregate source metrics")
            }
        }
    }

    /// Require a verified response evaluation to beat its exact matched-control
    /// baseline. Raw evaluations remain publishable for diagnosis; callers use
    /// this explicit gate before treating one as a qualified engineering result.
    public static func requireMatchedControlImprovement(_ evaluation: VivoCellResponseEvaluation) throws {
        try validateEvaluation(evaluation)
        guard evaluation.rmse < evaluation.matchedControlRMSE else {
            throw VivoCellResponseLearningError.invalid(
                "evaluation does not improve the exact matched-control baseline")
        }
        guard (evaluation.schemaVersion != 4 && evaluation.schemaVersion != 5) || evaluation.sourceMetrics.allSatisfy({
            $0.rmse < $0.matchedControlRMSE
        }) else {
            throw VivoCellResponseLearningError.invalid(
                "evaluation does not improve the exact matched-control baseline for every source")
        }
        guard evaluation.schemaVersion != 5 || evaluation.sourceTargetMetrics.allSatisfy({
            $0.rmse < $0.matchedControlRMSE
        }) else {
            throw VivoCellResponseLearningError.invalid(
                "evaluation does not improve the exact matched-control baseline for every raw source target")
        }
    }

    private static func validatePrediction(_ prediction: VivoCellResponsePrediction) throws {
        try prediction.featureAxis.validate()
        guard (prediction.schemaVersion == 2 && prediction.format == predictionFormat) ||
              (prediction.schemaVersion == 1 && prediction.format == legacyPredictionFormat) else {
            throw VivoCellResponseLearningError.invalid("prediction output")
        }
        guard vivoOmicsID(prediction.id),
              prediction.featureAxis.featureIDs.count == prediction.featureMask.count,
              prediction.featureMask.count == prediction.meanLogCPM.count,
              prediction.featureMask.count == prediction.contextBaselineLogCPM.count,
              prediction.featureMask.count == prediction.meanDeltaLogCPM.count,
              prediction.meanLogCPM.count == prediction.varianceLogCPM.count,
              (1...512).contains(prediction.contextSourceRows.count),
              prediction.contextSourceRows.allSatisfy({ $0 >= 0 }),
              prediction.featureMask.allSatisfy({ $0 == 0 || $0 == 1 }),
              prediction.featureMask.contains(where: { $0 > 0 }),
              prediction.contextBaselineLogCPM.allSatisfy(\.isFinite),
              prediction.meanDeltaLogCPM.allSatisfy(\.isFinite),
              prediction.meanLogCPM.allSatisfy(\.isFinite),
              prediction.varianceLogCPM.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoCellResponseLearningError.invalid("prediction output")
        }
        guard zip(zip(prediction.meanLogCPM, prediction.contextBaselineLogCPM), prediction.meanDeltaLogCPM).allSatisfy({ values in
            let expected = values.0.1 + values.1
            let scale = max(Float(1), max(abs(values.0.0), abs(expected)))
            return abs(values.0.0 - expected) <= 1e-5 * scale
        }) else {
            throw VivoCellResponseLearningError.invalid("prediction baseline and residual disagree")
        }
    }

    private static func approximatelyEqual(_ expected: [Float], _ actual: [Float]) -> Bool {
        guard expected.count == actual.count else { return false }
        return zip(expected, actual).allSatisfy { values in
            let scale = max(Float(1), max(abs(values.0), abs(values.1)))
            return abs(values.0 - values.1) <= 1e-5 * scale
        }
    }

    private static func approximatelyEqual(_ expected: [VivoCellResponseSourceEvaluation],
                                           _ actual: [VivoCellResponseSourceEvaluation]) -> Bool {
        guard expected.count == actual.count else { return false }
        return zip(expected, actual).allSatisfy { values in
            let left = values.0, right = values.1
            guard left.corpus == right.corpus, left.examples == right.examples,
                  left.observedFeatures == right.observedFeatures else {
                return false
            }
            func matches(_ expected: Double, _ actual: Double) -> Bool {
                abs(expected - actual) <= 1e-8 * max(1, abs(expected))
            }
            return matches(left.negativeLogLikelihood, right.negativeLogLikelihood) &&
                matches(left.rmse, right.rmse) && matches(left.matchedControlRMSE, right.matchedControlRMSE)
        }
    }

    private static func approximatelyEqual(_ expected: [VivoCellResponseSourceTargetEvaluation],
                                           _ actual: [VivoCellResponseSourceTargetEvaluation]) -> Bool {
        guard expected.count == actual.count else { return false }
        return zip(expected, actual).allSatisfy { values in
            let left = values.0, right = values.1
            guard left.source == right.source, left.targetID == right.targetID,
                  left.examples == right.examples, left.observedFeatures == right.observedFeatures else {
                return false
            }
            func matches(_ expected: Double, _ actual: Double) -> Bool {
                abs(expected - actual) <= 1e-8 * max(1, abs(expected))
            }
            return matches(left.negativeLogLikelihood, right.negativeLogLikelihood) &&
                matches(left.rmse, right.rmse) && matches(left.matchedControlRMSE, right.matchedControlRMSE)
        }
    }

    /// Fit a model only after replaying a source-bound biological cohort
    /// admission against the exact verified reader. Smaller development
    /// cohorts remain source-bound but cannot later qualify an evaluation.
    public static func train(corpus: VivoCellResponseCompositeCorpusReader, plan: VivoCellResponseTrainingPlan,
                             cohort: VivoCellResponseCohortAdmission,
                             implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseModelReceipt {
        try plan.validate(for: corpus)
        try VivoCellResponseCohort.verify(cohort, readers: corpus.sourceReaders)
        let cohortBytes = try VivoCanonicalJSON.encode(cohort)
        let cohortFingerprint = try VivoCanonicalJSON.fingerprint(cohortBytes)
        try VivoCellResponseArtifactIO.requireNew(destination)
        return try withGPUExecution {
            let planBytes = try VivoCanonicalJSON.encode(plan)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            let trainedTargetBindings = try corpus.trainingTargetBindings()
            let trainedTargetIndices = try trainedTargetIndices(corpus: corpus, bindings: trainedTargetBindings)
            MLXRandom.seed(plan.seed)
            let responsePrior = try trainingOnlyTargetResponsePrior(corpus: corpus,
                                                                     trainedTargetIndices: trainedTargetIndices)
            let model = VivoCellResponseMLXModel(featureCount: corpus.featureCount, targetCount: corpus.targetCount,
                                                 descriptorCount: corpus.descriptorCount, architecture: plan.architecture,
                                                 targetResponsePriorValues: responsePrior)
            let metrics = try runSteps(model: model, corpus: corpus, plan: plan, startingStep: 0,
                                       additionalSteps: plan.steps, trainedTargetIndices: trainedTargetIndices)
            let state = VivoCellResponseModelState(
                schemaVersion: 6, architecture: plan.architecture, featureAxis: corpus.featureAxis,
                targetIDs: corpus.targets.map(\.id), trainedTargetBindings: trainedTargetBindings,
                descriptorCount: corpus.descriptorCount,
                corpus: corpusFingerprint, trainingPlan: try VivoCanonicalJSON.fingerprint(planBytes),
                cohort: cohortFingerprint,
                step: metrics.step, seed: plan.seed, samplerVersion: samplerVersion,
                targetResponsePriorVersion: targetResponsePriorVersion,
                optimizerVersion: optimizerVersion,
                learningRate: plan.learningRate, weightDecay: plan.weightDecay, latestMetrics: metrics)
            return try VivoCellResponseArtifactIO.publish(model: model, state: state, trainingPlanBytes: planBytes,
                                                           implementation: implementation, parentCheckpoint: nil,
                                                           cohortAdmission: cohort, cohortAdmissionBytes: cohortBytes,
                                                           to: destination)
        }
    }

    /// Continue a checkpoint without changing model shape, optimizer semantics,
    /// source corpus identity, or deterministic batch sampler.
    public static func resume(model directory: URL, corpus: VivoCellResponseCompositeCorpusReader,
                              plan: VivoCellResponseResumePlan, implementation: VivoFingerprint,
                              to destination: URL) throws -> VivoCellResponseModelReceipt {
        return try withGPUExecution {
            try plan.validate()
            try VivoCellResponseArtifactIO.requireNew(destination)
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            guard let cohort = loaded.cohortAdmission, let cohortBytes = loaded.cohortAdmissionBytes,
                  loaded.state.cohort == loaded.receipt.cohort else {
                throw VivoCellResponseLearningError.incompatible("resume requires a cohort-bound model")
            }
            guard loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
                throw VivoCellResponseLearningError.incompatible("resume requires a v7 cohort-bound model")
            }
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("resume corpus differs from checkpoint")
            }
            try VivoCellResponseCohort.verify(cohort, readers: corpus.sourceReaders)
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            try loaded.trainingPlan.validate(for: corpus)
            let metrics = try runSteps(model: loaded.model, corpus: corpus, plan: loaded.trainingPlan,
                                       startingStep: loaded.state.step, additionalSteps: plan.additionalSteps,
                                       trainedTargetIndices: trainedTargetIndices)
            let state = VivoCellResponseModelState(
                schemaVersion: 6, architecture: loaded.state.architecture, featureAxis: loaded.state.featureAxis,
                targetIDs: loaded.state.targetIDs, trainedTargetBindings: loaded.state.trainedTargetBindings,
                descriptorCount: loaded.state.descriptorCount,
                corpus: loaded.state.corpus, trainingPlan: loaded.state.trainingPlan, cohort: loaded.state.cohort,
                step: metrics.step, seed: loaded.state.seed, samplerVersion: loaded.state.samplerVersion,
                targetResponsePriorVersion: loaded.state.targetResponsePriorVersion,
                optimizerVersion: loaded.state.optimizerVersion,
                learningRate: loaded.state.learningRate, weightDecay: loaded.state.weightDecay, latestMetrics: metrics)
            return try VivoCellResponseArtifactIO.publish(model: loaded.model, state: state,
                                                           trainingPlanBytes: loaded.trainingPlanBytes,
                                                           implementation: implementation,
                                                           parentCheckpoint: loaded.receipt.checkpoint,
                                                           cohortAdmission: cohort, cohortAdmissionBytes: cohortBytes,
                                                           to: destination)
        }
    }

    /// Verify receipt, checkpoint framing, safetensors metadata, and every
    /// model parameter shape without opening a source count store.
    public static func verifyModel(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseModelReceipt {
        try VivoCellResponseArtifactIO.load(directory, implementation: implementation).receipt
    }

    /// Score measured perturbed rows from a held partition. This is a software
    /// evaluation result; it is not a biological validation claim.
    public static func evaluate(model directory: URL, corpus: VivoCellResponseCompositeCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        return try withGPUExecution {
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            guard let cohort = loaded.cohortAdmission, loaded.state.cohort == loaded.receipt.cohort,
                  loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
                throw VivoCellResponseLearningError.incompatible("evaluation requires a v7 cohort-bound model")
            }
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("evaluation corpus differs from checkpoint")
            }
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let rawSourceByCorpus = try vivoCellResponseRawSourceByCorpus(cohort: cohort, corpus: corpus)
            let examples = try corpus.examples(in: partition)
            let result = try evaluate(model: loaded.model, corpus: corpus, examples: examples,
                                      architecture: loaded.state.architecture, trainedTargetIndices: trainedTargetIndices,
                                      maximumExamples: maximumExamples, seed: loaded.state.seed,
                                      salt: 0x4d8b9173, rawSourceByCorpus: rawSourceByCorpus)
            let evaluation = VivoCellResponseEvaluation(
                schemaVersion: 5, format: evaluationFormat, partition: partition,
                maximumExamples: maximumExamples, samplerVersion: samplerVersion, seed: loaded.state.seed,
                examples: result.4.count, observedFeatures: result.3,
                negativeLogLikelihood: result.0, rmse: result.1, matchedControlRMSE: result.2,
                selection: result.4, sourceMetrics: result.5, sourceTargetMetrics: result.6)
            try validateEvaluation(evaluation)
            return evaluation
        }
    }

    /// Publish an immutable evaluation receipt after scoring a verified model
    /// and corpus. The stored selection makes later review and replay exact.
    public static func evaluate(model directory: URL, corpus: VivoCellResponseCompositeCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseEvaluationReceipt {
        let evaluation = try evaluate(model: directory, corpus: corpus, partition: partition,
                                      maximumExamples: maximumExamples, implementation: implementation)
        let modelReceipt = try VivoCellResponseArtifactIO.load(directory, implementation: implementation).receipt
        return try publishEvaluation(evaluation, modelReceipt: modelReceipt, corpus: corpus,
                                     implementation: implementation, to: destination)
    }

    /// Predict a distribution from declared control samples. The output is a
    /// compact response vector artifact, never an H5AD count matrix.
    public static func predict(model directory: URL, corpus: VivoCellResponseCompositeCorpusReader,
                               plan: VivoCellResponsePredictionPlan, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoCellResponsePredictionReceipt {
        return try withGPUExecution {
            let loaded = try VivoCellResponseArtifactIO.load(directory, implementation: implementation)
            guard loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
                throw VivoCellResponseLearningError.incompatible("prediction requires a v7 cohort-bound model")
            }
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard loaded.state.corpus == corpusFingerprint else {
                throw VivoCellResponseLearningError.incompatible("prediction corpus differs from checkpoint")
            }
            _ = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let batch = try inferenceBatch(corpus: corpus, state: loaded.state, plan: plan)
            let output = loaded.model.outputs(context: batch.context, contextFeatureMask: batch.contextFeatureMask,
                                              targetIDs: batch.targetIDs,
                                              knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            let means = output.mean.asArray(Float.self), variances = output.variance.asArray(Float.self),
                baselines = output.baseline.asArray(Float.self), deltas = output.delta.asArray(Float.self)
            guard means.count == corpus.featureCount, variances.count == corpus.featureCount,
                  baselines.count == corpus.featureCount, deltas.count == corpus.featureCount,
                  means.allSatisfy(\.isFinite), baselines.allSatisfy(\.isFinite), deltas.allSatisfy(\.isFinite),
                  variances.allSatisfy({ $0.isFinite && $0 > 0 }) else {
                throw VivoCellResponseLearningError.invalid("prediction output")
            }
            let featureMask = batch.targetFeatureMask.asArray(Float.self)
            let prediction = VivoCellResponsePrediction(schemaVersion: 2, format: predictionFormat, id: plan.id,
                                                        featureAxis: corpus.featureAxis, featureMask: featureMask,
                                                        contextBaselineLogCPM: baselines, meanDeltaLogCPM: deltas,
                                                        contextSourceRows: batch.contextRows,
                                                        meanLogCPM: means, varianceLogCPM: variances)
            return try publishPrediction(prediction, query: plan, modelReceipt: loaded.receipt, corpus: corpus,
                                         implementation: implementation, to: destination)
        }
    }

    private static func loadVerifiedPrediction(_ directory: URL,
                                               implementation: VivoFingerprint) throws -> (VivoCellResponsePredictionReceipt, VivoCellResponsePredictionPlan, VivoCellResponsePrediction) {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("prediction artifact directory")
        }
        let receiptBytes = try VivoCellResponseArtifactIO.read(root, "receipt.json", maximumBytes: 1_048_576)
        let queryBytes = try VivoCellResponseArtifactIO.read(root, "query.json", maximumBytes: 4_194_304)
        let predictionBytes = try VivoCellResponseArtifactIO.read(root, "prediction.json", maximumBytes: 268_435_456)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponsePredictionReceipt.self, from: receiptBytes)
        let query = try VivoCanonicalJSON.decode(VivoCellResponsePredictionPlan.self, from: queryBytes)
        let prediction = try VivoCanonicalJSON.decode(VivoCellResponsePrediction.self, from: predictionBytes)
        guard ((receipt.schemaVersion == 2 && receipt.format == predictionFormat) ||
               (receipt.schemaVersion == 1 && receipt.format == legacyPredictionFormat)),
              receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes,
              try VivoCanonicalJSON.fingerprint(queryBytes) == receipt.query,
              try VivoCanonicalJSON.fingerprint(predictionBytes) == receipt.prediction,
              try VivoCanonicalJSON.encode(query) == queryBytes,
              try VivoCanonicalJSON.encode(prediction) == predictionBytes,
              prediction.id == query.id else {
            throw VivoCellResponseLearningError.invalid("prediction receipt or output")
        }
        try query.validate()
        try validatePrediction(prediction)
        return (receipt, query, prediction)
    }

    /// Verify a prediction artifact's self-consistency, including its own
    /// query and output fingerprints. This does not establish a model/corpus
    /// provenance chain; use the bound overload when those artifacts exist.
    public static func verifyPrediction(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponsePrediction {
        try loadVerifiedPrediction(directory, implementation: implementation).2
    }

    /// Verify a prediction against the exact model receipt and source-bound
    /// corpus used to create it. This is the production verification route.
    public static func verifyPrediction(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCompositeCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponsePrediction {
        return try withGPUExecution {
            let artifact = try loadVerifiedPrediction(directory, implementation: implementation)
            let loaded = try VivoCellResponseArtifactIO.load(modelDirectory, implementation: implementation)
            guard loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
                throw VivoCellResponseLearningError.incompatible("bound prediction verification requires a v7 cohort-bound model")
            }
            let modelFingerprint = try VivoCellResponseArtifactIO.modelFingerprint(loaded.receipt)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard artifact.0.model == modelFingerprint,
                  artifact.0.corpus == corpusFingerprint,
                  artifact.2.featureAxis == corpus.featureAxis else {
                throw VivoCellResponseLearningError.incompatible("prediction model or corpus provenance differs")
            }
            // Receipt hashes prove stored bytes. Rebuild the exact input and
            // output under the supplied corpus/model as well, so neither an
            // impossible query nor substituted response values can pass.
            _ = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let batch = try inferenceBatch(corpus: corpus, state: loaded.state, plan: artifact.1)
            let output = loaded.model.outputs(context: batch.context, contextFeatureMask: batch.contextFeatureMask,
                                              targetIDs: batch.targetIDs,
                                              knownTargetMask: batch.knownTargetMask, descriptors: batch.descriptors)
            guard artifact.2.contextSourceRows == batch.contextRows,
                  artifact.2.featureMask == batch.targetFeatureMask.asArray(Float.self),
                  approximatelyEqual(artifact.2.contextBaselineLogCPM, output.baseline.asArray(Float.self)),
                  approximatelyEqual(artifact.2.meanDeltaLogCPM, output.delta.asArray(Float.self)),
                  approximatelyEqual(artifact.2.meanLogCPM, output.mean.asArray(Float.self)),
                  approximatelyEqual(artifact.2.varianceLogCPM, output.variance.asArray(Float.self)) else {
                throw VivoCellResponseLearningError.incompatible("prediction control context or model output differs")
            }
            return artifact.2
        }
    }

    private static func loadVerifiedEvaluation(_ directory: URL,
                                               implementation: VivoFingerprint) throws -> (VivoCellResponseEvaluationReceipt, VivoCellResponseEvaluation) {
        let root: VivoRootedDocumentStore
        do {
            root = try .init(rootURL: directory)
        } catch {
            throw VivoCellResponseLearningError.invalid("evaluation artifact directory")
        }
        let receiptBytes = try VivoCellResponseArtifactIO.read(root, "receipt.json", maximumBytes: 1_048_576)
        let evaluationBytes = try VivoCellResponseArtifactIO.read(root, "evaluation.json", maximumBytes: 134_217_728)
        let receipt = try VivoCanonicalJSON.decode(VivoCellResponseEvaluationReceipt.self, from: receiptBytes)
        let evaluation = try VivoCanonicalJSON.decode(VivoCellResponseEvaluation.self, from: evaluationBytes)
        let pairedFormat =
            (receipt.schemaVersion == 4 && receipt.format == evaluationFormat &&
             evaluation.schemaVersion == 5 && evaluation.format == evaluationFormat) ||
            (receipt.schemaVersion == 3 && receipt.format == cohortV4EvaluationFormat &&
             evaluation.schemaVersion == 4 && evaluation.format == cohortV4EvaluationFormat) ||
            (receipt.schemaVersion == 2 && receipt.format == cohortV3EvaluationFormat &&
             evaluation.schemaVersion == 3 && evaluation.format == cohortV3EvaluationFormat) ||
            (receipt.schemaVersion == 1 && receipt.format == legacyEvaluationFormat &&
             evaluation.schemaVersion == 2 && evaluation.format == legacyEvaluationFormat)
        guard pairedFormat,
              receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == receiptBytes,
              try VivoCanonicalJSON.encode(evaluation) == evaluationBytes,
              try VivoCanonicalJSON.fingerprint(evaluationBytes) == receipt.evaluation else {
            throw VivoCellResponseLearningError.invalid("evaluation receipt")
        }
        try validateEvaluation(evaluation)
        return (receipt, evaluation)
    }

    /// Verify an evaluation artifact's internal canonical documents. Bind it to
    /// a model and corpus with the overload below before interpreting metrics.
    public static func verifyEvaluation(_ directory: URL, implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        try loadVerifiedEvaluation(directory, implementation: implementation).1
    }

    /// Recompute a published evaluation under the exact model and corpus. This
    /// verifies both its frozen target/control selection and reported metrics.
    public static func verifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCompositeCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        return try withGPUExecution {
            let artifact = try loadVerifiedEvaluation(directory, implementation: implementation)
            let loaded = try VivoCellResponseArtifactIO.load(modelDirectory, implementation: implementation)
            guard let cohort = loaded.cohortAdmission, loaded.state.cohort == loaded.receipt.cohort,
                  loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
                throw VivoCellResponseLearningError.incompatible("bound evaluation verification requires a v7 cohort-bound model")
            }
            guard artifact.0.schemaVersion == 4, artifact.0.format == evaluationFormat,
                  artifact.1.schemaVersion == 5, artifact.1.format == evaluationFormat else {
                throw VivoCellResponseLearningError.incompatible("bound evaluation verification requires a v5 evaluation")
            }
            let modelFingerprint = try VivoCellResponseArtifactIO.modelFingerprint(loaded.receipt)
            let corpusFingerprint = try VivoCellResponseArtifactIO.corpusFingerprint(corpus)
            guard artifact.0.model == modelFingerprint, artifact.0.corpus == corpusFingerprint,
                  artifact.1.seed == loaded.state.seed else {
                throw VivoCellResponseLearningError.incompatible("evaluation model or corpus provenance differs")
            }
            let trainedTargetIndices = try validateCheckpointCorpus(loaded.state, corpus: corpus)
            let rawSourceByCorpus = try vivoCellResponseRawSourceByCorpus(cohort: cohort, corpus: corpus)
            let examples = try corpus.examples(in: artifact.1.partition)
            let result = try evaluate(model: loaded.model, corpus: corpus, examples: examples,
                                      architecture: loaded.state.architecture, trainedTargetIndices: trainedTargetIndices,
                                      maximumExamples: artifact.1.maximumExamples, seed: artifact.1.seed,
                                      salt: 0x4d8b9173, rawSourceByCorpus: rawSourceByCorpus)
            guard artifact.1.examples == result.4.count, artifact.1.observedFeatures == result.3,
                  artifact.1.selection == result.4,
                  approximatelyEqual(artifact.1.sourceMetrics, result.5),
                  approximatelyEqual(artifact.1.sourceTargetMetrics, result.6),
                  abs(artifact.1.negativeLogLikelihood - result.0) <= 1e-8 * max(1, abs(result.0)),
                  abs(artifact.1.rmse - result.1) <= 1e-8 * max(1, abs(result.1)),
                  abs(artifact.1.matchedControlRMSE - result.2) <= 1e-8 * max(1, abs(result.2)) else {
                throw VivoCellResponseLearningError.incompatible("evaluation selection or metrics differ")
            }
            return artifact.1
        }
    }

    /// Qualification is a held-out result, not an arbitrary evaluation
    /// sample. Every declared test target/sample stratum must appear exactly
    /// once in the frozen evaluation selection. `verifyEvaluation` supplies
    /// the companion guarantee that each selected target and control row was
    /// recomputed from the bound model and corpus.
    static func requireQualificationEvaluationCoverage(_ evaluation: VivoCellResponseEvaluation,
                                                       corpus: VivoCellResponseCompositeCorpusReader,
                                                       cohort: VivoCellResponseCohortAdmission) throws {
        try validateEvaluation(evaluation)
        guard evaluation.schemaVersion == 5, evaluation.format == evaluationFormat,
              evaluation.sourceMetrics.map(\.corpus) == corpus.identity.sources.map(\.corpus),
              evaluation.partition == .test else {
            throw VivoCellResponseLearningError.invalid("qualification requires a held-out test evaluation")
        }
        let rawSourceByCorpus = try vivoCellResponseRawSourceByCorpus(cohort: cohort, corpus: corpus)
        let examples = try corpus.examples(in: .test)
        let expectedStrata = Set(examples.map(\.stratumIndex))
        guard !expectedStrata.isEmpty else {
            throw VivoCellResponseLearningError.invalid("qualification has no held-out test strata")
        }
        var stratumByTargetRow: [Int: Int] = [:]
        for example in examples {
            guard stratumByTargetRow.updateValue(example.stratumIndex, forKey: example.targetRow) == nil else {
                throw VivoCellResponseLearningError.invalid("qualification test target rows overlap")
            }
        }
        let selectedStrata = try evaluation.selection.map { example -> Int in
            guard let stratum = stratumByTargetRow[example.targetRow] else {
                throw VivoCellResponseLearningError.invalid("qualification selection is outside the held-out test")
            }
            return stratum
        }
        guard evaluation.maximumExamples >= expectedStrata.count,
              evaluation.examples == expectedStrata.count,
              selectedStrata.count == expectedStrata.count,
              Set(selectedStrata) == expectedStrata else {
            throw VivoCellResponseLearningError.invalid(
                "qualification requires exhaustive held-out test-stratum coverage")
        }
        var expectedSourceTargetKeys: Set<VivoCellResponseRawSourceTargetKey> = []
        for example in examples {
            guard let rawSource = rawSourceByCorpus[try corpus.sourceCorpus(forGlobalRow: example.targetRow).hex],
                  corpus.targets.indices.contains(example.targetIndex) else {
                throw VivoCellResponseLearningError.invalid("qualification source target binding")
            }
            expectedSourceTargetKeys.insert(.init(source: rawSource.hex,
                                                  targetID: corpus.targets[example.targetIndex].id))
        }
        let reportedSourceTargetKeys = Set(evaluation.sourceTargetMetrics.map {
            VivoCellResponseRawSourceTargetKey(source: $0.source.hex, targetID: $0.targetID)
        })
        guard reportedSourceTargetKeys == expectedSourceTargetKeys else {
            throw VivoCellResponseLearningError.invalid(
                "qualification requires exhaustive held-out raw source target coverage")
        }
    }

    /// Replay an evaluation under its bound model and corpus, then require a
    /// strict improvement over the recorded exact matched-control baseline.
    /// A failed qualification never changes the raw evaluation artifact.
    public static func qualifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                         corpus: VivoCellResponseCompositeCorpusReader,
                                         implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        let loaded = try VivoCellResponseArtifactIO.load(modelDirectory, implementation: implementation)
        guard let cohort = loaded.cohortAdmission, loaded.state.cohort == loaded.receipt.cohort else {
            throw VivoCellResponseLearningError.incompatible("qualification requires a cohort-bound model")
        }
        guard loaded.receipt.format == modelFormat, loaded.state.schemaVersion == 6 else {
            throw VivoCellResponseLearningError.incompatible("qualification requires a v7 cohort-bound model")
        }
        try VivoCellResponseCohort.verify(cohort, readers: corpus.sourceReaders)
        try VivoCellResponseCohort.requireQualificationCoverage(cohort)
        let evaluation = try verifyEvaluation(directory, model: modelDirectory, corpus: corpus,
                                              implementation: implementation)
        try requireQualificationEvaluationCoverage(evaluation, corpus: corpus, cohort: cohort)
        try requireMatchedControlImprovement(evaluation)
        return evaluation
    }

    // Keep the established one-source Swift surface while routing all current
    // execution through the same canonical composite namespace. This avoids a
    // compatibility fork between single-study and cross-study prediction.
    public static func train(corpus: VivoCellResponseCorpusReader, plan: VivoCellResponseTrainingPlan,
                             cohort: VivoCellResponseCohortAdmission,
                             implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseModelReceipt {
        try train(corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]), plan: plan, cohort: cohort,
                  implementation: implementation, to: destination)
    }

    public static func resume(model directory: URL, corpus: VivoCellResponseCorpusReader,
                              plan: VivoCellResponseResumePlan, implementation: VivoFingerprint,
                              to destination: URL) throws -> VivoCellResponseModelReceipt {
        try resume(model: directory, corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]), plan: plan,
                   implementation: implementation, to: destination)
    }

    public static func evaluate(model directory: URL, corpus: VivoCellResponseCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        try evaluate(model: directory, corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                     partition: partition, maximumExamples: maximumExamples, implementation: implementation)
    }

    public static func evaluate(model directory: URL, corpus: VivoCellResponseCorpusReader,
                                partition: VivoCellResponsePartition, maximumExamples: Int,
                                implementation: VivoFingerprint, to destination: URL) throws -> VivoCellResponseEvaluationReceipt {
        try evaluate(model: directory, corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                     partition: partition, maximumExamples: maximumExamples, implementation: implementation,
                     to: destination)
    }

    public static func predict(model directory: URL, corpus: VivoCellResponseCorpusReader,
                               plan: VivoCellResponsePredictionPlan, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoCellResponsePredictionReceipt {
        try predict(model: directory, corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]), plan: plan,
                    implementation: implementation, to: destination)
    }

    public static func verifyPrediction(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponsePrediction {
        try verifyPrediction(directory, model: modelDirectory,
                             corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                             implementation: implementation)
    }

    public static func verifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                        corpus: VivoCellResponseCorpusReader,
                                        implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        try verifyEvaluation(directory, model: modelDirectory,
                             corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                             implementation: implementation)
    }

    static func requireQualificationEvaluationCoverage(_ evaluation: VivoCellResponseEvaluation,
                                                       corpus: VivoCellResponseCorpusReader,
                                                       cohort: VivoCellResponseCohortAdmission) throws {
        try requireQualificationEvaluationCoverage(evaluation,
                                                   corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                                                   cohort: cohort)
    }

    public static func qualifyEvaluation(_ directory: URL, model modelDirectory: URL,
                                         corpus: VivoCellResponseCorpusReader,
                                         implementation: VivoFingerprint) throws -> VivoCellResponseEvaluation {
        try qualifyEvaluation(directory, model: modelDirectory,
                              corpus: VivoCellResponseCompositeCorpusReader(readers: [corpus]),
                              implementation: implementation)
    }
}
