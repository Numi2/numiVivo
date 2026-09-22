import Foundation
import MLX
import Testing
@testable import NumiVivoKit
@testable import NumiVivoLearning

@Suite("Cell response learning")
struct VivoCellResponseLearningTests {
    private static let implementation: VivoFingerprint = try! .init(bytes: Array(repeating: 7, count: 32))

    private static func sample(_ id: String, condition: String? = nil, batchID: String = "synthetic-batch",
                               replicateID: String? = nil) -> VivoOmicsSample {
        .init(id: id, biologicalReplicateID: replicateID ?? id + "-replicate", condition: condition ?? id,
              batchID: batchID, organism: "human")
    }

    /// This fixture is explicitly synthetic and only checks software behavior.
    /// It must not be interpreted as predictive or biological evidence.
    private static func dataset() -> VivoSingleCellDataset {
        let sampleIDs = [
            "control-training", "target-training", "target-b-training",
            "control-validation", "target-validation", "target-a-validation",
            "control-test", "target-a-test", "target-b-test", "empty-metadata-sample"
        ]
        let sourceSamples = [
            sample(sampleIDs[0], condition: sampleIDs[0], batchID: "context-training", replicateID: "pair-training"),
            sample(sampleIDs[1], condition: "target-a", batchID: "context-training", replicateID: "pair-training"),
            sample(sampleIDs[2], condition: "target-b", batchID: "context-training", replicateID: "pair-training"),
            sample(sampleIDs[3], condition: sampleIDs[3], batchID: "context-validation", replicateID: "pair-validation"),
            sample(sampleIDs[4], condition: "target-b", batchID: "context-validation", replicateID: "pair-validation"),
            sample(sampleIDs[5], condition: "target-a", batchID: "context-validation", replicateID: "pair-validation"),
            sample(sampleIDs[6], condition: sampleIDs[6], batchID: "context-test", replicateID: "pair-test"),
            sample(sampleIDs[7], condition: "target-a", batchID: "context-test", replicateID: "pair-test"),
            sample(sampleIDs[8], condition: "target-b", batchID: "context-test", replicateID: "pair-test"),
            sample(sampleIDs[9])
        ]
        let cells = [
            VivoOmicsCell(barcode: "ct0", sampleID: sampleIDs[0]), VivoOmicsCell(barcode: "ct1", sampleID: sampleIDs[0]),
            VivoOmicsCell(barcode: "tt0", sampleID: sampleIDs[1]), VivoOmicsCell(barcode: "tt1", sampleID: sampleIDs[1]),
            VivoOmicsCell(barcode: "tbtr0", sampleID: sampleIDs[2]), VivoOmicsCell(barcode: "tbtr1", sampleID: sampleIDs[2]),
            VivoOmicsCell(barcode: "cv0", sampleID: sampleIDs[3]), VivoOmicsCell(barcode: "cv1", sampleID: sampleIDs[3]),
            VivoOmicsCell(barcode: "tv0", sampleID: sampleIDs[4]), VivoOmicsCell(barcode: "tv1", sampleID: sampleIDs[4]),
            VivoOmicsCell(barcode: "tav0", sampleID: sampleIDs[5]), VivoOmicsCell(barcode: "tav1", sampleID: sampleIDs[5]),
            VivoOmicsCell(barcode: "cte0", sampleID: sampleIDs[6]), VivoOmicsCell(barcode: "cte1", sampleID: sampleIDs[6]),
            VivoOmicsCell(barcode: "tat0", sampleID: sampleIDs[7]), VivoOmicsCell(barcode: "tat1", sampleID: sampleIDs[7]),
            VivoOmicsCell(barcode: "tbt0", sampleID: sampleIDs[8]), VivoOmicsCell(barcode: "tbt1", sampleID: sampleIDs[8])
        ]
        // The first row deliberately omits gene-b. It is a measured zero,
        // while gene-absent is structurally absent from the source axis.
        let rows: [[(Int, UInt64)]] = [
            [(0, 8)], [(0, 7), (1, 1)],
            [(0, 1), (1, 8)], [(0, 2), (1, 7)],
            [(0, 2), (1, 7)], [(0, 1), (1, 8)],
            [(0, 8), (1, 2)], [(0, 7), (1, 2)],
            [(0, 2), (1, 8)], [(0, 1), (1, 7)],
            [(0, 1), (1, 8)], [(0, 2), (1, 7)],
            [(0, 8), (1, 1)], [(0, 7), (1, 2)],
            [(0, 1), (1, 8)], [(0, 2), (1, 7)],
            [(0, 2), (1, 8)], [(0, 1), (1, 7)]
        ]
        var offsets = [0], indices: [Int] = [], counts: [UInt64] = []
        for row in rows {
            indices.append(contentsOf: row.map(\.0)); counts.append(contentsOf: row.map(\.1)); offsets.append(counts.count)
        }
        // The final declared sample has zero cells. It models selection plans
        // that retain an empty source stratum; it must not require an invented
        // corpus role, while every observed cell remains explicitly assigned.
        return .init(id: "synthetic-cell-response", evidence: .synthetic,
                     sourceDescription: "Synthetic software-only response fixture", countUnit: .umiCount,
                     samples: sourceSamples,
                     features: [.init(id: "gene-a", name: "Gene A"), .init(id: "gene-b", name: "Gene B")],
                     cells: cells, matrix: .init(cellCount: cells.count, featureCount: 2, rowOffsets: offsets,
                                                  featureIndices: indices, counts: counts))
    }

    /// This fixture has one declared biological unit across distinct matched
    /// sample IDs. It exists only to prove that source-scoped split leakage is
    /// rejected even when plans and study labels differ.
    private static func sharedBiologicalUnitDataset() -> VivoSingleCellDataset {
        let base = dataset()
        let samples = base.samples.map { sample in
            VivoOmicsSample(id: sample.id, biologicalReplicateID: "shared-biological-unit",
                             donorID: sample.donorID, condition: sample.condition,
                             batchID: sample.batchID, organism: sample.organism)
        }
        return .init(id: "synthetic-cell-response-shared-biological-unit", evidence: .synthetic,
                     sourceDescription: "Synthetic source-split leakage fixture", countUnit: .umiCount,
                     samples: samples, features: base.features, cells: base.cells, matrix: base.matrix)
    }

    /// Wide synthetic output axis used only to catch decoder-gradient dilution.
    /// The added genes are measured zeros, not structural absences.
    private static func wideDataset(featureCount: Int = 4_096) -> VivoSingleCellDataset {
        let base = dataset()
        let additional = (2..<featureCount).map { index in
            VivoOmicsFeature(id: "gene-wide-\(index)", name: "gene-wide-\(index)")
        }
        return .init(id: "synthetic-cell-response-wide", evidence: .synthetic,
                     sourceDescription: "Synthetic wide-output software-only response fixture", countUnit: .umiCount,
                     samples: base.samples, features: base.features + additional, cells: base.cells,
                     matrix: .init(cellCount: base.matrix.cellCount, featureCount: featureCount,
                                   rowOffsets: base.matrix.rowOffsets, featureIndices: base.matrix.featureIndices,
                                   counts: base.matrix.counts))
    }

    private static func importPlan(samples: [VivoOmicsSample]? = nil) -> VivoH5ADImportPlan {
        .init(id: "synthetic-import", evidence: .synthetic,
              sourceDescription: "Synthetic software-only response fixture", countUnit: .umiCount,
              matrixPath: "X", samples: samples ?? dataset().samples, sampleColumn: "sample", barcodeColumn: "barcode",
              groupColumn: "group", featureNameColumn: "name")
    }

    private static let targets = [
        VivoCellResponseTarget(id: "target-a", descriptors: [0.25, -0.5]),
        VivoCellResponseTarget(id: "target-b", descriptors: [0.6, 0.4])
    ]
    private static let descriptorSource = VivoCellResponseDescriptorSource(
        sourceDescription: "Synthetic software-only descriptor fixture", targets: targets, sourceArtifacts: [implementation])
    private static let descriptorSourceBytes = try! VivoCanonicalJSON.encode(descriptorSource)
    private static let descriptorSourceFingerprint = try! VivoCanonicalJSON.fingerprint(descriptorSourceBytes)

    private static func corpusPlan(featureIDs: [String] = ["gene-a", "gene-b", "gene-absent"],
                                   allStrataTraining: Bool = false) -> VivoCellResponseCorpusPlan {
        let base: (String, VivoCellResponseRole, VivoCellResponsePartition, String, String, String?) -> VivoCellResponseAssignment = {
            sample, role, partition, context, pair, targetID in
            let guideID = role == .perturbed ? targetID! : sample
            return VivoCellResponseAssignment(
                sampleID: sample, sourceSample: Self.sample(sample, condition: guideID, batchID: context, replicateID: pair),
                role: role, targetID: role == .perturbed ? targetID : nil,
                guideID: guideID, modality: "rna", studyID: "synthetic-study",
                contextID: context, pairID: pair, partition: partition)
        }
        return .init(id: "synthetic-corpus", sourceDescription: "Synthetic software-only corpus",
                     featureAxis: .init(featureIDs: featureIDs), targets: targets,
                     descriptorSource: descriptorSourceFingerprint,
                     assignments: [
                        base("control-training", .control, .training, "context-training", "pair-training", nil),
                        base("target-training", .perturbed, .training, "context-training", "pair-training", "target-a"),
                        base("target-b-training", .perturbed, .training, "context-training", "pair-training", "target-b"),
                        base("control-validation", .control, allStrataTraining ? .training : .validation,
                             "context-validation", "pair-validation", nil),
                        base("target-validation", .perturbed, allStrataTraining ? .training : .validation,
                             "context-validation", "pair-validation", "target-b"),
                        base("target-a-validation", .perturbed, allStrataTraining ? .training : .validation,
                             "context-validation", "pair-validation", "target-a"),
                        base("control-test", .control, allStrataTraining ? .training : .test,
                             "context-test", "pair-test", nil),
                        base("target-a-test", .perturbed, allStrataTraining ? .training : .test,
                             "context-test", "pair-test", "target-a"),
                        base("target-b-test", .perturbed, allStrataTraining ? .training : .test,
                             "context-test", "pair-test", "target-b")
                     ], heldOutContextIDs: allStrataTraining ? [] : ["context-validation", "context-test"])
    }

    /// A hermetic schema-v4 artifact exercises real v5 receipt, checkpoint,
    /// safetensors, and state loading without relying on a mutable external
    /// model directory. It deliberately has no cohort admission.
    private static func writeLegacyV5Model(to destination: URL) throws -> VivoCellResponseModelReceipt {
        let architecture = VivoCellResponseArchitecture(hiddenWidth: 16, contextCells: 1,
                                                         maximumParameters: 100_000)
        let plan = VivoCellResponseTrainingPlan(id: "legacy-v5-plan", architecture: architecture,
                                                 steps: 1, batchSize: 1, learningRate: 0.01,
                                                 weightDecay: 0, seed: 71, validationEvery: 1,
                                                 validationExamples: 1)
        let planBytes = try VivoCanonicalJSON.encode(plan)
        let corpus = try VivoFingerprint(bytes: Array(repeating: 23, count: 32))
        let target = VivoCellResponseTarget(id: "legacy-target", descriptors: [0.25])
        let state = VivoCellResponseModelState(
            schemaVersion: 4, architecture: architecture,
            featureAxis: .init(featureIDs: ["legacy-gene"]), targetIDs: [target.id],
            trainedTargetBindings: [.init(id: target.id, descriptorFingerprint: try target.fingerprint())],
            descriptorCount: 1, corpus: corpus,
            trainingPlan: try VivoCanonicalJSON.fingerprint(planBytes), step: 1, seed: plan.seed,
            samplerVersion: VivoCellResponseLearning.samplerVersion,
            targetResponsePriorVersion: VivoCellResponseLearning.targetResponsePriorVersion,
            optimizerVersion: VivoCellResponseLearning.optimizerVersion,
            learningRate: plan.learningRate, weightDecay: plan.weightDecay,
            latestMetrics: .init(step: 1, trainNegativeLogLikelihood: 0.5,
                                 validationNegativeLogLikelihood: nil, validationRMSE: nil,
                                 validationMatchedControlRMSE: nil))
        let model = VivoCellResponseMLXModel(featureCount: 1, targetCount: 1, descriptorCount: 1,
                                              architecture: architecture)
        let weights = try MLX.saveToData(
            arrays: Dictionary(uniqueKeysWithValues: model.parameters().flattened()),
            metadata: [
                "format": VivoCellResponseLearning.legacyModelFormat,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1",
                "featureCount": "1",
                "step": "1"
            ])
        let checkpoint = try VivoCheckpointCodec.encode(.init(
            runtime: .cellResponseMLX, artifactFingerprint: corpus.hex, sourceFingerprint: corpus.hex,
            experimentFingerprint: state.trainingPlan.hex, stepIndex: state.step,
            logicalTime: Double(state.step), randomStreamVersion: state.samplerVersion,
            sections: [
                .init(id: "weights.safetensors", encoding: .rawBytes,
                      elementCount: UInt64(weights.count), elementStride: 1, data: weights),
                try .canonicalJSON(id: "state.json", value: state)
            ], metadata: [
                "format": VivoCellResponseLearning.legacyModelFormat,
                "optimizer": VivoCellResponseLearning.optimizerFormat,
                "sampler": "splitmix64-v1",
                "targetResponsePrior": "balanced-training-mean-v1"
            ]))
        let receipt = VivoCellResponseModelReceipt(
            schemaVersion: 1, format: VivoCellResponseLearning.legacyModelFormat,
            checkpoint: try VivoCanonicalJSON.fingerprint(checkpoint),
            trainingPlan: try VivoCanonicalJSON.fingerprint(planBytes), corpus: corpus,
            implementation: implementation)
        try FileManager.default.createDirectory(at: destination, withIntermediateDirectories: false)
        try planBytes.write(to: destination.appendingPathComponent("training-plan.json"), options: .withoutOverwriting)
        try checkpoint.write(to: destination.appendingPathComponent("checkpoint.nvckpt"), options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(receipt).write(to: destination.appendingPathComponent("receipt.json"),
                                                     options: .withoutOverwriting)
        return receipt
    }

    @Test("source-bound corpus trains, resumes, predicts, and rejects tampering")
    func sourceBoundLifecycle() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let h5ad = root.appendingPathComponent("synthetic.h5ad")
        let store = root.appendingPathComponent("count-store")
        let corpusDirectory = root.appendingPathComponent("corpus")
        let model = root.appendingPathComponent("model")
        let resumed = root.appendingPathComponent("resumed")
        let straightThrough = root.appendingPathComponent("straight-through")
        let prediction = root.appendingPathComponent("prediction")
        let renamedPrediction = root.appendingPathComponent("renamed-prediction")
        let straightPrediction = root.appendingPathComponent("straight-prediction")
        let knownPrediction = root.appendingPathComponent("known-prediction")
        let straightKnownPrediction = root.appendingPathComponent("straight-known-prediction")
        let persistedEvaluation = root.appendingPathComponent("persisted-evaluation")
        let foreignCorpusDirectory = root.appendingPathComponent("foreign-corpus")
        let foreignPrediction = root.appendingPathComponent("foreign-prediction")
        let encodedCorpusPlan = try VivoCanonicalJSON.encode(Self.corpusPlan())
        #expect(try VivoCanonicalJSON.decode(VivoCellResponseCorpusPlan.self, from: encodedCorpusPlan) == Self.corpusPlan())
        try VivoSingleCellH5AD.write(Self.dataset(), to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: Self.importPlan(), implementation: Self.implementation, to: store)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store, plan: Self.corpusPlan(), descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: corpusDirectory)
        let reader = try VivoCellResponseCorpus.open(corpusDirectory, sourceStore: store, implementation: Self.implementation)
        #expect(reader.featureMask == [1, 1, 0])
        #expect(try reader.trainingTargetBindings().map(\.id) == ["target-a", "target-b"])
        let first = try reader.normalizedRow(0)
        #expect(first[1] == 0)
        #expect(first[2] == -1)

        let plan = VivoCellResponseTrainingPlan(id: "synthetic-train",
                                                 architecture: .init(hiddenWidth: 16, contextCells: 2, maximumParameters: 100_000),
                                                 steps: 12, batchSize: 1, learningRate: 0.01, weightDecay: 0,
                                                 seed: 17, validationEvery: 4, validationExamples: 2)
        #expect(try VivoCellResponseLearning.trainingControlCacheMatchesUncachedBatchForTesting(
            corpus: reader, plan: plan, step: 1))
        let cohort = try VivoCellResponseCohort.admit(
            readers: [reader],
            requirements: .init(minimumTrainingBiologicalUnitsPerTarget: 1,
                                minimumValidationBiologicalUnitsPerTarget: 1,
                                minimumTestBiologicalUnitsPerTarget: 1,
                                minimumBiologicalUnits: 3))
        try VivoCellResponseCohort.verify(cohort, readers: [reader])
        let cohortFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(cohort))
        let receipt = try VivoCellResponseLearning.train(corpus: reader, plan: plan, cohort: cohort,
                                                          implementation: Self.implementation, to: model)
        #expect(receipt.format == VivoCellResponseLearning.modelFormat)
        #expect(receipt.cohort == cohortFingerprint)
        _ = try VivoCellResponseLearning.verifyModel(model, implementation: Self.implementation)
        let evaluation = try VivoCellResponseLearning.evaluate(model: model, corpus: reader, partition: .validation,
                                                                maximumExamples: 2, implementation: Self.implementation)
        #expect(evaluation.examples > 0)
        #expect(evaluation.observedFeatures > 0)
        #expect(evaluation.negativeLogLikelihood.isFinite)
        #expect(evaluation.rmse.isFinite)
        #expect(evaluation.matchedControlRMSE.isFinite)
        #expect(evaluation.selection.count == evaluation.examples)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireQualificationEvaluationCoverage(evaluation, corpus: reader)
        }
        let incompleteTestEvaluation = try VivoCellResponseLearning.evaluate(
            model: model, corpus: reader, partition: .test, maximumExamples: 1,
            implementation: Self.implementation)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireQualificationEvaluationCoverage(incompleteTestEvaluation, corpus: reader)
        }
        let exhaustiveTestEvaluation = try VivoCellResponseLearning.evaluate(
            model: model, corpus: reader, partition: .test, maximumExamples: 4,
            implementation: Self.implementation)
        try VivoCellResponseLearning.requireQualificationEvaluationCoverage(exhaustiveTestEvaluation, corpus: reader)
        let evaluationReceipt = try VivoCellResponseLearning.evaluate(model: model, corpus: reader, partition: .validation,
                                                                       maximumExamples: 2, implementation: Self.implementation,
                                                                       to: persistedEvaluation)
        #expect(evaluationReceipt.format == VivoCellResponseLearning.evaluationFormat)
        #expect(try VivoCellResponseLearning.verifyEvaluation(persistedEvaluation, model: model, corpus: reader,
                                                              implementation: Self.implementation) == evaluation)

        let legacy = root.appendingPathComponent("legacy-v5-model")
        let legacyReceipt = try Self.writeLegacyV5Model(to: legacy)
        #expect(legacyReceipt.cohort == nil)
        #expect(try VivoCellResponseLearning.verifyModel(legacy, implementation: Self.implementation) == legacyReceipt)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.resume(model: legacy, corpus: reader,
                                                 plan: .init(id: "legacy-resume", additionalSteps: 1),
                                                 implementation: Self.implementation,
                                                 to: root.appendingPathComponent("legacy-resumed"))
        }
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.qualifyEvaluation(persistedEvaluation, model: legacy,
                                                            corpus: reader, implementation: Self.implementation)
        }

        let resumedReceipt = try VivoCellResponseLearning.resume(model: model, corpus: reader,
                                                                  plan: .init(id: "synthetic-resume", additionalSteps: 3),
                                                                  implementation: Self.implementation, to: resumed)
        #expect(resumedReceipt.checkpoint != receipt.checkpoint)
        _ = try VivoCellResponseLearning.verifyModel(resumed, implementation: Self.implementation)
        let predictionPlan = VivoCellResponsePredictionPlan(
            id: "novel-target", target: .init(id: "novel-target", descriptors: [0.1, 0.2]),
            contextSampleIDs: ["control-validation"], useTrainedTargetEmbedding: false)
        let knownTargetPlan = VivoCellResponsePredictionPlan(
            id: "known-target", target: .init(id: "target-a", descriptors: [0.25, -0.5]),
            contextSampleIDs: ["control-training"], useTrainedTargetEmbedding: true)

        let originalCorpusPlan = Self.corpusPlan()
        let foreignPlan = VivoCellResponseCorpusPlan(
            id: "synthetic-foreign-corpus", sourceDescription: "Distinct synthetic software-only corpus",
            featureAxis: originalCorpusPlan.featureAxis, targets: originalCorpusPlan.targets,
            descriptorSource: originalCorpusPlan.descriptorSource,
            assignments: originalCorpusPlan.assignments,
            heldOutTargetIDs: originalCorpusPlan.heldOutTargetIDs,
            heldOutContextIDs: originalCorpusPlan.heldOutContextIDs)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store, plan: foreignPlan, descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: foreignCorpusDirectory)
        let foreignReader = try VivoCellResponseCorpus.open(foreignCorpusDirectory, sourceStore: store,
                                                             implementation: Self.implementation)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.predict(model: resumed, corpus: foreignReader, plan: predictionPlan,
                                                  implementation: Self.implementation, to: foreignPrediction)
        }

        let mismatchedKnownTarget = VivoCellResponsePredictionPlan(
            id: "mismatched-known-target", target: .init(id: "target-a", descriptors: [0.5, -0.5]),
            contextSampleIDs: ["control-validation"], useTrainedTargetEmbedding: true)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: mismatchedKnownTarget,
                                                  implementation: Self.implementation,
                                                  to: root.appendingPathComponent("mismatched-known-target"))
        }

        let descriptorOnlyKnownTarget = VivoCellResponsePredictionPlan(
            id: "descriptor-only-known-target", target: .init(id: "target-b", descriptors: [0.6, 0.4]),
            contextSampleIDs: ["control-validation"], useTrainedTargetEmbedding: false)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: descriptorOnlyKnownTarget,
                                                  implementation: Self.implementation,
                                                  to: root.appendingPathComponent("descriptor-only-known-target"))

        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: predictionPlan,
                                                  implementation: Self.implementation, to: prediction)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: knownTargetPlan,
                                                  implementation: Self.implementation, to: knownPrediction)
        let renamedPlan = VivoCellResponsePredictionPlan(id: "novel-target-renamed", target: predictionPlan.target,
                                                          contextSampleIDs: predictionPlan.contextSampleIDs)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: renamedPlan,
                                                  implementation: Self.implementation, to: renamedPrediction)
        let result = try VivoCellResponseLearning.verifyPrediction(prediction, implementation: Self.implementation)
        let renamedResult = try VivoCellResponseLearning.verifyPrediction(renamedPrediction, implementation: Self.implementation)
        let boundResult = try VivoCellResponseLearning.verifyPrediction(prediction, model: resumed, corpus: reader,
                                                                         implementation: Self.implementation)
        #expect(result.featureMask == [1, 1, 0])
        #expect(boundResult == result)
        #expect(result.contextSourceRows.count == 2)
        #expect(renamedResult.contextSourceRows == result.contextSourceRows)
        #expect(zip(renamedResult.meanLogCPM, result.meanLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(zip(result.meanLogCPM, result.contextBaselineLogCPM), result.meanDeltaLogCPM).allSatisfy {
            abs($0.0.0 - ($0.0.1 + $0.1)) <= 1e-5
        })
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.verifyPrediction(prediction, model: resumed, corpus: foreignReader,
                                                           implementation: Self.implementation)
        }
        #expect(result.meanLogCPM.count == 3 && result.varianceLogCPM.allSatisfy { $0 > 0 })

        let uninterruptedPlan = VivoCellResponseTrainingPlan(
            id: "synthetic-uninterrupted", architecture: plan.architecture,
            steps: 15, batchSize: plan.batchSize, learningRate: plan.learningRate, weightDecay: plan.weightDecay,
            seed: plan.seed, validationEvery: plan.validationEvery, validationExamples: plan.validationExamples)
        _ = try VivoCellResponseLearning.train(corpus: reader, plan: uninterruptedPlan,
                                                cohort: cohort, implementation: Self.implementation, to: straightThrough)
        _ = try VivoCellResponseLearning.predict(model: straightThrough, corpus: reader, plan: predictionPlan,
                                                  implementation: Self.implementation, to: straightPrediction)
        _ = try VivoCellResponseLearning.predict(model: straightThrough, corpus: reader, plan: knownTargetPlan,
                                                  implementation: Self.implementation, to: straightKnownPrediction)
        let straightResult = try VivoCellResponseLearning.verifyPrediction(straightPrediction, implementation: Self.implementation)
        let knownResult = try VivoCellResponseLearning.verifyPrediction(knownPrediction, implementation: Self.implementation)
        let straightKnownResult = try VivoCellResponseLearning.verifyPrediction(straightKnownPrediction,
                                                                                  implementation: Self.implementation)
        #expect(zip(result.meanLogCPM, straightResult.meanLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(result.varianceLogCPM, straightResult.varianceLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        // This path exercises the frozen training-only target response prior.
        // Equality after a checkpoint/resume proves its values were serialized
        // and restored rather than silently replaced with zeros.
        #expect(zip(knownResult.meanLogCPM, straightKnownResult.meanLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(knownResult.meanDeltaLogCPM, straightKnownResult.meanDeltaLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(knownResult.varianceLogCPM, straightKnownResult.varianceLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })

        var checkpoint = try Data(contentsOf: resumed.appendingPathComponent("checkpoint.nvckpt"))
        checkpoint[0] ^= 0x01
        try checkpoint.write(to: resumed.appendingPathComponent("checkpoint.nvckpt"))
        #expect(throws: (any Error).self) { try VivoCellResponseLearning.verifyModel(resumed, implementation: Self.implementation) }
    }

    @Test("plan rejects validation treated cells matched to training controls")
    func splitLeakageIsRejected() throws {
        let target = VivoCellResponseTarget(id: "target-a", descriptors: [1])
        let assignments = [
            VivoCellResponseAssignment(sampleID: "control", sourceSample: Self.sample("control", condition: "control", batchID: "context"), role: .control, guideID: "control", modality: "rna", studyID: "study",
                                        contextID: "context", pairID: "pair", partition: .training),
            VivoCellResponseAssignment(sampleID: "treated", sourceSample: Self.sample("treated", condition: "target-a", batchID: "context"), role: .perturbed, targetID: "target-a", guideID: "target-a",
                                        modality: "rna", studyID: "study", contextID: "context", pairID: "pair", partition: .validation)
        ]
        let plan = VivoCellResponseCorpusPlan(id: "leakage", sourceDescription: "Synthetic split guard",
                                               featureAxis: .init(featureIDs: ["gene-a"]), targets: [target],
                                               descriptorSource: try! VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
                                                    VivoCellResponseDescriptorSource(sourceDescription: "Synthetic split descriptor", targets: [target], sourceArtifacts: [Self.implementation]))),
                                               assignments: assignments)
        #expect(throws: (any Error).self) { try plan.validate() }
    }

    @Test("training-control cache replays source batches across strata and wraparound")
    func trainingControlCacheReplaysAcrossStrata() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-cache-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let h5ad = root.appendingPathComponent("synthetic.h5ad")
        let store = root.appendingPathComponent("count-store")
        let corpusDirectory = root.appendingPathComponent("corpus")
        try VivoSingleCellH5AD.write(Self.dataset(), to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: Self.importPlan(), implementation: Self.implementation, to: store)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store, plan: Self.corpusPlan(allStrataTraining: true),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: corpusDirectory)
        let reader = try VivoCellResponseCorpus.open(corpusDirectory, sourceStore: store, implementation: Self.implementation)
        let plan = VivoCellResponseTrainingPlan(id: "cache-equivalence",
                                                 architecture: .init(hiddenWidth: 16, contextCells: 3, maximumParameters: 100_000),
                                                 steps: 3, batchSize: 2, learningRate: 0.01, weightDecay: 0,
                                                 seed: 31, validationEvery: 3, validationExamples: 1)
        #expect(try VivoCellResponseLearning.trainingControlCacheMatchesUncachedBatchForTesting(
            corpus: reader, plan: plan, step: 1))
        #expect(try VivoCellResponseLearning.trainingControlCacheMatchesUncachedBatchForTesting(
            corpus: reader, plan: plan, step: 2))
    }

    @Test("cohort admission replays receipt-bound treated and control rows")
    func cohortAdmissionRejectsTamperedRows() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-cohort-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let h5ad = root.appendingPathComponent("synthetic.h5ad")
        let store = root.appendingPathComponent("count-store")
        let corpusDirectory = root.appendingPathComponent("corpus")
        try VivoSingleCellH5AD.write(Self.dataset(), to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: Self.importPlan(), implementation: Self.implementation, to: store)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store, plan: Self.corpusPlan(),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: corpusDirectory)
        let reader = try VivoCellResponseCorpus.open(corpusDirectory, sourceStore: store, implementation: Self.implementation)
        let requirements = VivoCellResponseCohortRequirements(minimumTrainingBiologicalUnitsPerTarget: 1,
                                                                minimumValidationBiologicalUnitsPerTarget: 1,
                                                                minimumTestBiologicalUnitsPerTarget: 1,
                                                                minimumBiologicalUnits: 3)
        let admitted = try VivoCellResponseCohort.admit(readers: [reader], requirements: requirements)
        #expect(admitted.contexts.count == 3)
        #expect(admitted.coverage.count == 6)
        #expect(admitted.coverage.allSatisfy { $0.biologicalUnits == 1 && $0.matchedContexts == 1 })
        #expect(try VivoCanonicalJSON.decode(VivoCellResponseCohortAdmission.self,
                                              from: VivoCanonicalJSON.encode(admitted)) == admitted)
        try VivoCellResponseCohort.verify(admitted, readers: [reader])
        #expect(throws: (any Error).self) {
            try VivoCellResponseCohort.requireQualificationCoverage(admitted)
        }

        var contexts = admitted.contexts
        let original = contexts[0]
        let target = original.treatedSamples[0]
        let substituted = target.targetID == "target-a" ? "target-b-training" : "target-training"
        let alteredTarget = VivoCellResponseCohortTargetSamples(targetID: target.targetID, sampleIDs: [substituted])
        contexts[0] = .init(corpus: original.corpus, studyID: original.studyID,
                            biologicalReplicateID: original.biologicalReplicateID, donorID: original.donorID,
                            contextID: original.contextID, pairID: original.pairID, modality: original.modality,
                            partition: original.partition,
                            treatedSamples: [alteredTarget] + Array(original.treatedSamples.dropFirst()),
                            controlSampleIDs: original.controlSampleIDs)
        let tampered = VivoCellResponseCohortAdmission(requirements: admitted.requirements,
                                                        sources: admitted.sources, featureAxis: admitted.featureAxis,
                                                        descriptorSource: admitted.descriptorSource,
                                                        targets: admitted.targets, contexts: contexts,
                                                        coverage: admitted.coverage)
        try VivoCellResponseCohort.validate(tampered)
        #expect(throws: (any Error).self) {
            try VivoCellResponseCohort.verify(tampered, readers: [reader])
        }

        var mixedDonorContexts = admitted.contexts
        let trainingIndex = mixedDonorContexts.firstIndex(where: { $0.partition == .training })!
        let testIndex = mixedDonorContexts.firstIndex(where: { $0.partition == .test })!
        let training = mixedDonorContexts[trainingIndex]
        mixedDonorContexts[trainingIndex] = .init(
            corpus: training.corpus, studyID: training.studyID,
            biologicalReplicateID: training.biologicalReplicateID, donorID: "donor-shared",
            contextID: training.contextID, pairID: training.pairID, modality: training.modality,
            partition: training.partition, treatedSamples: training.treatedSamples,
            controlSampleIDs: training.controlSampleIDs)
        let test = mixedDonorContexts[testIndex]
        mixedDonorContexts[testIndex] = .init(
            corpus: test.corpus, studyID: test.studyID,
            biologicalReplicateID: training.biologicalReplicateID, donorID: nil,
            contextID: test.contextID, pairID: test.pairID, modality: test.modality,
            partition: test.partition, treatedSamples: test.treatedSamples,
            controlSampleIDs: test.controlSampleIDs)
        mixedDonorContexts.sort {
            let left = [$0.corpus.hex, $0.studyID, $0.biologicalReplicateID, $0.donorID ?? "",
                        $0.contextID, $0.pairID, $0.modality]
            let right = [$1.corpus.hex, $1.studyID, $1.biologicalReplicateID, $1.donorID ?? "",
                         $1.contextID, $1.pairID, $1.modality]
            return left.lexicographicallyPrecedes(right)
        }
        let mixedDonor = VivoCellResponseCohortAdmission(
            requirements: admitted.requirements, sources: admitted.sources,
            featureAxis: admitted.featureAxis, descriptorSource: admitted.descriptorSource,
            targets: admitted.targets, contexts: mixedDonorContexts, coverage: admitted.coverage)
        #expect(throws: (any Error).self) {
            try VivoCellResponseCohort.validate(mixedDonor)
        }
    }

    @Test("cohort rejects one raw-source unit across split-labelled studies")
    func cohortRejectsSourceSplitLeakageDespiteStudyRelabeling() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-cohort-source-split-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let dataset = Self.sharedBiologicalUnitDataset()
        let h5ad = root.appendingPathComponent("shared-unit.h5ad")
        let store = root.appendingPathComponent("count-store")
        let corpusDirectory = root.appendingPathComponent("corpus")
        try VivoSingleCellH5AD.write(dataset, to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: Self.importPlan(samples: dataset.samples),
                                            implementation: Self.implementation, to: store)
        let sourceSamples = Dictionary(uniqueKeysWithValues: dataset.samples.map { ($0.id, $0) })
        let template = Self.corpusPlan()
        let assignments = template.assignments.map { assignment in
            let studyID = "synthetic-study-\(assignment.partition.rawValue)-view"
            return VivoCellResponseAssignment(sampleID: assignment.sampleID,
                                               sourceSample: sourceSamples[assignment.sampleID]!,
                                               role: assignment.role, targetID: assignment.targetID,
                                               guideID: assignment.guideID, modality: assignment.modality,
                                               studyID: studyID, contextID: assignment.contextID,
                                               pairID: assignment.pairID, partition: assignment.partition)
        }
        let plan = VivoCellResponseCorpusPlan(
            id: "shared-unit-source-split", sourceDescription: "Synthetic source-split cohort fixture",
            featureAxis: template.featureAxis, targets: template.targets,
            descriptorSource: template.descriptorSource, assignments: assignments,
            heldOutContextIDs: ["context-validation", "context-test"])
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store, plan: plan,
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: corpusDirectory)
        let reader = try VivoCellResponseCorpus.open(corpusDirectory, sourceStore: store,
                                                      implementation: Self.implementation)
        #expect(Set(reader.plan.assignments.map(\.studyID)).count == 3)
        #expect(throws: (any Error).self) {
            try VivoCellResponseCohort.admit(
                readers: [reader],
                requirements: .init(minimumTrainingBiologicalUnitsPerTarget: 1,
                                    minimumValidationBiologicalUnitsPerTarget: 1,
                                    minimumTestBiologicalUnitsPerTarget: 1,
                                    minimumBiologicalUnits: 3))
        }
    }

    @Test("output-axis optimizer trains a wide synthetic response residual")
    func outputAxisOptimizerTrainsWideResidual() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-wide-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let dataset = Self.wideDataset()
        let h5ad = root.appendingPathComponent("wide-synthetic.h5ad")
        let store = root.appendingPathComponent("count-store")
        let corpusDirectory = root.appendingPathComponent("corpus")
        try VivoSingleCellH5AD.write(dataset, to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: Self.importPlan(), implementation: Self.implementation, to: store)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: store,
                                                plan: Self.corpusPlan(featureIDs: dataset.features.map(\.id)),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: corpusDirectory)
        let reader = try VivoCellResponseCorpus.open(corpusDirectory, sourceStore: store, implementation: Self.implementation)
        // 4,096 observed outputs make the legacy mean-loss normalization
        // shrink an independent decoder-row update below the assertion while
        // the output-axis optimizer retains the calibrated update scale.
        let architecture = VivoCellResponseArchitecture(hiddenWidth: 16, contextCells: 2, maximumParameters: 250_000)
        let plan = VivoCellResponseTrainingPlan(id: "wide-output-axis", architecture: architecture,
                                                 steps: 8, batchSize: 1, learningRate: 0.01, weightDecay: 0,
                                                 seed: 23, validationEvery: 8, validationExamples: 1)
        let maximumDeltaDifference = try VivoCellResponseLearning.outputAxisRouteDifferenceForTesting(
            corpus: reader, plan: plan)
        #expect(maximumDeltaDifference > 1e-4)
    }

    @Test("qualification requires a strict matched-control improvement")
    func qualificationGateRejectsTiesAndLosses() throws {
        func evaluation(rmse: Double, baseline: Double) -> VivoCellResponseEvaluation {
            .init(schemaVersion: 2, format: VivoCellResponseLearning.evaluationFormat,
                  partition: .validation, maximumExamples: 1,
                  samplerVersion: VivoCellResponseLearning.samplerVersion, seed: 0,
                  examples: 1, observedFeatures: 1, negativeLogLikelihood: 0,
                  rmse: rmse, matchedControlRMSE: baseline,
                  selection: [.init(targetRow: 0, contextRows: [1])])
        }

        try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 0.9, baseline: 1))
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 1, baseline: 1))
        }
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 1.1, baseline: 1))
        }
    }
}
