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

    /// A second independent source deliberately reuses local row and sample
    /// labels. Its biological-replicate declarations and counts differ, so a
    /// composite must qualify identities by corpus receipt rather than by a
    /// bare local row or sample string.
    private static func secondDataset() -> VivoSingleCellDataset {
        let base = dataset()
        let samples = base.samples.map { sample in
            VivoOmicsSample(id: sample.id, biologicalReplicateID: "second-" + sample.biologicalReplicateID,
                             donorID: sample.donorID, condition: sample.condition, batchID: sample.batchID,
                             organism: sample.organism)
        }
        var counts = base.matrix.counts
        counts[0] = 16
        return .init(id: "synthetic-cell-response-second-source", evidence: .synthetic,
                     sourceDescription: "Second synthetic software-only response fixture", countUnit: .umiCount,
                     samples: samples, features: base.features, cells: base.cells,
                     matrix: .init(cellCount: base.matrix.cellCount, featureCount: base.matrix.featureCount,
                                   rowOffsets: base.matrix.rowOffsets, featureIndices: base.matrix.featureIndices,
                                   counts: counts))
    }

    /// Independent source with the same declared prediction axis but a smaller
    /// measured panel. `gene-b` is structurally absent here, unlike a measured
    /// zero in `dataset()`, so a composite learner must mask it before the
    /// shared cell encoder and likelihood.
    private static func secondPanelDataset() -> VivoSingleCellDataset {
        let base = secondDataset()
        var offsets = [0], indices: [Int] = [], counts: [UInt64] = []
        for row in 0..<base.matrix.cellCount {
            let lower = base.matrix.rowOffsets[row]
            let upper = base.matrix.rowOffsets[row + 1]
            for offset in lower..<upper where base.matrix.featureIndices[offset] == 0 {
                indices.append(0)
                counts.append(base.matrix.counts[offset])
            }
            offsets.append(counts.count)
        }
        return .init(id: "synthetic-cell-response-second-panel", evidence: .synthetic,
                     sourceDescription: "Second synthetic source with a narrower measured panel",
                     countUnit: .umiCount, samples: base.samples, features: [base.features[0]], cells: base.cells,
                     matrix: .init(cellCount: base.matrix.cellCount, featureCount: 1, rowOffsets: offsets,
                                   featureIndices: indices, counts: counts))
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
                                   allStrataTraining: Bool = false,
                                   replicatePrefix: String = "") -> VivoCellResponseCorpusPlan {
        let base: (String, VivoCellResponseRole, VivoCellResponsePartition, String, String, String?) -> VivoCellResponseAssignment = {
            sample, role, partition, context, pair, targetID in
            let guideID = role == .perturbed ? targetID! : sample
            return VivoCellResponseAssignment(
                sampleID: sample, sourceSample: Self.sample(sample, condition: guideID, batchID: context,
                                                            replicateID: replicatePrefix + pair),
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

    @Test("qualification requires held-out coverage for every raw-source target")
    func qualificationRequiresHeldOutSourceTargetCoverage() throws {
        let plan = Self.corpusPlan()
        let corpusA: VivoFingerprint = try .init(bytes: Array(repeating: 1, count: 32))
        let corpusB: VivoFingerprint = try .init(bytes: Array(repeating: 2, count: 32))
        let rawSourceA: VivoFingerprint = try .init(bytes: Array(repeating: 3, count: 32))
        let rawSourceB: VivoFingerprint = try .init(bytes: Array(repeating: 4, count: 32))
        let planA: VivoFingerprint = try .init(bytes: Array(repeating: 5, count: 32))
        let planB: VivoFingerprint = try .init(bytes: Array(repeating: 6, count: 32))
        let sources = [
            VivoCellResponseCohortSource(corpus: corpusA, plan: planA, source: rawSourceA),
            VivoCellResponseCohortSource(corpus: corpusB, plan: planB, source: rawSourceB)
        ]
        let developmentRequirements = VivoCellResponseCohortRequirements(
            minimumTrainingBiologicalUnitsPerTarget: 1,
            minimumValidationBiologicalUnitsPerTarget: 1,
            minimumTestBiologicalUnitsPerTarget: 1,
            minimumBiologicalUnits: 3)

        func context(corpus: VivoFingerprint, unit: Int, partition: VivoCellResponsePartition,
                     targetIDs: [String]) -> VivoCellResponseCohortContext {
            let suffix = String(format: "%02d", unit)
            return .init(
                corpus: corpus, studyID: "study", biologicalReplicateID: "unit-\(suffix)", donorID: nil,
                contextID: "context-\(suffix)", pairID: "pair-\(suffix)", modality: "rna",
                partition: partition,
                treatedSamples: targetIDs.map {
                    .init(targetID: $0, sampleIDs: ["treated-\($0)-\(suffix)"])
                },
                controlSampleIDs: ["control-\(suffix)"])
        }

        func admission(includeSourceATargetBTest: Bool) -> VivoCellResponseCohortAdmission {
            var contexts: [VivoCellResponseCohortContext] = []
            for unit in 0..<8 {
                contexts.append(context(corpus: corpusA, unit: unit, partition: .training,
                                        targetIDs: ["target-a", "target-b"]))
            }
            for unit in 8..<10 {
                contexts.append(context(corpus: corpusA, unit: unit, partition: .validation,
                                        targetIDs: ["target-a", "target-b"]))
            }
            for unit in 10..<12 {
                contexts.append(context(corpus: corpusA, unit: unit, partition: .test,
                                        targetIDs: includeSourceATargetBTest ? ["target-a", "target-b"] : ["target-a"]))
            }
            for unit in 12..<14 {
                contexts.append(context(corpus: corpusB, unit: unit, partition: .test,
                                        targetIDs: ["target-b"]))
            }
            let targetBTestCount = includeSourceATargetBTest ? 4 : 2
            let coverage = [
                VivoCellResponseCohortCoverage(targetID: "target-a", partition: .training,
                                                biologicalUnits: 8, matchedContexts: 8),
                VivoCellResponseCohortCoverage(targetID: "target-a", partition: .validation,
                                                biologicalUnits: 2, matchedContexts: 2),
                VivoCellResponseCohortCoverage(targetID: "target-a", partition: .test,
                                                biologicalUnits: 2, matchedContexts: 2),
                VivoCellResponseCohortCoverage(targetID: "target-b", partition: .training,
                                                biologicalUnits: 8, matchedContexts: 8),
                VivoCellResponseCohortCoverage(targetID: "target-b", partition: .validation,
                                                biologicalUnits: 2, matchedContexts: 2),
                VivoCellResponseCohortCoverage(targetID: "target-b", partition: .test,
                                                biologicalUnits: targetBTestCount, matchedContexts: targetBTestCount)
            ]
            return .init(requirements: developmentRequirements, sources: sources,
                         featureAxis: plan.featureAxis, descriptorSource: plan.descriptorSource,
                         targets: plan.targets, contexts: contexts, coverage: coverage)
        }

        let missingSourceATargetBTest = admission(includeSourceATargetBTest: false)
        try VivoCellResponseCohort.validate(missingSourceATargetBTest)
        #expect(throws: (any Error).self) {
            try VivoCellResponseCohort.requireQualificationCoverage(missingSourceATargetBTest)
        }

        let covered = admission(includeSourceATargetBTest: true)
        try VivoCellResponseCohort.validate(covered)
        try VivoCellResponseCohort.requireQualificationCoverage(covered)
    }

    @Test("source-balanced loss averages panel-normalized source lanes")
    func sourceBalancedLossWeightsPanelsEqually() {
        let terms = MLXArray([Float(2), 2, 777, 6, 777, 777], [2, 3])
        let mask = MLXArray([Float(1), 1, 0, 1, 0, 0], [2, 3])
        let values = VivoCellResponseLearning.sourceBalancedMaskedMean(terms, mask: mask).asArray(Float.self)
        #expect(values.count == 1)
        #expect(abs(values[0] - 4) <= 1e-6)
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

    @Test("composite corpus canonicalizes sources and qualifies colliding rows")
    func compositeCorpusCanonicalizesSourcesAndReplaysCohort() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cell-response-composite-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let first = Self.dataset(), second = Self.secondDataset()
        let firstH5AD = root.appendingPathComponent("first.h5ad")
        let secondH5AD = root.appendingPathComponent("second.h5ad")
        let firstStore = root.appendingPathComponent("first-store")
        let secondStore = root.appendingPathComponent("second-store")
        let firstCorpus = root.appendingPathComponent("first-corpus")
        let secondCorpus = root.appendingPathComponent("second-corpus")
        try VivoSingleCellH5AD.write(first, to: firstH5AD)
        try VivoSingleCellH5AD.write(second, to: secondH5AD)
        _ = try VivoH5ADCountStore.publish(source: firstH5AD, plan: Self.importPlan(samples: first.samples),
                                            implementation: Self.implementation, to: firstStore)
        _ = try VivoH5ADCountStore.publish(source: secondH5AD, plan: Self.importPlan(samples: second.samples),
                                            implementation: Self.implementation, to: secondStore)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: firstStore, plan: Self.corpusPlan(),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: firstCorpus)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: secondStore,
                                                plan: Self.corpusPlan(replicatePrefix: "second-"),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: secondCorpus)
        let firstReader = try VivoCellResponseCorpus.open(firstCorpus, sourceStore: firstStore,
                                                           implementation: Self.implementation)
        let secondReader = try VivoCellResponseCorpus.open(secondCorpus, sourceStore: secondStore,
                                                            implementation: Self.implementation)
        let firstReceipt = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(firstReader.receipt))
        let secondReceipt = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(secondReader.receipt))
        let forward = try VivoCellResponseCompositeCorpusReader(readers: [firstReader, secondReader])
        let reverse = try VivoCellResponseCompositeCorpusReader(readers: [secondReader, firstReader])
        #expect(forward.identity == reverse.identity)
        #expect(forward.fingerprint == reverse.fingerprint)
        #expect(forward.sourceCount == 2)
        let firstRows = try forward.rows(forSampleID: "target-training", corpus: firstReceipt)
        let secondRows = try forward.rows(forSampleID: "target-training", corpus: secondReceipt)
        #expect(!firstRows.isEmpty && !secondRows.isEmpty)
        #expect(Set(firstRows).isDisjoint(with: Set(secondRows)))
        #expect(try forward.sourceCorpus(forGlobalRow: firstRows[0]) == firstReceipt)
        #expect(try forward.sourceCorpus(forGlobalRow: secondRows[0]) == secondReceipt)
        let firstLocalRows = try firstReader.rows(forSampleID: "target-training")
        let secondLocalRows = try secondReader.rows(forSampleID: "target-training")
        let firstCompositeValues = try forward.normalizedRow(firstRows[0])
        let secondCompositeValues = try forward.normalizedRow(secondRows[0])
        let firstSourceValues = try firstReader.normalizedRow(firstLocalRows[0])
        let secondSourceValues = try secondReader.normalizedRow(secondLocalRows[0])
        #expect(firstCompositeValues == firstSourceValues)
        #expect(secondCompositeValues == secondSourceValues)
        let training = try forward.examples(in: .training)
        let expectedStrata = Set(try firstReader.examples(in: .training).map(\.stratumIndex)).count +
            Set(try secondReader.examples(in: .training).map(\.stratumIndex)).count
        #expect(Set(training.map(\.stratumIndex)).count == expectedStrata)
        #expect(Set(training.map(\.targetRow)).count == training.count)
        let requirements = VivoCellResponseCohortRequirements(minimumTrainingBiologicalUnitsPerTarget: 1,
                                                                minimumValidationBiologicalUnitsPerTarget: 1,
                                                                minimumTestBiologicalUnitsPerTarget: 1,
                                                                minimumBiologicalUnits: 3)
        let admitted = try VivoCellResponseCohort.admit(readers: [secondReader, firstReader], requirements: requirements)
        #expect(admitted.sources.count == 2)
        try VivoCellResponseCohort.verify(admitted, readers: [firstReader, secondReader])
        #expect(throws: (any Error).self) { try VivoCellResponseCohort.verify(admitted, readers: [firstReader]) }
        #expect(throws: (any Error).self) { try VivoCellResponseCompositeCorpusReader(readers: [firstReader, firstReader]) }
        let mismatchedCorpus = root.appendingPathComponent("mismatched-corpus")
        _ = try VivoCellResponseCorpus.prepare(sourceStore: secondStore,
                                                plan: Self.corpusPlan(featureIDs: ["gene-a", "gene-b"], replicatePrefix: "second-"),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: mismatchedCorpus)
        let mismatchedReader = try VivoCellResponseCorpus.open(mismatchedCorpus, sourceStore: secondStore,
                                                                implementation: Self.implementation)
        #expect(throws: (any Error).self) {
            try VivoCellResponseCompositeCorpusReader(readers: [firstReader, mismatchedReader])
        }
    }

    @Test("multi-source learner masks structural absences and requires a qualified prediction source")
    func multiSourceLifecycleMasksPanelsAndReplays() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(
            "numivivo-cell-response-multisource-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let first = Self.dataset(), second = Self.secondPanelDataset()
        let firstH5AD = root.appendingPathComponent("first.h5ad")
        let secondH5AD = root.appendingPathComponent("second-panel.h5ad")
        let firstStore = root.appendingPathComponent("first-store")
        let secondStore = root.appendingPathComponent("second-store")
        let firstCorpus = root.appendingPathComponent("first-corpus")
        let secondCorpus = root.appendingPathComponent("second-corpus")
        let model = root.appendingPathComponent("model")
        let resumed = root.appendingPathComponent("resumed")
        let uninterrupted = root.appendingPathComponent("uninterrupted")
        let evaluationDirectory = root.appendingPathComponent("evaluation")
        let firstPrediction = root.appendingPathComponent("first-prediction")
        let secondPrediction = root.appendingPathComponent("second-prediction")
        let resumedFirstPrediction = root.appendingPathComponent("resumed-first-prediction")
        let resumedSecondPrediction = root.appendingPathComponent("resumed-second-prediction")
        let uninterruptedFirstPrediction = root.appendingPathComponent("uninterrupted-first-prediction")
        let uninterruptedSecondPrediction = root.appendingPathComponent("uninterrupted-second-prediction")
        try VivoSingleCellH5AD.write(first, to: firstH5AD)
        try VivoSingleCellH5AD.write(second, to: secondH5AD)
        _ = try VivoH5ADCountStore.publish(source: firstH5AD, plan: Self.importPlan(samples: first.samples),
                                            implementation: Self.implementation, to: firstStore)
        _ = try VivoH5ADCountStore.publish(source: secondH5AD, plan: Self.importPlan(samples: second.samples),
                                            implementation: Self.implementation, to: secondStore)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: firstStore, plan: Self.corpusPlan(),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: firstCorpus)
        _ = try VivoCellResponseCorpus.prepare(sourceStore: secondStore,
                                                plan: Self.corpusPlan(replicatePrefix: "second-"),
                                                descriptorSourceBytes: Self.descriptorSourceBytes,
                                                implementation: Self.implementation, to: secondCorpus)
        let firstReader = try VivoCellResponseCorpus.open(firstCorpus, sourceStore: firstStore,
                                                           implementation: Self.implementation)
        let secondReader = try VivoCellResponseCorpus.open(secondCorpus, sourceStore: secondStore,
                                                            implementation: Self.implementation)
        #expect(firstReader.featureMask == [1, 1, 0])
        #expect(secondReader.featureMask == [1, 0, 0])
        let firstReceipt = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(firstReader.receipt))
        let secondReceipt = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(secondReader.receipt))
        let composite = try VivoCellResponseCompositeCorpusReader(readers: [firstReader, secondReader])
        let reversed = try VivoCellResponseCompositeCorpusReader(readers: [secondReader, firstReader])
        #expect(composite.identity == reversed.identity)
        #expect(composite.fingerprint == reversed.fingerprint)
        #expect(composite.featureMask == [1, 1, 0])
        let firstTargetRow = try composite.rows(forSampleID: "target-validation", corpus: firstReceipt)[0]
        let secondTargetRow = try composite.rows(forSampleID: "target-validation", corpus: secondReceipt)[0]
        #expect(try composite.featureMask(forGlobalRow: firstTargetRow) == [1, 1, 0])
        #expect(try composite.featureMask(forGlobalRow: secondTargetRow) == [1, 0, 0])

        let plan = VivoCellResponseTrainingPlan(
            id: "synthetic-multi-source", architecture: .init(hiddenWidth: 16, contextCells: 2,
                                                                  maximumParameters: 100_000),
            steps: 4, batchSize: 2, learningRate: 0.01, weightDecay: 0,
            seed: 37, validationEvery: 4, validationExamples: 4)
        let insufficientValidationBudget = VivoCellResponseTrainingPlan(
            id: "synthetic-multi-source-insufficient-validation", architecture: plan.architecture,
            steps: plan.steps, batchSize: plan.batchSize, learningRate: plan.learningRate,
            weightDecay: plan.weightDecay, seed: plan.seed,
            validationEvery: plan.validationEvery, validationExamples: 1)
        #expect(throws: (any Error).self) {
            try insufficientValidationBudget.validate(for: composite)
        }
        #expect(try VivoCellResponseLearning.trainingControlCacheMatchesUncachedBatchForTesting(
            corpus: composite, plan: plan, step: 1))
        let requirements = VivoCellResponseCohortRequirements(
            minimumTrainingBiologicalUnitsPerTarget: 1,
            minimumValidationBiologicalUnitsPerTarget: 1,
            minimumTestBiologicalUnitsPerTarget: 1,
            minimumBiologicalUnits: 3)
        let cohort = try VivoCellResponseCohort.admit(readers: [secondReader, firstReader], requirements: requirements)
        try VivoCellResponseCohort.verify(cohort, readers: composite.sourceReaders)
        let receipt = try VivoCellResponseLearning.train(corpus: composite, plan: plan, cohort: cohort,
                                                          implementation: Self.implementation, to: model)
        #expect(receipt.format == VivoCellResponseLearning.modelFormat)
        #expect(receipt.corpus == composite.fingerprint)
        _ = try VivoCellResponseLearning.verifyModel(model, implementation: Self.implementation)
        let evaluation = try VivoCellResponseLearning.evaluate(model: model, corpus: composite, partition: .validation,
                                                                maximumExamples: 4, implementation: Self.implementation)
        #expect(evaluation.examples == 4)
        #expect(evaluation.observedFeatures == 6)
        #expect(evaluation.sourceMetrics.map(\.corpus) == composite.identity.sources.map(\.corpus))
        #expect(evaluation.sourceMetrics.allSatisfy { $0.examples == 2 })
        #expect(evaluation.sourceMetrics.reduce(0) { $0 + $1.observedFeatures } == 6)
        let selectedFeatureCount = try evaluation.selection.reduce(into: 0) { total, selected in
            total += try composite.featureMask(forGlobalRow: selected.targetRow).reduce(0) { $0 + ($1 > 0 ? 1 : 0) }
        }
        #expect(selectedFeatureCount == evaluation.observedFeatures)
        let sourceBalancedEvaluation = try VivoCellResponseLearning.evaluate(
            model: model, corpus: composite, partition: .validation,
            maximumExamples: 2, implementation: Self.implementation)
        #expect(sourceBalancedEvaluation.examples == 2)
        #expect(sourceBalancedEvaluation.sourceMetrics.count == 2)
        #expect(sourceBalancedEvaluation.sourceMetrics.allSatisfy { $0.examples == 1 })
        #expect(sourceBalancedEvaluation.sourceMetrics.map(\.observedFeatures).sorted() == [1, 2])
        #expect(sourceBalancedEvaluation.observedFeatures == 3)
        let sourceMacroRMSE = sourceBalancedEvaluation.sourceMetrics.reduce(0) { $0 + $1.rmse } /
            Double(sourceBalancedEvaluation.sourceMetrics.count)
        #expect(abs(sourceBalancedEvaluation.rmse - sourceMacroRMSE) <= 1e-10)
        #expect(throws: (any Error).self) {
            _ = try VivoCellResponseLearning.evaluate(model: model, corpus: composite, partition: .validation,
                                                       maximumExamples: 1, implementation: Self.implementation)
        }
        _ = try VivoCellResponseLearning.evaluate(model: model, corpus: composite, partition: .validation,
                                                   maximumExamples: 4, implementation: Self.implementation,
                                                   to: evaluationDirectory)
        #expect(try VivoCellResponseLearning.verifyEvaluation(evaluationDirectory, model: model, corpus: reversed,
                                                               implementation: Self.implementation) == evaluation)

        let bareContext = VivoCellResponsePredictionPlan(
            id: "ambiguous-context", target: .init(id: "novel-target", descriptors: [0.1, 0.2]),
            contextSampleIDs: ["control-validation"])
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.predict(model: model, corpus: composite, plan: bareContext,
                                                  implementation: Self.implementation,
                                                  to: root.appendingPathComponent("ambiguous-prediction"))
        }
        let firstContext = VivoCellResponsePredictionPlan(
            id: "first-context", target: bareContext.target, contextSampleIDs: bareContext.contextSampleIDs,
            contextCorpus: firstReceipt)
        let secondContext = VivoCellResponsePredictionPlan(
            id: "second-context", target: bareContext.target, contextSampleIDs: bareContext.contextSampleIDs,
            contextCorpus: secondReceipt)
        #expect(try VivoCanonicalJSON.decode(VivoCellResponsePredictionPlan.self,
                                              from: VivoCanonicalJSON.encode(firstContext)) == firstContext)
        _ = try VivoCellResponseLearning.predict(model: model, corpus: composite, plan: firstContext,
                                                  implementation: Self.implementation, to: firstPrediction)
        _ = try VivoCellResponseLearning.predict(model: model, corpus: composite, plan: secondContext,
                                                  implementation: Self.implementation, to: secondPrediction)
        let firstResult = try VivoCellResponseLearning.verifyPrediction(firstPrediction, model: model, corpus: reversed,
                                                                         implementation: Self.implementation)
        let secondResult = try VivoCellResponseLearning.verifyPrediction(secondPrediction, model: model, corpus: reversed,
                                                                          implementation: Self.implementation)
        #expect(firstResult.featureMask == [1, 1, 0])
        #expect(secondResult.featureMask == [1, 0, 0])
        let secondPredictionSources = try secondResult.contextSourceRows.map {
            try reversed.sourceCorpus(forGlobalRow: $0)
        }
        #expect(secondPredictionSources.allSatisfy { $0 == secondReceipt })
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.resume(model: model, corpus: firstReader,
                                                 plan: .init(id: "missing-source", additionalSteps: 1),
                                                 implementation: Self.implementation,
                                                 to: root.appendingPathComponent("missing-source-resume"))
        }
        _ = try VivoCellResponseLearning.resume(model: model, corpus: reversed,
                                                 plan: .init(id: "multi-source-resume", additionalSteps: 1),
                                                 implementation: Self.implementation, to: resumed)
        _ = try VivoCellResponseLearning.verifyModel(resumed, implementation: Self.implementation)
        let uninterruptedPlan = VivoCellResponseTrainingPlan(
            id: "synthetic-multi-source-uninterrupted", architecture: plan.architecture,
            steps: 5, batchSize: plan.batchSize, learningRate: plan.learningRate,
            weightDecay: plan.weightDecay, seed: plan.seed,
            validationEvery: plan.validationEvery, validationExamples: plan.validationExamples)
        _ = try VivoCellResponseLearning.train(corpus: reversed, plan: uninterruptedPlan, cohort: cohort,
                                                implementation: Self.implementation, to: uninterrupted)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: composite, plan: firstContext,
                                                  implementation: Self.implementation, to: resumedFirstPrediction)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: composite, plan: secondContext,
                                                  implementation: Self.implementation, to: resumedSecondPrediction)
        _ = try VivoCellResponseLearning.predict(model: uninterrupted, corpus: reversed, plan: firstContext,
                                                  implementation: Self.implementation, to: uninterruptedFirstPrediction)
        _ = try VivoCellResponseLearning.predict(model: uninterrupted, corpus: reversed, plan: secondContext,
                                                  implementation: Self.implementation, to: uninterruptedSecondPrediction)
        let resumedFirst = try VivoCellResponseLearning.verifyPrediction(resumedFirstPrediction, implementation: Self.implementation)
        let resumedSecond = try VivoCellResponseLearning.verifyPrediction(resumedSecondPrediction, implementation: Self.implementation)
        let uninterruptedFirst = try VivoCellResponseLearning.verifyPrediction(uninterruptedFirstPrediction,
                                                                                implementation: Self.implementation)
        let uninterruptedSecond = try VivoCellResponseLearning.verifyPrediction(uninterruptedSecondPrediction,
                                                                                 implementation: Self.implementation)
        #expect(zip(resumedFirst.meanLogCPM, uninterruptedFirst.meanLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(resumedFirst.varianceLogCPM, uninterruptedFirst.varianceLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(resumedSecond.meanLogCPM, uninterruptedSecond.meanLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
        #expect(zip(resumedSecond.varianceLogCPM, uninterruptedSecond.varianceLogCPM).allSatisfy { abs($0.0 - $0.1) <= 1e-6 })
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
            .init(schemaVersion: 4, format: VivoCellResponseLearning.evaluationFormat,
                  partition: .validation, maximumExamples: 1,
                  samplerVersion: VivoCellResponseLearning.samplerVersion, seed: 0,
                  examples: 1, observedFeatures: 1, negativeLogLikelihood: 0,
                  rmse: rmse, matchedControlRMSE: baseline,
                  selection: [.init(targetRow: 0, contextRows: [1])],
                  sourceMetrics: [.init(corpus: Self.implementation, examples: 1, observedFeatures: 1,
                                        negativeLogLikelihood: 0, rmse: rmse, matchedControlRMSE: baseline)])
        }

        try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 0.9, baseline: 1))
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 1, baseline: 1))
        }
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireMatchedControlImprovement(evaluation(rmse: 1.1, baseline: 1))
        }
        let secondSource: VivoFingerprint = try! .init(bytes: Array(repeating: 8, count: 32))
        let sourceMetrics = [
            VivoCellResponseSourceEvaluation(corpus: Self.implementation, examples: 1, observedFeatures: 1,
                                              negativeLogLikelihood: 0, rmse: 0.5, matchedControlRMSE: 1),
            VivoCellResponseSourceEvaluation(corpus: secondSource, examples: 1, observedFeatures: 1,
                                              negativeLogLikelihood: 0, rmse: 1.3, matchedControlRMSE: 1)
        ]
        let hiddenSourceRegression = VivoCellResponseEvaluation(
            schemaVersion: 4, format: VivoCellResponseLearning.evaluationFormat,
            partition: .validation, maximumExamples: 2,
            samplerVersion: VivoCellResponseLearning.samplerVersion, seed: 0,
            examples: 2, observedFeatures: 2, negativeLogLikelihood: 0,
            rmse: 0.9, matchedControlRMSE: 1,
            selection: [.init(targetRow: 0, contextRows: [1]), .init(targetRow: 2, contextRows: [3])],
            sourceMetrics: sourceMetrics)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.requireMatchedControlImprovement(hiddenSourceRegression)
        }
    }
}
