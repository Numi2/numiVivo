import Foundation
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
        let sampleIDs = ["control-training", "target-training", "control-validation", "target-validation", "empty-metadata-sample"]
        let sourceSamples = [
            sample(sampleIDs[0], condition: sampleIDs[0], batchID: "context-training", replicateID: "pair-training"),
            sample(sampleIDs[1], condition: "target-a", batchID: "context-training", replicateID: "pair-training"),
            sample(sampleIDs[2], condition: sampleIDs[2], batchID: "context-validation", replicateID: "pair-validation"),
            sample(sampleIDs[3], condition: "target-b", batchID: "context-validation", replicateID: "pair-validation"),
            sample(sampleIDs[4])
        ]
        let cells = [
            VivoOmicsCell(barcode: "ct0", sampleID: sampleIDs[0]), VivoOmicsCell(barcode: "ct1", sampleID: sampleIDs[0]),
            VivoOmicsCell(barcode: "tt0", sampleID: sampleIDs[1]), VivoOmicsCell(barcode: "tt1", sampleID: sampleIDs[1]),
            VivoOmicsCell(barcode: "cv0", sampleID: sampleIDs[2]), VivoOmicsCell(barcode: "cv1", sampleID: sampleIDs[2]),
            VivoOmicsCell(barcode: "tv0", sampleID: sampleIDs[3]), VivoOmicsCell(barcode: "tv1", sampleID: sampleIDs[3])
        ]
        // The first row deliberately omits gene-b. It is a measured zero,
        // while gene-absent is structurally absent from the source axis.
        let rows: [[(Int, UInt64)]] = [
            [(0, 8)], [(0, 7), (1, 1)],
            [(0, 1), (1, 8)], [(0, 2), (1, 7)],
            [(0, 8), (1, 2)], [(0, 7), (1, 2)],
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

    private static func importPlan() -> VivoH5ADImportPlan {
        .init(id: "synthetic-import", evidence: .synthetic,
              sourceDescription: "Synthetic software-only response fixture", countUnit: .umiCount,
              matrixPath: "X", samples: dataset().samples, sampleColumn: "sample", barcodeColumn: "barcode",
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

    private static func corpusPlan() -> VivoCellResponseCorpusPlan {
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
                     featureAxis: .init(featureIDs: ["gene-a", "gene-b", "gene-absent"]), targets: targets,
                     descriptorSource: descriptorSourceFingerprint,
                     assignments: [
                        base("control-training", .control, .training, "context-training", "pair-training", nil),
                        base("target-training", .perturbed, .training, "context-training", "pair-training", "target-a"),
                        base("control-validation", .control, .validation, "context-validation", "pair-validation", nil),
                        base("target-validation", .perturbed, .validation, "context-validation", "pair-validation", "target-b")
                     ], heldOutTargetIDs: ["target-b"])
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
        #expect(try reader.trainingTargetBindings().map(\.id) == ["target-a"])
        let first = try reader.normalizedRow(0)
        #expect(first[1] == 0)
        #expect(first[2] == -1)

        let plan = VivoCellResponseTrainingPlan(id: "synthetic-train",
                                                 architecture: .init(hiddenWidth: 16, contextCells: 2, maximumParameters: 100_000),
                                                 steps: 12, batchSize: 1, learningRate: 0.01, weightDecay: 0,
                                                 seed: 17, validationEvery: 4, validationExamples: 2)
        let receipt = try VivoCellResponseLearning.train(corpus: reader, plan: plan, implementation: Self.implementation, to: model)
        #expect(receipt.format == VivoCellResponseLearning.modelFormat)
        _ = try VivoCellResponseLearning.verifyModel(model, implementation: Self.implementation)
        let evaluation = try VivoCellResponseLearning.evaluate(model: model, corpus: reader, partition: .validation,
                                                                maximumExamples: 2, implementation: Self.implementation)
        // Validation has one target×context stratum with two cells. Macro
        // sampling scores one deterministic cell rather than overweighting it.
        #expect(evaluation.examples == 1)
        #expect(evaluation.observedFeatures == 2)
        #expect(evaluation.negativeLogLikelihood.isFinite)
        #expect(evaluation.rmse.isFinite)
        #expect(evaluation.matchedControlRMSE.isFinite)
        #expect(evaluation.selection.count == evaluation.examples)
        let evaluationReceipt = try VivoCellResponseLearning.evaluate(model: model, corpus: reader, partition: .validation,
                                                                       maximumExamples: 2, implementation: Self.implementation,
                                                                       to: persistedEvaluation)
        #expect(evaluationReceipt.format == VivoCellResponseLearning.evaluationFormat)
        #expect(try VivoCellResponseLearning.verifyEvaluation(persistedEvaluation, model: model, corpus: reader,
                                                              implementation: Self.implementation) == evaluation)

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

        let heldOutKnownTarget = VivoCellResponsePredictionPlan(
            id: "heldout-known-target", target: .init(id: "target-b", descriptors: [0.6, 0.4]),
            contextSampleIDs: ["control-validation"], useTrainedTargetEmbedding: true)
        #expect(throws: (any Error).self) {
            try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: heldOutKnownTarget,
                                                  implementation: Self.implementation,
                                                  to: root.appendingPathComponent("heldout-known-target"))
        }
        let heldOutDescriptorOnly = VivoCellResponsePredictionPlan(
            id: "heldout-descriptor-only", target: .init(id: "target-b", descriptors: [0.6, 0.4]),
            contextSampleIDs: ["control-validation"], useTrainedTargetEmbedding: false)
        _ = try VivoCellResponseLearning.predict(model: resumed, corpus: reader, plan: heldOutDescriptorOnly,
                                                  implementation: Self.implementation,
                                                  to: root.appendingPathComponent("heldout-descriptor-only"))

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
                                                implementation: Self.implementation, to: straightThrough)
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
