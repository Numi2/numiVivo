import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileExpressionTests {
    let helpers = SingleCellFileAxisTests()
    let sourceImplementation = SingleCellFileAxisTests.implementation
    var analysisImplementation: VivoFingerprint { get throws { try VivoCanonicalJSON.fingerprint(Data("different-analysis-owner".utf8)) } }
    func counts(_ root: URL, _ data: VivoSingleCellDataset) throws -> URL {
        let axis = try helpers.axis(root, data)
        var bytes = Data()
        for row in data.cells.indices {
            for k in data.matrix.rowOffsets[row]..<data.matrix.rowOffsets[row + 1] {
                bytes.vivoAppendLE(UInt32(row)); bytes.vivoAppendLE(UInt32(data.matrix.featureIndices[k])); bytes.vivoAppendLE(data.matrix.counts[k])
            }
        }
        return try helpers.publish(root, axis, bytes)
    }
    func plan(_ directory: URL, _ contrast: VivoOmicsExpressionContrast) throws -> VivoFileExpressionPlan {
        .init(sourceReceipt: try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent("receipt.json"), maximumBytes: 65_536),
              sourceImplementation: sourceImplementation, contrast: contrast)
    }
    func checkStreamedBytes(_ report: VivoOmicsExpressionResult) throws {
        let expected = try VivoCanonicalJSON.encode(report)
        var streamed = Data()
        let digest = try VivoExpressionReportJSON.fingerprint(report, maximumBytes: expected.count) { streamed.append($0) }
        #expect(streamed == expected)
        #expect(digest == (try VivoCanonicalJSON.fingerprint(expected)))
        #expect(try VivoExpressionReportJSON.fingerprint(report, maximumBytes: expected.count) == digest)
    }
    func checkAllStatistics(_ actual: VivoOmicsExpressionResult, _ expected: VivoOmicsExpressionResult) throws {
        try checkStreamedBytes(actual); try checkStreamedBytes(expected)
        var a = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(actual)) as? [String: Any])
        let b = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(expected)) as? [String: Any])
        var design = try #require(a["design"] as? [String: Any])
        var observations = try #require(design["observations"] as? [[String: Any]])
        for i in actual.design.observations.indices {
            let row = actual.design.observations[i]
            #expect(row.sourceCellIndices == nil)
            #expect(row.sourceCellCount == expected.design.observations[i].sourceCellCount)
            if case .file(let reference) = row.membership { #expect(reference.groupIndex == actual.design.sourcePseudobulkIndices[i]) }
            else { Issue.record("File observation lost its explicit membership reference") }
            observations[i].removeValue(forKey: "fileMembership")
            observations[i]["sourceCellIndices"] = try #require(expected.design.observations[i].sourceCellIndices)
        }
        design["observations"] = observations; a["design"] = design
        #expect(NSDictionary(dictionary: a).isEqual(to: b))
    }
    @Test func legacyObservationEncodingIsExactAndFileEncodingIsExplicit() throws {
        let bulk = try VivoSingleCellAnalysis.pseudobulk(VivoSingleCellExamples.pairedCounts())
        for group in bulk.groups {
            let row = VivoOmicsDesignObservation(group), bytes = try VivoCanonicalJSON.encode(group)
            #expect(try VivoCanonicalJSON.encode(row) == bytes)
            #expect(try VivoCanonicalJSON.decode(VivoOmicsDesignObservation.self, from: bytes) == row)
        }
        let original = try #require(bulk.groups.first)
        let group = VivoFilePseudobulkGroup(biologicalReplicateID: original.biologicalReplicateID, donorID: original.donorID,
            condition: original.condition, organism: original.organism, cellGroup: original.cellGroup,
            sampleIDs: original.sampleIDs, batchIDs: original.batchIDs, sourceCellCount: 3)
        let row = VivoOmicsDesignObservation(group, receipt: sourceImplementation, index: 0), bytes = try VivoCanonicalJSON.encode(row)
        #expect(row.sourceCellIndices == nil && row.sourceCellCount == 3)
        #expect(try VivoCanonicalJSON.decode(VivoOmicsDesignObservation.self, from: bytes) == row)
        var object = try #require(JSONSerialization.jsonObject(with: bytes) as? [String: Any])
        object["sourceCellIndices"] = []
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoOmicsDesignObservation.self, from: JSONSerialization.data(withJSONObject: object)) }
        object.removeValue(forKey: "sourceCellIndices"); object.removeValue(forKey: "fileMembership")
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoOmicsDesignObservation.self, from: JSONSerialization.data(withJSONObject: object)) }
        object["fileMembership"] = ["countBundleReceipt": ["bytes": Array(sourceImplementation.bytes)], "groupIndex": 0, "sourceCellCount": 3, "discardedField": true]
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoOmicsDesignObservation.self, from: JSONSerialization.data(withJSONObject: object)) }
    }
    @Test func allWaldLikelihoodRatioAndQLResultsMatchResidentOwner() throws {
        try helpers.workspace { root in
            let data = try SingleCellNBCohortTests.fixture(), directory = try counts(root, data)
            let methods: [VivoOmicsNBTestMethod?] = [nil, .likelihoodRatio, .quasiLikelihoodAdjusted]
            for method in methods {
                var contrast = VivoOmicsExpressionContrast(id: "paired", controlCondition: "ctrl", treatmentCondition: "stim", design: .pairedDonors)
                contrast.model = .negativeBinomial; contrast.minimumCellsPerPseudobulk = 1
                var options = VivoOmicsNBCohortOptions(); options.trend = .mean; options.testMethod = method
                contrast.negativeBinomialOptions = options
                let expected = try VivoPseudobulkDifferentialExpression.run(data, contrast: contrast)
                let actual = try VivoFileExpression.run(source: directory, plan: plan(directory, contrast))
                try checkAllStatistics(actual, expected)
                #expect(actual.testedFeatures > 80)
            }
        }
    }
    @Test func emptyCellsCountForFilteringAndDistinctAnalysisOwnerReplays() throws {
        try helpers.workspace { root in
            let data = try VivoSingleCellExamples.pairedCounts(), directory = try counts(root, data)
            var contrast = VivoSingleCellExamples.pairedPlan().contrasts[0]; contrast.minimumCellsPerPseudobulk = 3
            let declared = try plan(directory, contrast), output = root.appendingPathComponent("analysis")
            let receipt = try VivoFileExpression.publish(source: directory, plan: declared, implementation: analysisImplementation, to: output)
            #expect(receipt.sourceImplementation == sourceImplementation && receipt.implementation != sourceImplementation)
            let replay = try VivoFileExpression.verify(output, implementation: analysisImplementation)
            try checkAllStatistics(replay, VivoPseudobulkDifferentialExpression.run(data, contrast: contrast))
            #expect(throws: (any Error).self) { try VivoFileExpression.verify(output, implementation: sourceImplementation) }
            contrast.minimumCellsPerPseudobulk = 4
            #expect(throws: (any Error).self) { try VivoFileExpression.run(source: directory, plan: plan(directory, contrast)) }
            #expect(throws: (any Error).self) { try VivoPseudobulkDifferentialExpression.run(data, contrast: contrast) }
        }
    }
    @Test func rehashedAlteredStatisticalReportDoesNotReconstruct() throws {
        try helpers.workspace { root in
            let data = try VivoSingleCellExamples.pairedCounts(), directory = try counts(root, data), output = root.appendingPathComponent("analysis")
            let receipt = try VivoFileExpression.publish(source: directory, plan: plan(directory, VivoSingleCellExamples.pairedPlan().contrasts[0]), implementation: analysisImplementation, to: output)
            var object = try #require(JSONSerialization.jsonObject(with: Data(contentsOf: output.appendingPathComponent("report.json"))) as? [String: Any])
            var features = try #require(object["features"] as? [[String: Any]])
            features[0]["log2FoldChange"] = 123.0; object["features"] = features
            let raw = try JSONSerialization.data(withJSONObject: object, options: [.sortedKeys]); try raw.write(to: output.appendingPathComponent("report.json"))
            let changed = try VivoFileExpressionReceipt(schemaVersion: receipt.schemaVersion, method: receipt.method, sourceReceipt: receipt.sourceReceipt,
                sourceImplementation: receipt.sourceImplementation, plan: receipt.plan, report: VivoCanonicalJSON.fingerprint(raw), implementation: receipt.implementation)
            try VivoCanonicalJSON.encode(changed).write(to: output.appendingPathComponent("receipt.json"))
            #expect(throws: (any Error).self) { try VivoFileExpression.verify(output, implementation: analysisImplementation) }
        }
    }
    @Test func snapshotCopyPreservesOriginalInputAfterSourceMutation() throws {
        try helpers.workspace { root in
            let data = try VivoSingleCellExamples.pairedCounts(), directory = try counts(root, data)
            let declared = try plan(directory, VivoSingleCellExamples.pairedPlan().contrasts[0])
            let source = try VivoFileCountSnapshot.open(directory, implementation: sourceImplementation)
            let expected = try VivoFileExpression.evaluate(source, plan: declared)
            try Data([0]).write(to: directory.appendingPathComponent("quality.bin"))
            let copy = root.appendingPathComponent("copy"); try source.copy(to: copy)
            #expect(try VivoFileExpression.run(source: copy, plan: declared) == expected)
            #expect(throws: (any Error).self) { try VivoFileExpression.run(source: directory, plan: declared) }
        }
    }
    @Test func sourceIdentityMismatchFailsBeforePublication() throws {
        try helpers.workspace { root in
            let directory = try counts(root, VivoSingleCellExamples.pairedCounts()), output = root.appendingPathComponent("wrong")
            let contrast = VivoSingleCellExamples.pairedPlan().contrasts[0], correct = try plan(directory, contrast)
            for incorrect in [VivoFileExpressionPlan(sourceReceipt: sourceImplementation, sourceImplementation: sourceImplementation, contrast: contrast),
                              VivoFileExpressionPlan(sourceReceipt: correct.sourceReceipt, sourceImplementation: try analysisImplementation, contrast: contrast)] {
                #expect(throws: (any Error).self) { try VivoFileExpression.publish(source: directory, plan: incorrect, implementation: analysisImplementation, to: output) }
                #expect(!FileManager.default.fileExists(atPath: output.path))
            }
        }
    }
    @Test func streamedReportsPreserveOptionalPriorsAndFailedQLDiagnostics() throws {
        let data = try SingleCellNBCohortTests.fixture()
        var contrast = VivoOmicsExpressionContrast(id: "paired", controlCondition: "ctrl", treatmentCondition: "stim", design: .pairedDonors)
        contrast.model = .negativeBinomial; contrast.minimumCellsPerPseudobulk = 1
        for mode in 0..<3 {
            var options = VivoOmicsNBCohortOptions(); options.trend = .mean
            if mode == 0 { options.effectPriorStandardDeviationLog2 = 1 }
            if mode == 1 { options.effectPriorEstimation = .weightedUpperQuantile }
            if mode == 2 { options.testMethod = .quasiLikelihoodAdjusted }
            contrast.negativeBinomialOptions = options
            var result = try VivoPseudobulkDifferentialExpression.run(data, contrast: contrast)
            try checkStreamedBytes(result)
            if mode == 2 {
                var nb = try #require(result.negativeBinomial)
                let ql = try #require(nb.quasiLikelihood)
                let failure = VivoOmicsNBQLFailure(featureIndex: nil, stage: "é/β", message: "quoted \"failure\"\n")
                nb.effectPriorEstimationError = "retained diagnostic"
                nb.quasiLikelihood = .init(featureIndices: ql.featureIndices, averageQLDispersion: ql.averageQLDispersion,
                    abundanceFits: ql.abundanceFits, moderation: ql.moderation, inference: ql.inference, failures: [failure],
                    failedUpstreamFit: .init(completed: false, abundanceFits: [nil], globalFit: .init(completed: false,
                        averageQuasiDispersion: 0.5, initialFits: [nil], updates: [], refittedFits: [nil], adjustedResiduals: [nil], failures: [failure]), failures: [failure]))
                result.negativeBinomial = nb
                try checkStreamedBytes(result)
            }
        }
    }
    @Test func streamedReportRejectsOverflowSinkFailureAndCancellation() async throws {
        let result = try VivoPseudobulkDifferentialExpression.run(VivoSingleCellExamples.pairedCounts(), contrast: VivoSingleCellExamples.pairedPlan().contrasts[0])
        let expected = try VivoCanonicalJSON.encode(result)
        #expect(throws: (any Error).self) { try VivoExpressionReportJSON.fingerprint(result, maximumBytes: expected.count - 1) }
        enum SinkFailure: Error { case stopped }
        #expect(throws: SinkFailure.self) {
            try VivoExpressionReportJSON.fingerprint(result, maximumBytes: expected.count) { _ in throw SinkFailure.stopped }
        }
        let cancelled = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try VivoExpressionReportJSON.fingerprint(result, maximumBytes: expected.count)
        }
        do { _ = try await cancelled.value; Issue.record("Cancelled encoder succeeded") }
        catch is CancellationError { }
    }

}
