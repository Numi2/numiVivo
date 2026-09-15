import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct QuantitativeAssayTests {
    private static func fixture() -> VivoQuantitativeAssayDataset {
        let sample = VivoOmicsSample(id: "sample", biologicalReplicateID: "rep",
                                     condition: "control", batchID: "batch", organism: "NCBITaxon:9606")
        let observations: [VivoAssayObservation] = [
            .init(identity: .init(barcode: "b0", sampleID: "sample"), kind: .cell, position: nil),
            .init(identity: .init(barcode: "b1", sampleID: "sample"), kind: .cell, position: nil),
            .init(identity: .init(barcode: "b2", sampleID: "sample"), kind: .spot,
                  position: .init(frameID: "slide", coordinates: [4, 8]))
        ]
        let protein = VivoQuantitativeAssaySpace(
            id: "protein", kind: .proteomics, featureNamespace: "uniprot",
            unit: "log2-intensity", sourceDescription: "Measured protein panel",
            features: [.init(id: "P1", name: "protein-1", interval: nil)],
            observationIndices: [0, 2],
            matrix: .init(observationCount: 2, featureCount: 1,
                          rowOffsets: [0, 1, 2], featureIndices: [0, 0], values: [0, 1.25]))
        let metabolite = VivoQuantitativeAssaySpace(
            id: "metabolite", kind: .metabolomics, featureNamespace: "hmdb",
            unit: "log10-concentration", sourceDescription: "Measured metabolite panel",
            features: [.init(id: "HMDB1", name: "metabolite-1", interval: nil)],
            observationIndices: [1],
            matrix: .init(observationCount: 1, featureCount: 1,
                          rowOffsets: [0, 1], featureIndices: [0], values: [-2.5]))
        return .init(id: "quantitative", evidence: .measured,
                     sourceDescription: "Bounded quantitative-assay fixture",
                     samples: [sample], observations: observations,
                     spatialFrames: [.init(id: "slide", unit: .pixel, axes: ["x", "y"],
                                           sourceDescription: "Measured image frame")],
                     assays: [protein, metabolite])
    }

    @Test func preservesMissingAndMeasuredZeroAcrossJSONRoundTrip() throws {
        let dataset = Self.fixture()
        try dataset.validate()
        #expect(try dataset.value(assayID: "protein", observationIndex: 0, featureID: "P1") == 0)
        #expect(try dataset.value(assayID: "protein", observationIndex: 1, featureID: "P1") == nil)
        #expect(try dataset.value(assayID: "protein", observationIndex: 2, featureID: "P1") == 1.25)
        #expect(try dataset.value(assayID: "metabolite", observationIndex: 1, featureID: "HMDB1") == -2.5)
        #expect(try VivoCanonicalJSON.decode(VivoQuantitativeAssayDataset.self,
                                              from: VivoCanonicalJSON.encode(dataset)) == dataset)
    }

    @Test func descriptiveSummaryBindsCoverageMomentsAndMissingness() throws {
        let dataset = Self.fixture()
        let summary = try dataset.descriptiveSummary()
        #expect(summary.method == VivoQuantitativeAssaySummaries.method)
        let expectedFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(dataset))
        #expect(summary.datasetFingerprint == expectedFingerprint)
        #expect(summary.observationCount == 3)
        #expect(summary.assayCount == 2)
        let protein = try #require(summary.coverage.first { $0.assayID == "protein" })
        #expect(protein.assayRowCount == 2)
        #expect(protein.featureCount == 1)
        #expect(protein.measuredValueCount == 2)
        #expect(protein.missingValueCount == 1)
        #expect(protein.measuredZeroCount == 1)
        let proteinFeature = try #require(summary.features.first { $0.featureID == "P1" })
        #expect(proteinFeature.measuredValueCount == 2)
        #expect(proteinFeature.missingValueCount == 1)
        #expect(proteinFeature.measuredZeroCount == 1)
        #expect(proteinFeature.minimum == 0)
        #expect(proteinFeature.maximum == 1.25)
        #expect(proteinFeature.mean == 0.625)
        #expect(proteinFeature.variance == 0.78125)
        let metabolite = try #require(summary.features.first { $0.featureID == "HMDB1" })
        #expect(metabolite.measuredValueCount == 1)
        #expect(metabolite.missingValueCount == 2)
        #expect(metabolite.variance == nil)
        #expect(try VivoCanonicalJSON.decode(VivoQuantitativeAssaySummary.self,
                                              from: VivoCanonicalJSON.encode(summary)) == summary)
        try summary.validate()
    }

    @Test func pairwiseAssociationUsesCompleteObservationsAndKeepsZerosMeasured() throws {
        let base = Self.fixture()
        let metabolite = VivoQuantitativeAssaySpace(
            id: "metabolite", kind: .metabolomics, featureNamespace: "hmdb",
            unit: "log10-concentration", sourceDescription: "Measured metabolite panel",
            features: [
                .init(id: "HMDB1", name: "metabolite-1", interval: nil),
                .init(id: "HMDB0", name: "constant-zero", interval: nil)
            ], observationIndices: [0, 1, 2],
            matrix: .init(observationCount: 3, featureCount: 2,
                          rowOffsets: [0, 2, 4, 6],
                          featureIndices: [0, 1, 0, 1, 0, 1],
                          values: [0, 1, 4, 1, 8, 1]))
        let dataset = VivoQuantitativeAssayDataset(
            id: base.id, evidence: base.evidence, sourceDescription: base.sourceDescription,
            samples: base.samples, observations: base.observations,
            spatialFrames: base.spatialFrames, assays: [base.assays[0], metabolite])
        let pair = VivoQuantitativeAssayFeaturePair(leftAssayID: "protein", leftFeatureID: "P1",
                                                     rightAssayID: "metabolite", rightFeatureID: "HMDB1")
        let constantPair = VivoQuantitativeAssayFeaturePair(leftAssayID: "protein", leftFeatureID: "P1",
                                                             rightAssayID: "metabolite", rightFeatureID: "HMDB0")
        let summary = try dataset.pairwiseAssociationSummary([pair, constantPair])
        let result = try #require(summary.associations.first { $0.rightFeatureID == "HMDB1" })
        #expect(result.overlapCount == 2)
        #expect(result.missingPairCount == 1)
        #expect(result.leftMean == 0.625)
        #expect(result.rightMean == 4)
        #expect(result.leftVariance == 0.78125)
        #expect(result.rightVariance == 32)
        #expect(result.covariance == 5)
        #expect(result.pearsonCorrelation == 1)
        #expect(result.status == .computed)
        let constant = try #require(summary.associations.first { $0.rightFeatureID == "HMDB0" })
        #expect(constant.overlapCount == 2)
        #expect(constant.status == .zeroVariance)
        #expect(constant.rightVariance == 0)
        #expect(constant.pearsonCorrelation == nil)
        try summary.validate()
        #expect(try VivoCanonicalJSON.decode(VivoQuantitativeAssayAssociationSummary.self,
                                              from: VivoCanonicalJSON.encode(summary)) == summary)
    }

    @Test func pairwiseAssociationRetainsInsufficientOverlapAndRejectsDuplicatePairs() throws {
        let pair = VivoQuantitativeAssayFeaturePair(leftAssayID: "protein", leftFeatureID: "P1",
                                                     rightAssayID: "metabolite", rightFeatureID: "HMDB1")
        let insufficient = try Self.fixture().pairwiseAssociationSummary([pair])
        let result = try #require(insufficient.associations.first)
        #expect(result.overlapCount == 0)
        #expect(result.missingPairCount == 3)
        #expect(result.leftMean == nil)
        #expect(result.rightMean == nil)
        #expect(result.status == .insufficientOverlap)

        let reversed = VivoQuantitativeAssayFeaturePair(leftAssayID: pair.rightAssayID,
                                                         leftFeatureID: pair.rightFeatureID,
                                                         rightAssayID: pair.leftAssayID,
                                                         rightFeatureID: pair.leftFeatureID)
        #expect(throws: (any Error).self) {
            try Self.fixture().pairwiseAssociationSummary([pair, reversed])
        }
        var object = try #require(JSONSerialization.jsonObject(
            with: VivoCanonicalJSON.encode(insufficient)) as? [String: Any])
        object["unexpected"] = true
        #expect(throws: (any Error).self) {
            try VivoCanonicalJSON.decode(VivoQuantitativeAssayAssociationSummary.self,
                                         from: JSONSerialization.data(withJSONObject: object))
        }
    }

    @Test func rejectsNonfiniteValuesUnsortedIndicesAndDuplicateRows() throws {
        let base = Self.fixture()
        let badValue = VivoSparseValues(observationCount: 1, featureCount: 1,
                                        rowOffsets: [0, 1], featureIndices: [0], values: [.nan])
        let nonfinite = VivoQuantitativeAssaySpace(
            id: "bad", kind: .proteomics, featureNamespace: "p", unit: "intensity",
            sourceDescription: "invalid", features: [.init(id: "p", name: "p", interval: nil)],
            observationIndices: [0], matrix: badValue)
        let badDataset = VivoQuantitativeAssayDataset(id: "bad", evidence: .synthetic,
            sourceDescription: "invalid", samples: base.samples, observations: base.observations,
            spatialFrames: base.spatialFrames, assays: [nonfinite])
        #expect(throws: (any Error).self) { try badDataset.validate() }

        let duplicate = VivoQuantitativeAssaySpace(
            id: "bad", kind: .metabolomics, featureNamespace: "m", unit: "intensity",
            sourceDescription: "invalid", features: [.init(id: "m", name: "m", interval: nil)],
            observationIndices: [0],
            matrix: .init(observationCount: 1, featureCount: 1,
                          rowOffsets: [0, 2], featureIndices: [0, 0], values: [1, 2]))
        let duplicateDataset = VivoQuantitativeAssayDataset(id: "bad", evidence: .synthetic,
            sourceDescription: "invalid", samples: base.samples, observations: base.observations,
            spatialFrames: base.spatialFrames, assays: [duplicate])
        #expect(throws: (any Error).self) { try duplicateDataset.validate() }
    }

    @Test func rejectsUnmeasuredObservationReferencesAndEmptyUnits() throws {
        var dataset = Self.fixture()
        let badRows = VivoQuantitativeAssaySpace(
            id: "bad", kind: .spatialImaging, featureNamespace: "pixels", unit: "intensity",
            sourceDescription: "invalid", features: [.init(id: "px", name: "px", interval: nil)],
            observationIndices: [3],
            matrix: .init(observationCount: 1, featureCount: 1,
                          rowOffsets: [0, 1], featureIndices: [0], values: [1]))
        dataset = .init(id: dataset.id, evidence: dataset.evidence,
                        sourceDescription: dataset.sourceDescription, samples: dataset.samples,
                        observations: dataset.observations, spatialFrames: dataset.spatialFrames,
                        assays: [badRows])
        #expect(throws: (any Error).self) { try dataset.validate() }

        let emptyUnit = VivoQuantitativeAssaySpace(
            id: "bad", kind: .proteomics, featureNamespace: "p", unit: "",
            sourceDescription: "invalid", features: [.init(id: "p", name: "p", interval: nil)],
            observationIndices: [0],
            matrix: .init(observationCount: 1, featureCount: 1,
                          rowOffsets: [0, 1], featureIndices: [0], values: [1]))
        let invalid = VivoQuantitativeAssayDataset(id: "bad", evidence: .synthetic,
            sourceDescription: "invalid", samples: Self.fixture().samples,
            observations: Self.fixture().observations, spatialFrames: Self.fixture().spatialFrames,
            assays: [emptyUnit])
        #expect(throws: (any Error).self) { try invalid.validate() }
    }
}
