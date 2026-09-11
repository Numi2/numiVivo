import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellDurationTests {
    static let source = try! VivoFingerprint(bytes: Array(repeating: 0, count: 32))
    static func mapping() -> VivoH5ADImportPlan {
        .init(id: "duration-fixture", evidence: .synthetic, sourceDescription: "Numerical control only", countUnit: .umiCount, matrixPath: "X", samples: [], sampleColumn: "sample")
    }
    static func plan(panel: [String]? = nil, exposures: [VivoDurationExposure]? = nil) -> VivoDurationPlan {
        .init(mapping: mapping(), featureNamespace: "fixture", perturbationID: "stim", controlCondition: "control",
              exposures: exposures ?? [.init(condition: "early", hours: 1), .init(condition: "late", hours: 8)],
              responseFeatureIDs: panel, provenance: "Numerical fixture, no biological claim")
    }
    static func bulk(_ values: [[UInt64]], donors: [String], conditions: [String], features: [String] = ["A", "B", "C"]) -> VivoPseudobulkCounts {
        let groups = values.indices.map { i in
            VivoPseudobulkGroup(biologicalReplicateID: donors[i], donorID: donors[i], condition: conditions[i], organism: "human", cellGroup: "PBMC",
                               sampleIDs: ["s\(i)"], batchIDs: [], sourceCellIndices: [i])
        }
        return .init(method: "fixture", countUnit: .umiCount, groups: groups, featureIDs: features,
                     matrix: .init(cellCount: values.count, featureCount: features.count,
                                   rowOffsets: (0...values.count).map { $0 * features.count },
                                   featureIndices: values.flatMap { _ in Array(features.indices) }, counts: values.flatMap { $0 }))
    }
    static func training() -> VivoPseudobulkCounts {
        bulk([[10,20,70],[30,30,40],[60,20,20],[20,30,50],[40,20,40],[50,40,10]],
             donors: ["a","a","a","b","b","b"], conditions: ["control","early","late","control","early","late"])
    }
    static func query(_ hours: [Double], donor: String = "q", condition: String = "control", features: [String] = ["A","B","C"], counts: [UInt64] = [15,25,60], model: VivoDurationModel) throws -> VivoDurationReport {
        let p = VivoDurationQueryPlan(mapping: mapping(), featureNamespace: "fixture", perturbationID: "stim", hours: hours)
        return try VivoDurationPerturbation.evaluate(bulk([counts], donors: [donor], conditions: [condition], features: features), plan: p, model: model)
    }
    @Test func observedKnotsZeroAnchorAndLogTimeInterpolation() throws {
        let c = VivoDurationCurve(donorID: "a", hours: [0,1,8], responses: [[0,0],[2,-4],[6,8]])
        #expect(try VivoDurationPerturbation.interpolate(c, hours: 0) == [0,0])
        #expect(try VivoDurationPerturbation.interpolate(c, hours: 1) == [2,-4])
        #expect(try VivoDurationPerturbation.interpolate(c, hours: 8) == [6,8])
        let midpoint = sqrt(2.0 * 9.0) - 1
        let result = try VivoDurationPerturbation.interpolate(c, hours: midpoint)
        #expect(abs(result[0]-4) < 1e-12 && abs(result[1]-2) < 1e-12)
        for t in [-1.0, 9, Double.nan, Double.infinity] {
            #expect(throws: (any Error).self) { try VivoDurationPerturbation.interpolate(c, hours: t) }
        }
    }
    @Test func fitsDonorCurvesAndPredictionIsEqualDonorWeighted() throws {
        let b = Self.training(), p = Self.plan(), model = try VivoDurationPerturbation.model(b, plan: p, source: Self.source)
        #expect(model == (try VivoDurationPerturbation.model(b, plan: p, source: Self.source)))
        #expect(model.maximumSolveResidual < 1e-12 && model.maximumSupportedHours == 8)
        let (_, logs, _) = try VivoPerturbation.dense(b)
        for j in 0..<3 {
            #expect(model.curves[0].responses[1][j] == logs[1][j]-logs[0][j])
            #expect(model.curves[1].responses[2][j] == logs[5][j]-logs[3][j])
        }
        let report = try Self.query([0,1,3,8], model: model)
        let zero = report.predictions[0].prediction
        #expect(zero.estimates.allSatisfy { $0.predictedTreated == zero.control && $0.unclippedResponse == [0,0,0] })
        #expect(zero.meanResponsePredictiveInterval?.unavailableFeatureIndices == [0,1,2])
        let early = report.predictions[1].prediction, mean = early.estimates[1]
        for j in 0..<3 { #expect(abs(mean.unclippedResponse[j]-(logs[1][j]-logs[0][j]+logs[4][j]-logs[3][j])/2) < 1e-12) }
        let interval = try #require(early.meanResponsePredictiveInterval)
        #expect(interval.degreesOfFreedom == 1 && abs(interval.studentCriticalValue-12.7062047361747) < 1e-9)
    }
    @Test func unequalSchedulesAreInterpolatedBeforeDonorAverage() throws {
        let b = Self.bulk([[10,20,70],[30,30,40],[60,20,20],[20,30,50],[40,20,40],[50,40,10],[70,20,10]],
                          donors: ["a","a","a","b","b","b","b"], conditions: ["control","early","late","control","early","middle","late"])
        let p = Self.plan(exposures: [.init(condition: "early", hours: 1), .init(condition: "middle", hours: 3), .init(condition: "late", hours: 8)])
        let model = try VivoDurationPerturbation.model(b, plan: p, source: Self.source)
        let a = try VivoDurationPerturbation.interpolate(model.curves[0], hours: 3), c = model.curves[1].responses[2]
        let mean = try Self.query([3], model: model).predictions[0].prediction.estimates[1].unclippedResponse
        for j in mean.indices { #expect(abs(mean[j]-(a[j]+c[j])/2) < 1e-12) }
    }
    @Test func panelRetainsFullDenominatorsAndFeatureOrdering() throws {
        let model = try VivoDurationPerturbation.model(Self.training(), plan: Self.plan(panel: ["B","A"]), source: Self.source)
        let a = try Self.query([1], model: model), b = try Self.query([1], counts: [15,25,960], model: model)
        #expect(a.predictions[0].prediction.libraryCounts == 100 && b.predictions[0].prediction.libraryCounts == 1000)
        #expect(abs(a.predictions[0].prediction.control[0]-log1p(250_000.0)) < 1e-12)
        #expect(abs(b.predictions[0].prediction.control[0]-log1p(25_000.0)) < 1e-12)
        #expect(try Self.query([1], features: ["C","B","A"], counts: [60,25,15], model: model) == a)
        #expect(throws: (any Error).self) { try Self.query([1], features: ["C","A"], counts: [60,15], model: model) }
    }
    @Test func rejectsUnqualifiedQueriesAndPlans() throws {
        let model = try VivoDurationPerturbation.model(Self.training(), plan: Self.plan(), source: Self.source)
        for hours in [[9.0],[-1.0],[1.0,1.0],[],[Double.nan]] {
            #expect(throws: (any Error).self) { try Self.query(hours, model: model) }
        }
        #expect(throws: (any Error).self) { try Self.query([1], donor: "a", model: model) }
        #expect(throws: (any Error).self) { try Self.query([1], condition: "early", model: model) }
        for exposures in [[VivoDurationExposure(condition: "early", hours: 1), .init(condition: "late", hours: 1)],
                          [.init(condition: "early", hours: 0), .init(condition: "late", hours: 8)],
                          [.init(condition: "early", hours: 1), .init(condition: "missing", hours: 8)]] {
            #expect(throws: (any Error).self) { try VivoDurationPerturbation.model(Self.training(), plan: Self.plan(exposures: exposures), source: Self.source) }
        }
        #expect(throws: (any Error).self) { try VivoDurationPerturbation.model(Self.training(), plan: Self.plan(panel: ["missing"]), source: Self.source) }
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoDurationExposure.self, from: Data("{\"condition\":\"early\",\"hours\":1,\"dose\":10}".utf8)) }
    }
}
