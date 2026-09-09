import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct MNNIntegrationTests {
    struct Pair: Decodable { let firstLevel: Int; let secondLevel: Int; let pairs: [[Int]] }
    struct Fixture: Decodable {
        let cells: Int; let components: Int; let neighbors: Int; let scale: Double
        let batch: [Int]; let scores: [[Double]]; let expected: [[Double]]
        let alignmentOrder: [[Int]]; let anchors: [Pair]
    }
    static func fixture() throws -> Fixture {
        let path = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures/MNNIntegration.json")
        return try JSONDecoder().decode(Fixture.self, from: Data(contentsOf: path))
    }
    static func inputs(_ f: Fixture) -> ([VivoOmicsCellIdentity],[VivoOmicsSample],VivoMNNIntegrationOptions) {
        let cells = f.batch.enumerated().map { VivoOmicsCellIdentity(sampleID: "s\($0.element)", barcode: "c\($0.offset)") }
        let samples = (0..<3).map { VivoOmicsSample(id: "s\($0)", biologicalReplicateID: "d\($0)", donorID: "d\($0)", condition: "shared", batchID: "unreported", organism: "human") }
        var options = VivoMNNIntegrationOptions(); options.neighbors = f.neighbors; options.minimumAlignment = 0
        return (cells,samples,options)
    }
    @Test func matchesIndependentReferenceAndUniformScale() throws {
        let f = try Self.fixture(), (cells,samples,o) = Self.inputs(f), x = f.scores.flatMap { $0 }
        let result = try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: f.components, samples: samples, options: o)
        #expect(abs(result.report.medianRowNorm-f.scale)<1e-12)
        let expected = f.expected.flatMap { $0 }
        for (a,b) in zip(result.scores,expected) { #expect(abs(a-b)<1e-10) }
        #expect(result.report.assemblyOrder.map { [result.report.alignments[$0].firstLevel,result.report.alignments[$0].secondLevel] } == f.alignmentOrder)
        for pair in f.anchors {
            let alignment = try #require(result.report.alignments.first { $0.firstLevel == pair.firstLevel && $0.secondLevel == pair.secondLevel })
            let actual = result.anchors[alignment.anchorOffset..<(alignment.anchorOffset+alignment.anchors)].map { [$0.first,$0.second] }
            #expect(actual == pair.pairs)
        }
        let scaled = try VivoMNNIntegration.run(cells: cells, scores: x.map { $0*17 }, dimensions: f.components, samples: samples, options: o)
        for (a,b) in zip(scaled.scores,result.scores) { #expect(abs(a/17-b)<1e-10) }
        #expect(result.report.steps.contains { $0.targetLevels.count>1 })
        #expect(result.report.steps.allSatisfy { $0.zeroWeightCells==0 })
    }
    @Test func rejectsUnavailableMetadataBudgetsAndDegenerateAxes() throws {
        let f = try Self.fixture(), (cells,samples,base) = Self.inputs(f), x = f.scores.flatMap { $0 }
        var o = base; o.maximumWork = 1
        #expect(throws: (any Error).self) { try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: f.components, samples: samples, options: o) }
        o = base; o.maximumResidentBytes = 1
        #expect(throws: (any Error).self) { try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: f.components, samples: samples, options: o) }
        o = base; o.covariate = .batch
        #expect(throws: (any Error).self) { try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: f.components, samples: samples, options: o) }
        #expect(throws: (any Error).self) { try VivoMNNIntegration.run(cells: cells, scores: x.map { _ in 0 }, dimensions: f.components, samples: samples, options: base) }
        let confounded = (0..<3).map { VivoOmicsSample(id: "s\($0)", biologicalReplicateID: "d\($0)", donorID: "d\($0)", condition: "c\($0)", batchID: "unreported", organism: "human") }
        #expect(throws: (any Error).self) { try VivoMNNIntegration.run(cells: cells, scores: x, dimensions: f.components, samples: confounded, options: base) }
    }
    @Test func methodSelectionPreservesLegacyEncodingAndNoForcedMatches() throws {
        let legacy = VivoPCAIntegrationPlan()
        let data = try VivoCanonicalJSON.encode(legacy)
        #expect(!String(decoding: data,as: UTF8.self).contains("mnn"))
        #expect(try JSONDecoder().decode(VivoPCAIntegrationPlan.self, from: data)==legacy)
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoPCAIntegrationPlan.self, from: Data("{\"schemaVersion\":1,\"integration\":{},\"mnn\":{}}".utf8)).validate() }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoMNNIntegrationOptions.self, from: Data("{\"sigm\":15}".utf8)) }
        let f = try Self.fixture(), (cells,samples,base) = Self.inputs(f)
        var o = base; o.minimumAlignment = 1
        let result = try VivoMNNIntegration.run(cells: cells, scores: f.scores.flatMap { $0 }, dimensions: f.components, samples: samples, options: o)
        #expect(result.report.steps.isEmpty && result.report.panoramas.count == 3)
        for (a,b) in zip(result.scores,f.scores.flatMap({ $0 })) { #expect(abs(a-b)<1e-12) }
    }
}
