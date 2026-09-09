import Foundation
@testable import NumiVivoKit

@main struct MNNInvariantMain {
    static func main() throws {
        let samples = (0..<2).map { VivoOmicsSample(id: "s\($0)", biologicalReplicateID: "d\($0)", donorID: "d\($0)", condition: "shared", batchID: "unreported", organism: "numerical-fixture") }
        let cells = (0..<6).map { VivoOmicsCellIdentity(sampleID: "s\($0%2)", barcode: "c\($0)") }
        var options = VivoMNNIntegrationOptions(); options.neighbors = 1
        let tied = try VivoMNNIntegration.run(cells: cells, scores: (0..<6).flatMap { _ in [1.0,0.0] }, dimensions: 2, samples: samples, options: options)
        guard tied.anchors.count == 1, tied.anchors[0].first == 0, tied.anchors[0].second == 1 else { throw VivoOmicsError.invalid("MNN grouped source-index tie order") }
        let spread = (0..<4).map { VivoOmicsCellIdentity(sampleID: "s\($0/2)", barcode: "c\($0)") }
        options.sigma = 1000
        let x = [1.0,0,2,0,1.1,0,1e9,0]
        let underflow = try VivoMNNIntegration.run(cells: spread, scores: x, dimensions: 2, samples: samples, options: options)
        guard underflow.report.steps.count == 1, underflow.report.steps[0].zeroWeightCells == 1,
              abs(underflow.scores[6]-x[6]) <= 1e-6,
              abs(underflow.scores[4]-x[0]) < 1e-12 else { throw VivoOmicsError.invalid("MNN zero-weight preservation") }
        var rejected = 0
        for mode in 0..<4 {
            var o = options
            if mode == 0 { o.maximumWork = 8 }
            if mode == 1 { o.maximumResidentBytes = try VivoMNNIntegration.admittedBytes(rows: 4, dimensions: 2, options: options) }
            if mode == 2 { o.sigma = .nan }
            let identities = mode == 3 ? [spread[0],spread[0],spread[2],spread[3]] : spread
            do { _ = try VivoMNNIntegration.run(cells: identities, scores: x, dimensions: 2, samples: samples, options: o) }
            catch { rejected += 1 }
        }
        guard rejected == 4 else { throw VivoOmicsError.invalid("MNN late budget or input rejection") }
        print("{\"status\":\"passed\",\"deterministicTies\":true,\"zeroWeightRowsPreserved\":true,\"dynamicWorkAndMemoryLimitsReject\":true,\"nonfiniteOptionsAndDuplicateIdentitiesReject\":true,\"rejections\":4,\"qualification\":\"Numerical invariants only\"}")
    }
}
