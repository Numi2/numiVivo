import Foundation
import NumiVivoKit

@main struct LegacyEmbeddingMain {
    static func main() throws {
        let args = CommandLine.arguments
        guard args.count == 5 else { throw VivoOmicsError.invalid("expected graph JSON, scores.bin, embedding options JSON, output JSON") }
        let graph = try VivoCanonicalJSON.decode(VivoSingleCellNeighborGraph.self, from: VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: args[1]), maximumBytes: 536_870_912))
        guard graph.neighborIndices.count <= 4_000_000 else { throw VivoOmicsError.limit("legacy oracle graph bound") }
        let raw = try Data(contentsOf: URL(fileURLWithPath: args[2]))
        guard raw.count == graph.cells.count * graph.dimensions * 16 else { throw VivoOmicsError.invalid("legacy score axes") }
        var scores = [[Double]](repeating: [], count: graph.cells.count)
        try raw.withUnsafeBytes { bytes in
            for row in graph.cells.indices { for column in 0..<graph.dimensions {
                let at = (row * graph.dimensions + column) * 16
                let packed = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: at, as: UInt64.self))
                let value = Double(bitPattern: UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: at + 8, as: UInt64.self)))
                guard Int(packed & 0xffff_ffff) == row, Int(packed >> 32) == column, value.isFinite else { throw VivoOmicsError.invalid("legacy score coordinates") }
                scores[row].append(value)
            } }
        }
        let options = try VivoCanonicalJSON.decode(VivoSingleCellEmbeddingOptions.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
        let result = try VivoLegacyEmbeddingReference.run(graph, scores: scores, options: options)
        try VivoCanonicalJSON.encode(result).write(to: URL(fileURLWithPath: args[4]), options: .withoutOverwriting)
    }
}
