import Foundation
import NumiVivoKit

@main struct LegacyClusteringMain {
    static func main() throws {
        let args = CommandLine.arguments
        guard args.count == 4 else { throw VivoOmicsError.invalid("expected graph JSON, clustering options JSON, output JSON") }
        let graph = try VivoCanonicalJSON.decode(VivoSingleCellNeighborGraph.self, from: VivoSingleCellCampaignIO.readDocument(URL(fileURLWithPath: args[1]), maximumBytes: 536_870_912))
        guard graph.neighborIndices.count <= 4_000_000 else { throw VivoOmicsError.limit("legacy oracle graph bound") }
        let options = try VivoCanonicalJSON.decode(VivoSingleCellClusteringOptions.self, from: Data(contentsOf: URL(fileURLWithPath: args[2])))
        let result = try VivoLegacyClusteringReference.run(graph, options: options)
        try VivoCanonicalJSON.encode(result).write(to: URL(fileURLWithPath: args[3]), options: .withoutOverwriting)
    }
}
