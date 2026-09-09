import Foundation
import NumiVivoKit

@main struct LegacyIntegrationMain {
    static func main() throws {
        let args = CommandLine.arguments
        guard args.count == 5 else { throw VivoOmicsError.invalid("expected PCA bundle, dimensions, integration options JSON, output JSON") }
        let root = URL(fileURLWithPath: args[1])
        let metadata = try VivoCanonicalJSON.decode(VivoSingleCellCountMetadata.self, from: Data(contentsOf: root.appendingPathComponent("metadata.json")))
        guard let dimensions = Int(args[2]), (1...64).contains(dimensions), metadata.cells.count <= 100_000 else { throw VivoOmicsError.invalid("legacy oracle axes") }
        let raw = try Data(contentsOf: root.appendingPathComponent("scores.bin"))
        guard raw.count == metadata.cells.count * dimensions * 16 else { throw VivoOmicsError.invalid("legacy score size") }
        var scores = [[Double]](repeating: [], count: metadata.cells.count)
        try raw.withUnsafeBytes { bytes in
            for row in scores.indices { for column in 0..<dimensions {
                let at = (row * dimensions + column) * 16
                let packed = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: at, as: UInt64.self))
                let value = Double(bitPattern: UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: at + 8, as: UInt64.self)))
                guard Int(packed & 0xffff_ffff) == row, Int(packed >> 32) == column, value.isFinite else { throw VivoOmicsError.invalid("legacy score coordinates") }
                scores[row].append(value)
            } }
        }
        let options = try VivoCanonicalJSON.decode(VivoSingleCellIntegrationOptions.self, from: Data(contentsOf: URL(fileURLWithPath: args[3])))
        let input = VivoLegacyIntegrationInput(cells: metadata.cells.map { .init(sampleID: $0.sampleID, barcode: $0.barcode) }, scores: scores)
        let result = try VivoLegacyIntegrationReference.run(input, samples: metadata.samples, options: options)
        try VivoCanonicalJSON.encode(result).write(to: URL(fileURLWithPath: args[4]), options: .withoutOverwriting)
    }
}
