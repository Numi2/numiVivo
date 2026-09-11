import Foundation
@testable import NumiVivoKit

@main struct KernelCheck {
    static func main() throws {
        let root = URL(fileURLWithPath: CommandLine.arguments[1])
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        let raw = root.appendingPathComponent("counts.bin")
        let writer = try VivoCountRecordWriter(raw)
        let regular = VivoMetalCountNormalization.batchEntries + 1
        for _ in 0..<regular { try writer.append(row: 0, feature: 0, bits: 1) }
        for (row, count) in [(1, UInt64(1)), (1, UInt64(1) << 53), (2, UInt64(1)), (2, UInt64.max - 1)] {
            try writer.append(row: row, feature: row, bits: count)
        }
        _ = try writer.finish()
        let records = try VivoWindowedCountRecords(raw, entries: regular + 4)
        let quality = VivoCountStoreQuality(rowTotals: [UInt64(regular), (UInt64(1) << 53) + 1, UInt64.max],
            rowNonzeros: [UInt64(regular), 2, 2], featureTotals: [UInt64(regular), (UInt64(1) << 53) + 1, UInt64.max],
            featureNonzeros: [UInt64(regular), 2, 2])
        var results: [[String: String]] = []
        for target in [1.0, 10_000.0, 1_000_000_000.0] {
            var first: VivoFingerprint?
            for attempt in 0..<2 {
                let url = root.appendingPathComponent("values-\(Int(target))-\(attempt).bin")
                let output = try VivoCountRecordWriter(url)
                let execution = try VivoMetalCountNormalization.run(records: records, quality: quality, target: target, writer: output)
                let hash = try output.finish()
                if let first { guard hash == first else { fatalError("GPU repeat differs") } } else { first = hash }
                let normalized = try VivoWindowedCountRecords(url, entries: regular + 4)
                var maximum = 0.0
                for i in 0..<records.count {
                    let r = try records.record(i), o = try normalized.record(i)
                    let expected = log1p(Double(r.bits) * (target / Double(quality.rowTotals[r.row])))
                    let actual = Double(bitPattern: o.bits)
                    guard o.row == r.row, o.feature == r.feature, actual > 0, actual.isFinite,
                          abs(actual - expected) <= 3e-6 + 3e-6 * abs(expected) else { fatalError("kernel mismatch at \(i)") }
                    maximum = max(maximum, abs(actual - expected))
                }
                results.append(["target": String(target), "attempt": String(attempt), "entries": String(records.count),
                                "maximumAbsoluteError": String(maximum), "device": execution.deviceName])
            }
        }
        let encoded = try JSONSerialization.data(withJSONObject: ["status": "passed", "cases": results], options: [.sortedKeys, .prettyPrinted])
        try encoded.write(to: root.appendingPathComponent("checks.json"))
        print(String(data: encoded, encoding: .utf8)!)
    }
}
