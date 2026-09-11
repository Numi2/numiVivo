import Foundation
@testable import NumiVivoKit

@main struct WriterCheck {
    static func main() throws {
        let root = URL(fileURLWithPath: CommandLine.arguments[1])
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        var results: [[String: Int]] = []
        for count in [0, 1, 65_535, 65_536, 65_537, 131_075] {
            let url = root.appendingPathComponent("\(count).bin")
            let writer = try VivoCountRecordWriter(url), hashOnly = try VivoCountRecordWriter(nil)
            for (row, feature) in [(-1, 0), (0, -1), (Int(UInt32.max) + 1, 0), (0, Int(UInt32.max) + 1)] {
                do { try writer.append(row: row, feature: feature, bits: 1); fatalError("invalid coordinate accepted") }
                catch let e as VivoOmicsError { guard case .limit = e else { throw e } }
            }
            guard writer.entries == 0 else { fatalError("rejection changed writer") }
            for i in 0..<count {
                let bits = UInt64.max - UInt64(i) * 17
                try writer.append(row: i, feature: Int(UInt32.max) - i, bits: bits)
                try hashOnly.append(row: i, feature: Int(UInt32.max) - i, bits: bits)
            }
            let stored = try writer.finish(), reconstructed = try hashOnly.finish()
            guard stored == reconstructed, writer.entries == count, hashOnly.entries == count,
                  try VivoOmicsFileSnapshot.fingerprint(url, maximumBytes: 4_194_304) == stored else {
                fatalError("writer size, hash or reconstruction differs")
            }
            do { _ = try VivoCountRecordWriter(url); fatalError("overwrite accepted") }
            catch let e as VivoOmicsError { guard case .invalid = e else { throw e } }
            results.append(["entries": count, "invalidCoordinatesRejected": 4, "overwriteRejected": 1])
        }
        // A failed flush must never allow the next append past allocated storage.
        let closed = try VivoCountRecordWriter(root.appendingPathComponent("closed.bin"))
        _ = try closed.finish()
        var flushRejected = false
        do { for i in 0..<65_536 { try closed.append(row: i, feature: 0, bits: 1) } }
        catch { flushRejected = true }
        guard flushRejected else { fatalError("closed output accepted flush") }
        do { try closed.append(row: 0, feature: 0, bits: 1); fatalError("failed flush allowed append") }
        catch let e as VivoOmicsError { guard case .limit = e else { throw e } }
        let data = try JSONSerialization.data(withJSONObject: ["status": "passed", "cases": results, "failedFlushAppendRejected": true], options: [.sortedKeys, .prettyPrinted])
        try data.write(to: root.appendingPathComponent("checks.json"));print(String(decoding: data, as: UTF8.self))
    }
}
