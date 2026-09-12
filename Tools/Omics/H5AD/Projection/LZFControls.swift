import Foundation

@main struct LZFControls {
    static func main() throws {
        func decode(_ input: [UInt8], _ capacity: Int) -> (Int, [UInt8]) {
            var result = [UInt8](repeating: 0xa5, count: capacity + 16)
            let count = input.withUnsafeBufferPointer { source in
                result.withUnsafeMutableBufferPointer { target in
                    VivoHDF5.decodeLZF(input: source.baseAddress!, count: input.count,
                                      output: target.baseAddress!, capacity: capacity)
                }
            }
            precondition(result.suffix(16).allSatisfy { $0 == 0xa5 })
            return (count, Array(result.prefix(max(0,count))))
        }
        let short = decode([0,65,32,0],4)
        precondition(short.0 == 4 && short.1 == [65,65,65,65])
        let long = decode([0,65,224,255,0],265)
        precondition(long.0 == 265 && long.1 == Array(repeating:65,count:265))
        let invalids: [[UInt8]] = [[0],[31,1],[32,0],[224],[224,0],[0,65,32,255]]
        for invalid in invalids {
            precondition(decode(invalid,512).0 == 0)
        }
        precondition(decode([0,65,32,0],3).0 == -1)
        var state: UInt64 = 0x3194d218
        func random() -> UInt64 {
            state = state &* 6364136223846793005 &+ 1442695040888963407
            return state
        }
        for _ in 0..<20_000 {
            let length = Int(random() % 64) + 1, capacity = Int(random() % 256)
            let input = (0..<length).map { _ in UInt8(truncatingIfNeeded: random() >> 32) }
            let result = decode(input,capacity)
            precondition(result.0 >= -1 && result.0 <= capacity)
        }
        let root = URL(fileURLWithPath: CommandLine.arguments[1])
        let manifest = try JSONSerialization.jsonObject(with: Data(contentsOf: root.appendingPathComponent("chunks.json"))) as! [[String: String]]
        for entry in manifest {
            let input = try Data(contentsOf: root.appendingPathComponent(entry["compressed"]!))
            let expected = try Data(contentsOf: root.appendingPathComponent(entry["expected"]!))
            let result = decode(Array(input),expected.count)
            precondition(result.0 == expected.count && result.1 == Array(expected))
        }
        print("Passed literal, overlapping, long-reference, truncation and bounds controls; 20000 deterministic streams; \(manifest.count) independent h5py compressed chunks")
    }
}
