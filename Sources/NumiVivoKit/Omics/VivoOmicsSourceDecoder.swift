import Foundation
import CNumiVivoZlib

public enum VivoOmicsCompressionError: Error, LocalizedError, Sendable {
    case invalid(String)
    public var errorDescription: String? { switch self { case .invalid(let value): return "omics compression: \(value)" } }
}
/// Gzip is detected by its header, never by an untrusted filename. CRC/trailer
/// validation, concatenated members and total expanded-byte bounds are enforced.
/// The campaign keeps compressed originals in its existing input snapshot.
public enum VivoOmicsSourceDecoder {
    public static func decode(_ input: Data, maximumExpandedBytes: Int) throws -> Data {
        guard maximumExpandedBytes >= 0, input.count <= Int(UInt32.max) else {
            throw VivoOmicsCompressionError.invalid("invalid expansion limit or oversized compressed input")
        }
        let magic = input.prefix(2)
        guard magic.elementsEqual([UInt8(0x1f), 0x8b]) else {
            guard input.count <= maximumExpandedBytes else { throw VivoOmicsCompressionError.invalid("plain input exceeds expanded-byte allowance") }
            return input
        }
        return try input.withUnsafeBytes { source in
            var stream = z_stream()
            guard inflateInit2_(&stream, 31, ZLIB_VERSION, Int32(MemoryLayout<z_stream>.size)) == Z_OK else {
                throw VivoOmicsCompressionError.invalid("could not initialize native gzip decoder")
            }
            defer { inflateEnd(&stream) }
            stream.next_in = UnsafeMutablePointer(mutating: source.bindMemory(to: Bytef.self).baseAddress)
            stream.avail_in = uInt(input.count)
            var output = Data(), chunk = [UInt8](repeating: 0, count: 65_536), members = 0
            while true {
                try Task.checkCancellation()
                let before = stream.avail_in
                stream.avail_out = uInt(chunk.count)
                let code = chunk.withUnsafeMutableBytes { buffer -> Int32 in
                    stream.next_out = buffer.bindMemory(to: Bytef.self).baseAddress
                    return inflate(&stream, Z_NO_FLUSH)
                }
                let produced = chunk.count - Int(stream.avail_out)
                guard produced <= maximumExpandedBytes - output.count else {
                    throw VivoOmicsCompressionError.invalid("gzip expansion exceeds aggregate source limit")
                }
                output.append(contentsOf: chunk.prefix(produced))
                if code == Z_STREAM_END {
                    members += 1
                    guard members <= 1024 else { throw VivoOmicsCompressionError.invalid("too many concatenated gzip members") }
                    if stream.avail_in == 0 { return output }
                    guard stream.avail_in >= 2, let next = stream.next_in, next[0] == 0x1f, next[1] == 0x8b else {
                        throw VivoOmicsCompressionError.invalid("trailing non-gzip data")
                    }
                    let remaining = stream.avail_in
                    guard inflateReset2(&stream, 31) == Z_OK else { throw VivoOmicsCompressionError.invalid("gzip member reset failed") }
                    stream.avail_in = remaining; stream.next_in = next
                } else {
                    guard code == Z_OK, produced > 0 || stream.avail_in < before else {
                        throw VivoOmicsCompressionError.invalid("truncated or corrupt gzip stream (native code \(code))")
                    }
                }
            }
        }
    }
}
