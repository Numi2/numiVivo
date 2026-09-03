import Foundation

struct VivoArenaCheckpointData: Sendable {
    let state: [Float]
    let parameters: [Float]
    let temporalState: [Float]
    let transport: [VivoSpeciesTransportABI]
}

public enum VivoTransportRecordLE {
    public static let stride = 16

    public static func encode(_ values: [VivoSpeciesTransportABI]) -> Data {
        var data = Data(capacity: values.count * stride)
        for value in values {
            append(value.diffusion.bitPattern, to: &data)
            append(value.membranePermeability.bitPattern, to: &data)
            append(value.decayRate.bitPattern, to: &data)
            append(value.flags, to: &data)
        }
        return data
    }

    public static func decode(_ data: Data) throws -> [VivoSpeciesTransportABI] {
        guard data.count.isMultiple(of: stride) else {
            throw VivoArtifactValidationError.invalid("transport payload length is not divisible by 16")
        }
        var result: [VivoSpeciesTransportABI] = []
        result.reserveCapacity(data.count / stride)
        var offset = 0
        while offset < data.count {
            let diffusion = Float(bitPattern: try readUInt32(data, offset: offset))
            let permeability = Float(bitPattern: try readUInt32(data, offset: offset + 4))
            let decay = Float(bitPattern: try readUInt32(data, offset: offset + 8))
            let flags = try readUInt32(data, offset: offset + 12)
            guard diffusion.isFinite, permeability.isFinite, decay.isFinite else {
                throw VivoArtifactValidationError.invalid("transport payload contains a non-finite coefficient")
            }
            result.append(.init(
                diffusion: diffusion,
                membranePermeability: permeability,
                decayRate: decay,
                flags: flags
            ))
            offset += stride
        }
        return result
    }

    private static func append(_ value: UInt32, to data: inout Data) {
        var littleEndian = value.littleEndian
        withUnsafeBytes(of: &littleEndian) { data.append(contentsOf: $0) }
    }

    private static func readUInt32(_ data: Data, offset: Int) throws -> UInt32 {
        guard offset >= 0, offset <= data.count - 4 else {
            throw VivoArtifactValidationError.invalid("transport payload is truncated")
        }
        return data.withUnsafeBytes { raw in
            var value: UInt32 = 0
            memcpy(&value, raw.baseAddress!.advanced(by: offset), 4)
            return UInt32(littleEndian: value)
        }
    }
}
