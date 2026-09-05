import Foundation

public indirect enum VivoJSONValue: Codable, Sendable, Hashable {
    case null
    case boolean(Bool)
    case integer(Int64)
    case number(Double)
    case string(String)
    case array([VivoJSONValue])
    case object([String: VivoJSONValue])

    public init(from decoder: Decoder) throws {
        let container = try decoder.singleValueContainer()
        if container.decodeNil() {
            self = .null
        } else if let value = try? container.decode(Bool.self) {
            self = .boolean(value)
        } else if let value = try? container.decode(Int64.self) {
            self = .integer(value)
        } else if let value = try? container.decode(Double.self) {
            self = .number(value)
        } else if let value = try? container.decode(String.self) {
            self = .string(value)
        } else if let value = try? container.decode([VivoJSONValue].self) {
            self = .array(value)
        } else if let value = try? container.decode([String: VivoJSONValue].self) {
            self = .object(value)
        } else {
            throw DecodingError.dataCorruptedError(
                in: container,
                debugDescription: "Unsupported JSON value."
            )
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.singleValueContainer()
        switch self {
        case .null:
            try container.encodeNil()
        case .boolean(let value):
            try container.encode(value)
        case .integer(let value):
            try container.encode(value)
        case .number(let value):
            guard value.isFinite else {
                throw EncodingError.invalidValue(
                    value,
                    .init(
                        codingPath: encoder.codingPath,
                        debugDescription: "JSON numbers must be finite."
                    )
                )
            }
            try container.encode(value)
        case .string(let value):
            try container.encode(value)
        case .array(let value):
            try container.encode(value)
        case .object(let value):
            try container.encode(value)
        }
    }

    public var foundationValue: Any {
        switch self {
        case .null:
            return NSNull()
        case .boolean(let value):
            return value
        case .integer(let value):
            return value
        case .number(let value):
            return value
        case .string(let value):
            return value
        case .array(let values):
            return values.map(\.foundationValue)
        case .object(let values):
            return values.mapValues(\.foundationValue)
        }
    }

    public static func decode(data: Data) throws -> VivoJSONValue {
        try JSONDecoder().decode(VivoJSONValue.self, from: data)
    }

    public func canonicalData() throws -> Data {
        guard JSONSerialization.isValidJSONObject(foundationValue) else {
            if case .null = self { return Data("null".utf8) }
            // JSON serialization is a Foundation concern, independent of Metal
            // runtime construction. Preserve the existing permitted root shapes.
            throw EncodingError.invalidValue(
                self, .init(codingPath: [], debugDescription: "Value is not valid JSON.")
            )
        }
        return try JSONSerialization.data(
            withJSONObject: foundationValue,
            options: [.sortedKeys, .withoutEscapingSlashes]
        )
    }
}
