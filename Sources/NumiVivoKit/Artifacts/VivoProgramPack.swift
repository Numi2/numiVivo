import Foundation

public enum VivoProgramPackError: Error, Sendable, CustomStringConvertible {
    case truncated(String)
    case invalidMagic
    case unsupportedVersion(major: UInt16, minor: UInt16)
    case incompatibleCompilerABI(UInt32)
    case invalidHeader(String)
    case invalidSection(String)
    case missingSection(VivoProgramPack.SectionKind)
    case invalidRuntimeContract(String)
    case invalidNumericValue(String)
    public var description: String {
        switch self {
        case .truncated(let field): return "ProgramPack is truncated at \(field)."
        case .invalidMagic: return "ProgramPack magic is invalid."
        case .unsupportedVersion(let major, let minor): return "Unsupported ProgramPack \(major).\(minor)."
        case .incompatibleCompilerABI(let value): return "Unsupported compiler ABI \(value)."
        case .invalidHeader(let value): return "Invalid ProgramPack header: \(value)"
        case .invalidSection(let value): return "Invalid ProgramPack section: \(value)"
        case .missingSection(let value): return "Missing ProgramPack section: \(value)"
        case .invalidRuntimeContract(let value): return "Invalid ProgramPack runtime contract: \(value)"
        case .invalidNumericValue(let value): return "Invalid ProgramPack numeric value: \(value)"
        }
    }
}

/// Construction is the single trust boundary. A pack cannot bypass native
/// semantic/hash validation by being decoded directly instead of via the CLI.
public struct VivoProgramPack: Sendable {
    public static let supportedMajor: UInt16 = 1
    public static let compilerABI: UInt32 = 1
    public static let fixedHeaderSize = 128
    public static let sectionDescriptorSize = 72
    public enum SectionKind: UInt32, CaseIterable, Codable, Sendable, CustomStringConvertible {
        case strings = 1, species, parameters, reactionParameterIndices, stoichiometry
        case reactions, expressions, actions, rules, monitors, cohorts
        case speciesIncidenceOffsets, speciesIncidence, runtimeContract
        public var description: String {
            switch self {
            case .strings: return "strings"
            case .species: return "species"
            case .parameters: return "parameters"
            case .reactionParameterIndices: return "reaction-parameter-indices"
            case .stoichiometry: return "stoichiometry"
            case .reactions: return "reactions"
            case .expressions: return "expressions"
            case .actions: return "actions"
            case .rules: return "rules"
            case .monitors: return "monitors"
            case .cohorts: return "cohorts"
            case .speciesIncidenceOffsets: return "species-incidence-offsets"
            case .speciesIncidence: return "species-incidence"
            case .runtimeContract: return "runtime-contract"
            }
        }
    }
    public struct Header: Sendable, Codable, Equatable {
        public let major: UInt16, minor: UInt16
        public let headerBytes: UInt32, compilerABI: UInt32, flags: UInt32
        public let fidelity: VivoFidelity
        public let sectionCount: UInt32
        public let totalBytes: UInt64
        public let sourceFingerprint: VivoFingerprint, contentFingerprint: VivoFingerprint
    }
    public struct Section: Sendable, Codable, Equatable {
        public let kind: SectionKind
        public let flags: UInt32
        public let offset: UInt64, size: UInt64
        public let stride: UInt32, count: UInt32, alignment: UInt32
        public let fingerprint: VivoFingerprint
        public var range: Range<Int> { Int(offset)..<Int(offset + size) }
    }
    public struct RuntimeContract: Sendable, Codable, Equatable {
        public let speciesCount: UInt32, parameterCount: UInt32, reactionCount: UInt32
        public let ruleCount: UInt32, monitorCount: UInt32, cohortCount: UInt32
        public let temporalStateCount: UInt32, maximumExpressionStack: UInt32
        public let featureFlags: UInt32, authoritativeScalarBytes: UInt32, randomStreamVersion: UInt32
    }
    public struct SpeciesMetadata: Sendable, Codable, Equatable {
        public let identifier: String, compartment: String, unit: String
        public let flags: UInt32
        public let initialValue: Float, minimum: Float, maximum: Float
        public var isExternallyOwned: Bool { flags & (1 << 1) != 0 }
        public var isInput: Bool { flags & (1 << 2) != 0 }
        public var isState: Bool { flags & (1 << 3) != 0 }
        public var isOutput: Bool { flags & (1 << 4) != 0 }
        public var isCountValued: Bool { flags & (1 << 5) != 0 }
    }
    public struct ParameterMetadata: Sendable, Codable, Equatable {
        public let identifier: String, unit: String, evidenceSource: String
        public let flags: UInt32
        public let value: Double, minimum: Double, maximum: Double
        public let evidenceClass: UInt32
    }
    public let data: Data
    public let header: Header
    public let sections: [SectionKind: Section]
    public let runtimeContract: RuntimeContract

    public init(data: Data) throws {
        guard data.count >= Self.fixedHeaderSize, data.count <= 256 * 1024 * 1024 else {
            throw VivoProgramPackError.invalidHeader("pack size must be between 128 bytes and 256 MiB")
        }
        // Data slices may retain a nonzero startIndex. Own zero-based bytes
        // before using wire offsets as collection indices.
        let data = data.withUnsafeBytes { Data($0) }
        let inspection = VivoNativeCompilerBridge.inspectProgramPack(data, verifySectionHashes: true)
        guard inspection.invocation.succeeded else {
            let report = String(decoding: inspection.invocation.primary.prefix(8192), as: UTF8.self)
            throw VivoProgramPackError.invalidRuntimeContract(report)
        }
        let reader = PackReader(data: data)
        guard let fidelity = VivoFidelity(rawValue: try reader.u32(24)) else {
            throw VivoProgramPackError.invalidHeader("unknown fidelity")
        }
        let count = try reader.u32(28)
        var sections: [SectionKind: Section] = [:]
        for index in 0..<Int(count) {
            let base = Self.fixedHeaderSize + index * Self.sectionDescriptorSize
            guard let kind = SectionKind(rawValue: try reader.u32(base)) else { continue }
            sections[kind] = try Section(kind: kind, flags: reader.u32(base + 4),
                                         offset: reader.u64(base + 8), size: reader.u64(base + 16),
                                         stride: reader.u32(base + 24), count: reader.u32(base + 28),
                                         alignment: reader.u32(base + 32), fingerprint: reader.fingerprint(base + 40))
        }
        for kind in SectionKind.allCases where sections[kind] == nil { throw VivoProgramPackError.missingSection(kind) }
        self.data = data
        self.sections = sections
        header = try Header(major: reader.u16(8), minor: reader.u16(10), headerBytes: reader.u32(12),
                             compilerABI: reader.u32(16), flags: reader.u32(20), fidelity: fidelity,
                             sectionCount: count, totalBytes: reader.u64(32),
                             sourceFingerprint: reader.fingerprint(40), contentFingerprint: reader.fingerprint(72))
        let base = Int(sections[.runtimeContract]!.offset)
        runtimeContract = try RuntimeContract(speciesCount: reader.u32(base), parameterCount: reader.u32(base + 4),
                                              reactionCount: reader.u32(base + 8), ruleCount: reader.u32(base + 12),
                                              monitorCount: reader.u32(base + 16), cohortCount: reader.u32(base + 20),
                                              temporalStateCount: reader.u32(base + 24), maximumExpressionStack: reader.u32(base + 28),
                                              featureFlags: reader.u32(base + 32), authoritativeScalarBytes: reader.u32(base + 36),
                                              randomStreamVersion: reader.u32(base + 40))
    }
    public func section(_ kind: SectionKind) throws -> Section {
        guard let value = sections[kind] else { throw VivoProgramPackError.missingSection(kind) }
        return value
    }
    public func sectionData(_ kind: SectionKind) throws -> Data { try data.subdata(in: section(kind).range) }
    public func speciesMetadata() throws -> [SpeciesMetadata] {
        let records = try section(.species)
        let reader = PackReader(data: data)
        return try (0..<Int(records.count)).map { index in
            let base = Int(records.offset) + index * Int(records.stride)
            return try SpeciesMetadata(identifier: string(reader.u32(base)), compartment: string(reader.u32(base + 4)),
                                        unit: string(reader.u32(base + 8)), flags: reader.u32(base + 12),
                                        initialValue: reader.f32(base + 16), minimum: reader.f32(base + 20), maximum: reader.f32(base + 24))
        }
    }
    public func parameterMetadata() throws -> [ParameterMetadata] {
        let records = try section(.parameters)
        let reader = PackReader(data: data)
        return try (0..<Int(records.count)).map { index in
            let base = Int(records.offset) + index * Int(records.stride)
            return try ParameterMetadata(identifier: string(reader.u32(base)), unit: string(reader.u32(base + 4)),
                                          evidenceSource: string(reader.u32(base + 8)), flags: reader.u32(base + 12),
                                          value: reader.f64(base + 16), minimum: reader.f64(base + 24),
                                          maximum: reader.f64(base + 32), evidenceClass: reader.u32(base + 40))
        }
    }
    public func initialState(laneCount: Int) throws -> [Float] {
        let metadata = try speciesMetadata()
        var values = [Float](repeating: 0, count: try Self.elements(metadata.count, laneCount))
        for (index, item) in metadata.enumerated() {
            values.replaceSubrange((index * laneCount)..<((index + 1) * laneCount),
                                   with: repeatElement(item.initialValue, count: laneCount))
        }
        return values
    }
    public func parameterValues(environmentCount: Int) throws -> [Float] {
        let metadata = try parameterMetadata()
        var values = [Float](repeating: 0, count: try Self.elements(metadata.count, environmentCount))
        for (index, item) in metadata.enumerated() {
            let value = Float(item.value)
            guard value.isFinite, item.value == 0 || value != 0 else {
                throw VivoProgramPackError.invalidNumericValue("parameter \(item.identifier) cannot be represented in FP32")
            }
            values.replaceSubrange((index * environmentCount)..<((index + 1) * environmentCount),
                                   with: repeatElement(value, count: environmentCount))
        }
        return values
    }
    public func speciesIndex(named identifier: String) throws -> UInt32? {
        try speciesMetadata().firstIndex { $0.identifier == identifier }.map(UInt32.init)
    }
    private func string(_ offset: UInt32) throws -> String {
        let strings = try section(.strings)
        guard UInt64(offset) < strings.size else { throw VivoProgramPackError.invalidSection("string offset") }
        let start = Int(strings.offset) + Int(offset)
        let end = Int(strings.offset + strings.size)
        let suffix = data[start..<end]
        guard let terminator = suffix.firstIndex(of: 0),
              let result = String(data: data[start..<terminator], encoding: .utf8) else {
            throw VivoProgramPackError.invalidSection("unterminated or invalid UTF-8 string")
        }
        return result
    }
    private static func elements(_ rows: Int, _ columns: Int) throws -> Int {
        let value = rows.multipliedReportingOverflow(by: columns)
        guard rows >= 0, columns > 0, !value.overflow, value.partialValue <= Int.max / MemoryLayout<Float>.stride else {
            throw VivoProgramPackError.invalidNumericValue("state shape overflows host address space")
        }
        return value.partialValue
    }
}
private struct PackReader {
    let data: Data
    private func integer<T: FixedWidthInteger>(_ offset: Int, _: T.Type) throws -> T {
        guard offset >= 0, offset <= data.count, MemoryLayout<T>.size <= data.count - offset else {
            throw VivoProgramPackError.truncated("offset \(offset)")
        }
        return data.withUnsafeBytes { T(littleEndian: $0.loadUnaligned(fromByteOffset: offset, as: T.self)) }
    }
    func u16(_ offset: Int) throws -> UInt16 { try integer(offset, UInt16.self) }
    func u32(_ offset: Int) throws -> UInt32 { try integer(offset, UInt32.self) }
    func u64(_ offset: Int) throws -> UInt64 { try integer(offset, UInt64.self) }
    func f32(_ offset: Int) throws -> Float { try Float(bitPattern: u32(offset)) }
    func f64(_ offset: Int) throws -> Double { try Double(bitPattern: u64(offset)) }
    func fingerprint(_ offset: Int) throws -> VivoFingerprint {
        guard offset >= 0, offset <= data.count, 32 <= data.count - offset else {
            throw VivoProgramPackError.truncated("fingerprint")
        }
        return try VivoFingerprint(bytes: Array(data[offset..<(offset + 32)]))
    }
}
