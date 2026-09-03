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
        case .truncated(let field): "ProgramPack is truncated while reading \(field)."
        case .invalidMagic: "ProgramPack magic is invalid."
        case .unsupportedVersion(let major, let minor): "ProgramPack version \(major).\(minor) is unsupported."
        case .incompatibleCompilerABI(let version): "ProgramPack compiler ABI \(version) is incompatible."
        case .invalidHeader(let message): "ProgramPack header is invalid: \(message)"
        case .invalidSection(let message): "ProgramPack section is invalid: \(message)"
        case .missingSection(let section): "ProgramPack is missing required section \(section)."
        case .invalidRuntimeContract(let message): "ProgramPack runtime contract is invalid: \(message)"
        case .invalidNumericValue(let message): "ProgramPack contains an invalid numeric value: \(message)"
        }
    }
}

public struct VivoFingerprint: Hashable, Codable, Sendable, CustomStringConvertible {
    public static let byteCount = 32
    public let bytes: [UInt8]

    public init(bytes: [UInt8]) throws {
        guard bytes.count == Self.byteCount else {
            throw VivoProgramPackError.invalidHeader("fingerprint must contain exactly 32 bytes")
        }
        self.bytes = bytes
    }

    public var hex: String {
        bytes.map { String(format: "%02x", $0) }.joined()
    }

    public var description: String { hex }
}

public struct VivoProgramPack: Sendable {
    public static let supportedMajor: UInt16 = 1
    public static let compilerABI: UInt32 = 1
    public static let fixedHeaderSize = 128
    public static let sectionDescriptorSize = 72

    public enum SectionKind: UInt32, CaseIterable, Codable, Sendable, CustomStringConvertible {
        case strings = 1
        case species = 2
        case parameters = 3
        case reactionParameterIndices = 4
        case stoichiometry = 5
        case reactions = 6
        case expressions = 7
        case actions = 8
        case rules = 9
        case monitors = 10
        case cohorts = 11
        case speciesIncidenceOffsets = 12
        case speciesIncidence = 13
        case runtimeContract = 14

        public var description: String {
            switch self {
            case .strings: "strings"
            case .species: "species"
            case .parameters: "parameters"
            case .reactionParameterIndices: "reaction-parameter-indices"
            case .stoichiometry: "stoichiometry"
            case .reactions: "reactions"
            case .expressions: "expressions"
            case .actions: "actions"
            case .rules: "rules"
            case .monitors: "monitors"
            case .cohorts: "cohorts"
            case .speciesIncidenceOffsets: "species-incidence-offsets"
            case .speciesIncidence: "species-incidence"
            case .runtimeContract: "runtime-contract"
            }
        }
    }

    public struct Header: Sendable, Codable, Equatable {
        public let major: UInt16
        public let minor: UInt16
        public let headerBytes: UInt32
        public let compilerABI: UInt32
        public let flags: UInt32
        public let fidelity: VivoFidelity
        public let sectionCount: UInt32
        public let totalBytes: UInt64
        public let sourceFingerprint: VivoFingerprint
        public let contentFingerprint: VivoFingerprint
    }

    public struct Section: Sendable, Codable, Equatable {
        public let kind: SectionKind
        public let flags: UInt32
        public let offset: UInt64
        public let size: UInt64
        public let stride: UInt32
        public let count: UInt32
        public let alignment: UInt32
        public let fingerprint: VivoFingerprint

        public var range: Range<Int> {
            Int(offset)..<Int(offset + size)
        }
    }

    public struct RuntimeContract: Sendable, Codable, Equatable {
        public let speciesCount: UInt32
        public let parameterCount: UInt32
        public let reactionCount: UInt32
        public let ruleCount: UInt32
        public let monitorCount: UInt32
        public let cohortCount: UInt32
        public let temporalStateCount: UInt32
        public let maximumExpressionStack: UInt32
        public let featureFlags: UInt32
        public let authoritativeScalarBytes: UInt32
        public let randomStreamVersion: UInt32
    }

    public struct SpeciesMetadata: Sendable, Codable, Equatable {
        public let identifier: String
        public let compartment: String
        public let unit: String
        public let flags: UInt32
        public let initialValue: Float
        public let minimum: Float
        public let maximum: Float

        public var isExternallyOwned: Bool { flags & (1 << 1) != 0 }
        public var isInput: Bool { flags & (1 << 2) != 0 }
        public var isState: Bool { flags & (1 << 3) != 0 }
        public var isOutput: Bool { flags & (1 << 4) != 0 }
        public var isCountValued: Bool { flags & (1 << 5) != 0 }
    }

    public struct ParameterMetadata: Sendable, Codable, Equatable {
        public let identifier: String
        public let unit: String
        public let evidenceSource: String
        public let flags: UInt32
        public let value: Double
        public let minimum: Double
        public let maximum: Double
        public let evidenceClass: UInt32
    }

    public let data: Data
    public let header: Header
    public let sections: [SectionKind: Section]
    public let runtimeContract: RuntimeContract

    public init(data: Data) throws {
        let reader = BinaryReader(data: data)
        guard data.count >= Self.fixedHeaderSize else {
            throw VivoProgramPackError.truncated("fixed header")
        }
        let magic = try reader.bytes(at: 0, count: 8)
        guard magic == [78, 86, 73, 86, 79, 80, 75, 0] else {
            throw VivoProgramPackError.invalidMagic
        }

        let major = try reader.u16(at: 8)
        let minor = try reader.u16(at: 10)
        guard major == Self.supportedMajor else {
            throw VivoProgramPackError.unsupportedVersion(major: major, minor: minor)
        }
        let headerBytes = try reader.u32(at: 12)
        let compilerABI = try reader.u32(at: 16)
        guard compilerABI == Self.compilerABI else {
            throw VivoProgramPackError.incompatibleCompilerABI(compilerABI)
        }
        let flags = try reader.u32(at: 20)
        let fidelityRaw = try reader.u32(at: 24)
        guard let fidelity = VivoFidelity(rawValue: fidelityRaw) else {
            throw VivoProgramPackError.invalidHeader("unknown fidelity value \(fidelityRaw)")
        }
        let sectionCount = try reader.u32(at: 28)
        let totalBytes = try reader.u64(at: 32)
        guard totalBytes == UInt64(data.count) else {
            throw VivoProgramPackError.invalidHeader("totalBytes \(totalBytes) does not match data length \(data.count)")
        }
        let sourceFingerprint = try VivoFingerprint(bytes: reader.bytes(at: 40, count: 32))
        let contentFingerprint = try VivoFingerprint(bytes: reader.bytes(at: 72, count: 32))

        let descriptorBytes = try Self.checkedMultiply(UInt64(sectionCount), UInt64(Self.sectionDescriptorSize), field: "descriptor table")
        let minimumHeader = try Self.checkedAdd(UInt64(Self.fixedHeaderSize), descriptorBytes, field: "header")
        guard UInt64(headerBytes) >= minimumHeader, UInt64(headerBytes) <= totalBytes else {
            throw VivoProgramPackError.invalidHeader("headerBytes is outside the descriptor table bounds")
        }

        var parsedSections: [SectionKind: Section] = [:]
        parsedSections.reserveCapacity(Int(sectionCount))
        var occupied: [Range<UInt64>] = []
        occupied.reserveCapacity(Int(sectionCount))

        for index in 0..<Int(sectionCount) {
            let base = Self.fixedHeaderSize + index * Self.sectionDescriptorSize
            let typeRaw = try reader.u32(at: base)
            guard let kind = SectionKind(rawValue: typeRaw) else {
                throw VivoProgramPackError.invalidSection("unknown type \(typeRaw) at descriptor \(index)")
            }
            guard parsedSections[kind] == nil else {
                throw VivoProgramPackError.invalidSection("duplicate section \(kind)")
            }
            let sectionFlags = try reader.u32(at: base + 4)
            let offset = try reader.u64(at: base + 8)
            let size = try reader.u64(at: base + 16)
            let stride = try reader.u32(at: base + 24)
            let count = try reader.u32(at: base + 28)
            let alignment = try reader.u32(at: base + 32)
            let fingerprint = try VivoFingerprint(bytes: reader.bytes(at: base + 40, count: 32))

            guard alignment > 0, offset % UInt64(alignment) == 0 else {
                throw VivoProgramPackError.invalidSection("\(kind) violates declared alignment \(alignment)")
            }
            let expectedSize = try Self.checkedMultiply(UInt64(stride), UInt64(count), field: "\(kind) size")
            guard expectedSize == size else {
                throw VivoProgramPackError.invalidSection("\(kind) size \(size) does not equal stride × count \(expectedSize)")
            }
            let end = try Self.checkedAdd(offset, size, field: "\(kind) range")
            guard offset >= UInt64(headerBytes), end <= totalBytes else {
                throw VivoProgramPackError.invalidSection("\(kind) lies outside the pack payload")
            }
            let range = offset..<end
            for existing in occupied where range.overlaps(existing) {
                throw VivoProgramPackError.invalidSection("\(kind) overlaps another section")
            }
            occupied.append(range)
            parsedSections[kind] = Section(
                kind: kind,
                flags: sectionFlags,
                offset: offset,
                size: size,
                stride: stride,
                count: count,
                alignment: alignment,
                fingerprint: fingerprint
            )
        }

        for required in [SectionKind.strings, .species, .parameters, .reactions, .runtimeContract] {
            guard parsedSections[required] != nil else {
                throw VivoProgramPackError.missingSection(required)
            }
        }

        let parsedHeader = Header(
            major: major,
            minor: minor,
            headerBytes: headerBytes,
            compilerABI: compilerABI,
            flags: flags,
            fidelity: fidelity,
            sectionCount: sectionCount,
            totalBytes: totalBytes,
            sourceFingerprint: sourceFingerprint,
            contentFingerprint: contentFingerprint
        )
        let contract = try Self.parseRuntimeContract(reader: reader, section: parsedSections[.runtimeContract]!)
        guard contract.speciesCount == parsedSections[.species]?.count else {
            throw VivoProgramPackError.invalidRuntimeContract("species count disagrees with the species section")
        }
        guard contract.parameterCount == parsedSections[.parameters]?.count else {
            throw VivoProgramPackError.invalidRuntimeContract("parameter count disagrees with the parameters section")
        }
        guard contract.reactionCount == parsedSections[.reactions]?.count else {
            throw VivoProgramPackError.invalidRuntimeContract("reaction count disagrees with the reactions section")
        }
        guard contract.maximumExpressionStack <= 256 else {
            throw VivoProgramPackError.invalidRuntimeContract("expression stack exceeds the Metal runtime limit")
        }
        guard contract.authoritativeScalarBytes == 4 else {
            throw VivoProgramPackError.invalidRuntimeContract("authoritative scalar width must be FP32")
        }

        self.data = data
        self.header = parsedHeader
        self.sections = parsedSections
        self.runtimeContract = contract
    }

    public func section(_ kind: SectionKind) throws -> Section {
        guard let section = sections[kind] else {
            throw VivoProgramPackError.missingSection(kind)
        }
        return section
    }

    public func sectionData(_ kind: SectionKind) throws -> Data {
        let descriptor = try section(kind)
        return data.subdata(in: descriptor.range)
    }

    public func speciesMetadata() throws -> [SpeciesMetadata] {
        let speciesSection = try section(.species)
        guard speciesSection.stride == 32 else {
            throw VivoProgramPackError.invalidSection("species stride must be 32 bytes")
        }
        let strings = try section(.strings)
        let reader = BinaryReader(data: data)
        var result: [SpeciesMetadata] = []
        result.reserveCapacity(Int(speciesSection.count))
        for index in 0..<Int(speciesSection.count) {
            let base = Int(speciesSection.offset) + index * Int(speciesSection.stride)
            result.append(SpeciesMetadata(
                identifier: try readString(offset: try reader.u32(at: base), strings: strings, reader: reader),
                compartment: try readString(offset: try reader.u32(at: base + 4), strings: strings, reader: reader),
                unit: try readString(offset: try reader.u32(at: base + 8), strings: strings, reader: reader),
                flags: try reader.u32(at: base + 12),
                initialValue: try reader.f32(at: base + 16),
                minimum: try reader.f32(at: base + 20),
                maximum: try reader.f32(at: base + 24)
            ))
        }
        return result
    }

    public func parameterMetadata() throws -> [ParameterMetadata] {
        let parameterSection = try section(.parameters)
        guard parameterSection.stride == 48 else {
            throw VivoProgramPackError.invalidSection("parameter stride must be 48 bytes")
        }
        let strings = try section(.strings)
        let reader = BinaryReader(data: data)
        var result: [ParameterMetadata] = []
        result.reserveCapacity(Int(parameterSection.count))
        for index in 0..<Int(parameterSection.count) {
            let base = Int(parameterSection.offset) + index * Int(parameterSection.stride)
            result.append(ParameterMetadata(
                identifier: try readString(offset: try reader.u32(at: base), strings: strings, reader: reader),
                unit: try readString(offset: try reader.u32(at: base + 4), strings: strings, reader: reader),
                evidenceSource: try readString(offset: try reader.u32(at: base + 8), strings: strings, reader: reader),
                flags: try reader.u32(at: base + 12),
                value: try reader.f64(at: base + 16),
                minimum: try reader.f64(at: base + 24),
                maximum: try reader.f64(at: base + 32),
                evidenceClass: try reader.u32(at: base + 40)
            ))
        }
        return result
    }

    public func initialState(laneCount: Int) throws -> [Float] {
        guard laneCount > 0 else {
            throw VivoProgramPackError.invalidNumericValue("laneCount must be positive")
        }
        let metadata = try speciesMetadata()
        let total = try Self.checkedElementCount(metadata.count, laneCount, field: "initial state")
        var state = [Float](repeating: 0, count: total)
        for (speciesIndex, species) in metadata.enumerated() {
            let start = speciesIndex * laneCount
            state.replaceSubrange(start..<(start + laneCount), with: repeatElement(species.initialValue, count: laneCount))
        }
        return state
    }

    public func parameterValues(environmentCount: Int) throws -> [Float] {
        guard environmentCount > 0 else {
            throw VivoProgramPackError.invalidNumericValue("environmentCount must be positive")
        }
        let metadata = try parameterMetadata()
        let total = try Self.checkedElementCount(metadata.count, environmentCount, field: "parameter values")
        var values = [Float](repeating: 0, count: total)
        for (parameterIndex, parameter) in metadata.enumerated() {
            guard parameter.value.isFinite,
                  abs(parameter.value) <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoProgramPackError.invalidNumericValue("parameter \(parameter.identifier) cannot be represented as FP32")
            }
            let value = Float(parameter.value)
            let start = parameterIndex * environmentCount
            values.replaceSubrange(start..<(start + environmentCount), with: repeatElement(value, count: environmentCount))
        }
        return values
    }

    public func speciesIndex(named identifier: String) throws -> UInt32? {
        let species = try speciesMetadata()
        guard let index = species.firstIndex(where: { $0.identifier == identifier }) else { return nil }
        return UInt32(index)
    }

    private static func parseRuntimeContract(reader: BinaryReader, section: Section) throws -> RuntimeContract {
        guard section.count == 1, section.stride == 80, section.size == 80 else {
            throw VivoProgramPackError.invalidRuntimeContract("runtime-contract section must contain one 80-byte record")
        }
        let base = Int(section.offset)
        return RuntimeContract(
            speciesCount: try reader.u32(at: base),
            parameterCount: try reader.u32(at: base + 4),
            reactionCount: try reader.u32(at: base + 8),
            ruleCount: try reader.u32(at: base + 12),
            monitorCount: try reader.u32(at: base + 16),
            cohortCount: try reader.u32(at: base + 20),
            temporalStateCount: try reader.u32(at: base + 24),
            maximumExpressionStack: try reader.u32(at: base + 28),
            featureFlags: try reader.u32(at: base + 32),
            authoritativeScalarBytes: try reader.u32(at: base + 36),
            randomStreamVersion: try reader.u32(at: base + 40)
        )
    }

    private func readString(offset: UInt32, strings: Section, reader: BinaryReader) throws -> String {
        guard UInt64(offset) < strings.size else {
            throw VivoProgramPackError.invalidSection("string offset \(offset) is outside the string table")
        }
        let start = Int(strings.offset) + Int(offset)
        let limit = Int(strings.offset + strings.size)
        var bytes: [UInt8] = []
        var cursor = start
        while cursor < limit {
            let byte = try reader.u8(at: cursor)
            if byte == 0 { return String(decoding: bytes, as: UTF8.self) }
            bytes.append(byte)
            cursor += 1
        }
        throw VivoProgramPackError.invalidSection("unterminated string at offset \(offset)")
    }

    private static func checkedAdd(_ lhs: UInt64, _ rhs: UInt64, field: String) throws -> UInt64 {
        let (result, overflow) = lhs.addingReportingOverflow(rhs)
        guard !overflow else { throw VivoProgramPackError.invalidHeader("\(field) overflows UInt64") }
        return result
    }

    private static func checkedMultiply(_ lhs: UInt64, _ rhs: UInt64, field: String) throws -> UInt64 {
        let (result, overflow) = lhs.multipliedReportingOverflow(by: rhs)
        guard !overflow else { throw VivoProgramPackError.invalidHeader("\(field) overflows UInt64") }
        return result
    }

    private static func checkedElementCount(_ lhs: Int, _ rhs: Int, field: String) throws -> Int {
        let (result, overflow) = lhs.multipliedReportingOverflow(by: rhs)
        guard !overflow else { throw VivoProgramPackError.invalidNumericValue("\(field) element count overflow") }
        return result
    }
}

private struct BinaryReader: Sendable {
    let data: Data

    func u8(at offset: Int) throws -> UInt8 {
        guard offset >= 0, offset < data.count else {
            throw VivoProgramPackError.truncated("byte at offset \(offset)")
        }
        return data[data.index(data.startIndex, offsetBy: offset)]
    }

    func bytes(at offset: Int, count: Int) throws -> [UInt8] {
        guard offset >= 0, count >= 0, offset <= data.count, count <= data.count - offset else {
            throw VivoProgramPackError.truncated("\(count) bytes at offset \(offset)")
        }
        let start = data.index(data.startIndex, offsetBy: offset)
        let end = data.index(start, offsetBy: count)
        return Array(data[start..<end])
    }

    func u16(at offset: Int) throws -> UInt16 {
        UInt16(try u8(at: offset)) |
        UInt16(try u8(at: offset + 1)) << 8
    }

    func u32(at offset: Int) throws -> UInt32 {
        UInt32(try u8(at: offset)) |
        UInt32(try u8(at: offset + 1)) << 8 |
        UInt32(try u8(at: offset + 2)) << 16 |
        UInt32(try u8(at: offset + 3)) << 24
    }

    func u64(at offset: Int) throws -> UInt64 {
        UInt64(try u32(at: offset)) | UInt64(try u32(at: offset + 4)) << 32
    }

    func f32(at offset: Int) throws -> Float {
        Float(bitPattern: try u32(at: offset))
    }

    func f64(at offset: Int) throws -> Double {
        Double(bitPattern: try u64(at: offset))
    }
}
