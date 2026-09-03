import Foundation
import NumiVivoCore

public struct VivoCompilerConfiguration: Sendable, Hashable {
    public var fidelity: VivoFidelity
    public var strictUnits: Bool
    public var strictSafety: Bool
    public var deterministicPack: Bool
    public var permitHypotheticalParameters: Bool
    public var requireTermination: Bool

    public init(
        fidelity: VivoFidelity = .f2Stochastic,
        strictUnits: Bool = true,
        strictSafety: Bool = true,
        deterministicPack: Bool = true,
        permitHypotheticalParameters: Bool = true,
        requireTermination: Bool = true
    ) {
        self.fidelity = fidelity
        self.strictUnits = strictUnits
        self.strictSafety = strictSafety
        self.deterministicPack = deterministicPack
        self.permitHypotheticalParameters = permitHypotheticalParameters
        self.requireTermination = requireTermination
    }
}

public struct VivoCompilerOutput: Sendable {
    public var programPack: VivoProgramPack
    public var reportJSON: Data

    public var reportText: String {
        String(decoding: reportJSON, as: UTF8.self)
    }
}

public enum VivoCompiler {
    public static func compile(
        json source: Data,
        configuration: VivoCompilerConfiguration = .init()
    ) throws -> VivoCompilerOutput {
        var options = makeOptions(configuration)
        var pack = NVivoByteBuffer(data: nil, size: 0)
        var report = NVivoByteBuffer(data: nil, size: 0)
        defer {
            nvivo_buffer_release(&pack)
            nvivo_buffer_release(&report)
        }

        let status = source.withUnsafeBytes { rawBuffer in
            nvivo_compile_program_json(
                rawBuffer.bindMemory(to: UInt8.self).baseAddress,
                rawBuffer.count,
                &options,
                &pack,
                &report
            )
        }
        let reportData = copy(buffer: report)
        guard status.rawValue == 0 else {
            throw VivoRuntimeError.compiler(
                status: Int32(truncatingIfNeeded: status.rawValue),
                report: String(decoding: reportData, as: UTF8.self)
            )
        }
        return VivoCompilerOutput(
            programPack: try VivoProgramPack(data: copy(buffer: pack)),
            reportJSON: reportData
        )
    }

    public static func compile(
        file sourceURL: URL,
        configuration: VivoCompilerConfiguration = .init()
    ) throws -> VivoCompilerOutput {
        try compile(json: Data(contentsOf: sourceURL), configuration: configuration)
    }

    @discardableResult
    public static func validate(
        json source: Data,
        configuration: VivoCompilerConfiguration = .init()
    ) throws -> Data {
        var options = makeOptions(configuration)
        var report = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&report) }

        let status = source.withUnsafeBytes { rawBuffer in
            nvivo_validate_program_json(
                rawBuffer.bindMemory(to: UInt8.self).baseAddress,
                rawBuffer.count,
                &options,
                &report
            )
        }
        let reportData = copy(buffer: report)
        guard status.rawValue == 0 else {
            throw VivoRuntimeError.compiler(
                status: Int32(truncatingIfNeeded: status.rawValue),
                report: String(decoding: reportData, as: UTF8.self)
            )
        }
        return reportData
    }

    public static func validate(
        file sourceURL: URL,
        configuration: VivoCompilerConfiguration = .init()
    ) throws -> Data {
        try validate(json: Data(contentsOf: sourceURL), configuration: configuration)
    }

    private static func makeOptions(_ configuration: VivoCompilerConfiguration) -> NVivoCompileOptions {
        var options = NVivoCompileOptions()
        nvivo_default_compile_options(&options)
        options.requested_fidelity = configuration.fidelity.rawValue
        var flags: UInt32 = 0
        if configuration.strictUnits { flags |= 1 << 0 }
        if configuration.strictSafety { flags |= 1 << 1 }
        if configuration.deterministicPack { flags |= 1 << 2 }
        if configuration.permitHypotheticalParameters { flags |= 1 << 3 }
        if configuration.requireTermination { flags |= 1 << 4 }
        options.flags = flags
        return options
    }

    private static func copy(buffer: NVivoByteBuffer) -> Data {
        guard let pointer = buffer.data, buffer.size > 0 else { return Data() }
        return Data(bytes: pointer, count: buffer.size)
    }
}

public struct VivoProgramPack: Sendable, Hashable {
    public enum SectionType: UInt32, CaseIterable, Sendable {
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
    }

    public struct Header: Sendable, Hashable {
        public var major: UInt16
        public var minor: UInt16
        public var headerBytes: UInt32
        public var compilerABI: UInt32
        public var flags: UInt32
        public var fidelity: VivoFidelity
        public var sectionCount: UInt32
        public var totalBytes: UInt64
        public var sourceFingerprint: String
        public var contentFingerprint: String
    }

    public struct Section: Sendable, Hashable {
        public var type: SectionType
        public var flags: UInt32
        public var offset: Int
        public var size: Int
        public var stride: Int
        public var count: Int
        public var alignment: Int
        public var fingerprint: String
    }

    public struct RuntimeContract: Sendable, Hashable {
        public var speciesCount: Int
        public var parameterCount: Int
        public var reactionCount: Int
        public var ruleCount: Int
        public var monitorCount: Int
        public var cohortCount: Int
        public var temporalStateCount: Int
        public var maximumExpressionStack: Int
        public var featureFlags: UInt32
        public var authoritativeScalarBytes: Int
        public var randomStreamVersion: Int
        public var manifestStringOffset: Int
    }

    public struct Species: Sendable, Hashable {
        public var index: UInt32
        public var id: String
        public var compartment: String
        public var unit: String
        public var flags: UInt32
        public var initialValue: Float
        public var minimum: Float
        public var maximum: Float

        public var isConserved: Bool { flags & (1 << 0) != 0 }
        public var isExternallyOwned: Bool { flags & (1 << 1) != 0 }
        public var isInput: Bool { flags & (1 << 2) != 0 }
        public var isState: Bool { flags & (1 << 3) != 0 }
        public var isOutput: Bool { flags & (1 << 4) != 0 }
        public var isCountValued: Bool { flags & (1 << 5) != 0 }
    }

    public struct Reaction: Sendable, Hashable {
        public var index: UInt32
        public var id: String
        public var flags: UInt32
        public var delaySeconds: Float
        public var characteristicRate: Float
        public var cohortIndex: UInt32
    }

    public struct Cohort: Sendable, Hashable {
        public var index: UInt32
        public var reactionOffset: UInt32
        public var reactionCount: UInt32
        public var rateLaw: UInt32
        public var flags: UInt32
        public var maximumStableStep: Float
        public var stiffnessEstimate: Float
        public var preferredThreads: UInt32
    }

    public let data: Data
    public let header: Header
    public let sections: [SectionType: Section]
    public let runtimeContract: RuntimeContract
    public let species: [Species]
    public let reactions: [Reaction]
    public let cohorts: [Cohort]
    public let manifestJSON: Data?

    public init(data: Data) throws {
        guard data.count >= 128 else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack is smaller than its fixed header.")
        }

        var summary = NVivoPackSummary()
        var report = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&report) }
        let status = data.withUnsafeBytes { rawBuffer in
            nvivo_inspect_program_pack(
                rawBuffer.bindMemory(to: UInt8.self).baseAddress,
                rawBuffer.count,
                1,
                &summary,
                &report
            )
        }
        guard status.rawValue == 0 else {
            let reportData: Data
            if let pointer = report.data, report.size > 0 {
                reportData = Data(bytes: pointer, count: report.size)
            } else {
                reportData = Data()
            }
            throw VivoRuntimeError.invalidProgramPack(String(decoding: reportData, as: UTF8.self))
        }

        let major = data.nvivoUInt16(at: 8)
        let minor = data.nvivoUInt16(at: 10)
        let headerBytes = data.nvivoUInt32(at: 12)
        let compilerABI = data.nvivoUInt32(at: 16)
        let flags = data.nvivoUInt32(at: 20)
        let fidelityRaw = data.nvivoUInt32(at: 24)
        let sectionCount = data.nvivoUInt32(at: 28)
        let totalBytes = data.nvivoUInt64(at: 32)
        guard let fidelity = VivoFidelity(rawValue: fidelityRaw) else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack fidelity value \(fidelityRaw) is invalid.")
        }
        guard totalBytes == UInt64(data.count) else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack byte count does not match totalBytes.")
        }
        let expectedHeader = 128 + Int(sectionCount) * 72
        guard Int(headerBytes) == expectedHeader, expectedHeader <= data.count else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack descriptor table has an invalid size.")
        }

        self.data = data
        self.header = Header(
            major: major,
            minor: minor,
            headerBytes: headerBytes,
            compilerABI: compilerABI,
            flags: flags,
            fidelity: fidelity,
            sectionCount: sectionCount,
            totalBytes: totalBytes,
            sourceFingerprint: data.nvivoHex(range: 40..<72),
            contentFingerprint: data.nvivoHex(range: 72..<104)
        )

        var decodedSections: [SectionType: Section] = [:]
        decodedSections.reserveCapacity(Int(sectionCount))
        for index in 0..<Int(sectionCount) {
            let base = 128 + index * 72
            let typeRaw = data.nvivoUInt32(at: base)
            guard let type = SectionType(rawValue: typeRaw) else { continue }
            let offsetValue = data.nvivoUInt64(at: base + 8)
            let sizeValue = data.nvivoUInt64(at: base + 16)
            guard offsetValue <= UInt64(Int.max), sizeValue <= UInt64(Int.max) else {
                throw VivoRuntimeError.invalidProgramPack("ProgramPack section range exceeds the host address space.")
            }
            let offset = Int(offsetValue)
            let size = Int(sizeValue)
            guard offset >= 0, size >= 0, offset <= data.count, size <= data.count - offset else {
                throw VivoRuntimeError.invalidProgramPack("ProgramPack section \(type) is out of bounds.")
            }
            decodedSections[type] = Section(
                type: type,
                flags: data.nvivoUInt32(at: base + 4),
                offset: offset,
                size: size,
                stride: Int(data.nvivoUInt32(at: base + 24)),
                count: Int(data.nvivoUInt32(at: base + 28)),
                alignment: Int(data.nvivoUInt32(at: base + 32)),
                fingerprint: data.nvivoHex(range: (base + 40)..<(base + 72))
            )
        }
        self.sections = decodedSections

        guard let contractSection = decodedSections[.runtimeContract], contractSection.size >= 80 else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack is missing a complete runtime contract.")
        }
        let contractBase = contractSection.offset
        self.runtimeContract = RuntimeContract(
            speciesCount: Int(data.nvivoUInt32(at: contractBase)),
            parameterCount: Int(data.nvivoUInt32(at: contractBase + 4)),
            reactionCount: Int(data.nvivoUInt32(at: contractBase + 8)),
            ruleCount: Int(data.nvivoUInt32(at: contractBase + 12)),
            monitorCount: Int(data.nvivoUInt32(at: contractBase + 16)),
            cohortCount: Int(data.nvivoUInt32(at: contractBase + 20)),
            temporalStateCount: Int(data.nvivoUInt32(at: contractBase + 24)),
            maximumExpressionStack: Int(data.nvivoUInt32(at: contractBase + 28)),
            featureFlags: data.nvivoUInt32(at: contractBase + 32),
            authoritativeScalarBytes: Int(data.nvivoUInt32(at: contractBase + 36)),
            randomStreamVersion: Int(data.nvivoUInt32(at: contractBase + 40)),
            manifestStringOffset: Int(data.nvivoUInt64(at: contractBase + 48))
        )

        guard runtimeContract.maximumExpressionStack <= 64 else {
            throw VivoRuntimeError.invalidProgramPack(
                "Program requires an expression stack of \(runtimeContract.maximumExpressionStack), exceeding the Metal ABI limit of 64."
            )
        }
        guard runtimeContract.authoritativeScalarBytes == MemoryLayout<Float>.size else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack authoritative scalar width is unsupported.")
        }

        let stringsSection = try Self.requireSection(.strings, in: decodedSections)
        func string(at relativeOffset: UInt32) throws -> String {
            let start = stringsSection.offset + Int(relativeOffset)
            guard start >= stringsSection.offset, start < stringsSection.offset + stringsSection.size else {
                throw VivoRuntimeError.invalidProgramPack("String-table offset is out of bounds.")
            }
            let endLimit = stringsSection.offset + stringsSection.size
            var end = start
            while end < endLimit, data[end] != 0 { end += 1 }
            guard end < endLimit else {
                throw VivoRuntimeError.invalidProgramPack("String-table entry is not NUL terminated.")
            }
            return String(decoding: data[start..<end], as: UTF8.self)
        }

        let speciesSection = try Self.requireSection(.species, in: decodedSections)
        guard speciesSection.stride == 32, speciesSection.count == runtimeContract.speciesCount else {
            throw VivoRuntimeError.invalidProgramPack("Species table does not match the runtime contract.")
        }
        var decodedSpecies: [Species] = []
        decodedSpecies.reserveCapacity(speciesSection.count)
        for index in 0..<speciesSection.count {
            let base = speciesSection.offset + index * speciesSection.stride
            decodedSpecies.append(
                Species(
                    index: UInt32(index),
                    id: try string(at: data.nvivoUInt32(at: base)),
                    compartment: try string(at: data.nvivoUInt32(at: base + 4)),
                    unit: try string(at: data.nvivoUInt32(at: base + 8)),
                    flags: data.nvivoUInt32(at: base + 12),
                    initialValue: data.nvivoFloat32(at: base + 16),
                    minimum: data.nvivoFloat32(at: base + 20),
                    maximum: data.nvivoFloat32(at: base + 24)
                )
            )
        }
        self.species = decodedSpecies

        let reactionsSection = try Self.requireSection(.reactions, in: decodedSections)
        guard reactionsSection.stride == 64, reactionsSection.count == runtimeContract.reactionCount else {
            throw VivoRuntimeError.invalidProgramPack("Reaction table does not match the runtime contract.")
        }
        var decodedReactions: [Reaction] = []
        decodedReactions.reserveCapacity(reactionsSection.count)
        for index in 0..<reactionsSection.count {
            let base = reactionsSection.offset + index * reactionsSection.stride
            decodedReactions.append(
                Reaction(
                    index: UInt32(index),
                    id: try string(at: data.nvivoUInt32(at: base)),
                    flags: data.nvivoUInt32(at: base + 44),
                    delaySeconds: data.nvivoFloat32(at: base + 48),
                    characteristicRate: data.nvivoFloat32(at: base + 52),
                    cohortIndex: data.nvivoUInt32(at: base + 56)
                )
            )
        }
        self.reactions = decodedReactions

        let cohortsSection = try Self.requireSection(.cohorts, in: decodedSections)
        guard cohortsSection.stride == 32, cohortsSection.count == runtimeContract.cohortCount else {
            throw VivoRuntimeError.invalidProgramPack("Cohort table does not match the runtime contract.")
        }
        var decodedCohorts: [Cohort] = []
        decodedCohorts.reserveCapacity(cohortsSection.count)
        for index in 0..<cohortsSection.count {
            let base = cohortsSection.offset + index * cohortsSection.stride
            decodedCohorts.append(
                Cohort(
                    index: UInt32(index),
                    reactionOffset: data.nvivoUInt32(at: base),
                    reactionCount: data.nvivoUInt32(at: base + 4),
                    rateLaw: data.nvivoUInt32(at: base + 8),
                    flags: data.nvivoUInt32(at: base + 12),
                    maximumStableStep: data.nvivoFloat32(at: base + 16),
                    stiffnessEstimate: data.nvivoFloat32(at: base + 20),
                    preferredThreads: data.nvivoUInt32(at: base + 24)
                )
            )
        }
        self.cohorts = decodedCohorts

        if runtimeContract.manifestStringOffset > 0 {
            self.manifestJSON = try string(at: UInt32(runtimeContract.manifestStringOffset)).data(using: .utf8)
        } else {
            self.manifestJSON = nil
        }
    }

    public init(contentsOf url: URL) throws {
        try self.init(data: Data(contentsOf: url, options: [.mappedIfSafe]))
    }

    public static func == (left: VivoProgramPack, right: VivoProgramPack) -> Bool {
        left.header.contentFingerprint == right.header.contentFingerprint
    }

    public func hash(into hasher: inout Hasher) {
        hasher.combine(header.contentFingerprint)
    }

    public func data(for section: SectionType) throws -> Data {
        let descriptor = try Self.requireSection(section, in: sections)
        return data.subdata(in: descriptor.offset..<(descriptor.offset + descriptor.size))
    }

    public func withUnsafeSectionBytes<Result>(
        _ section: SectionType,
        _ body: (UnsafeRawBufferPointer) throws -> Result
    ) throws -> Result {
        let descriptor = try Self.requireSection(section, in: sections)
        return try data.withUnsafeBytes { rawBuffer in
            guard let baseAddress = rawBuffer.baseAddress else {
                return try body(UnsafeRawBufferPointer(start: nil, count: 0))
            }
            return try body(
                UnsafeRawBufferPointer(
                    start: baseAddress.advanced(by: descriptor.offset),
                    count: descriptor.size
                )
            )
        }
    }

    public func materializeGPUParameters() throws -> Data {
        var output = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&output) }
        let status = data.withUnsafeBytes { rawBuffer in
            nvivo_materialize_gpu_parameters(
                rawBuffer.bindMemory(to: UInt8.self).baseAddress,
                rawBuffer.count,
                &output
            )
        }
        guard status.rawValue == 0 else {
            throw VivoRuntimeError.invalidProgramPack("Unable to materialize the Metal parameter table.")
        }
        guard let pointer = output.data, output.size > 0 else { return Data() }
        return Data(bytes: pointer, count: output.size)
    }

    public var speciesByID: [String: Species] {
        Dictionary(uniqueKeysWithValues: species.map { ($0.id, $0) })
    }

    public var maximumReactionDelay: Float {
        reactions.lazy.map(\.delaySeconds).max() ?? 0
    }

    private static func requireSection(
        _ type: SectionType,
        in sections: [SectionType: Section]
    ) throws -> Section {
        guard let section = sections[type] else {
            throw VivoRuntimeError.invalidProgramPack("ProgramPack is missing section \(type).")
        }
        return section
    }
}

private extension Data {
    func nvivoUInt16(at offset: Int) -> UInt16 {
        withUnsafeBytes { rawBuffer in
            UInt16(littleEndian: rawBuffer.loadUnaligned(fromByteOffset: offset, as: UInt16.self))
        }
    }

    func nvivoUInt32(at offset: Int) -> UInt32 {
        withUnsafeBytes { rawBuffer in
            UInt32(littleEndian: rawBuffer.loadUnaligned(fromByteOffset: offset, as: UInt32.self))
        }
    }

    func nvivoUInt64(at offset: Int) -> UInt64 {
        withUnsafeBytes { rawBuffer in
            UInt64(littleEndian: rawBuffer.loadUnaligned(fromByteOffset: offset, as: UInt64.self))
        }
    }

    func nvivoFloat32(at offset: Int) -> Float {
        Float(bitPattern: nvivoUInt32(at: offset))
    }

    func nvivoHex(range: Range<Int>) -> String {
        let alphabet = Array("0123456789abcdef".utf8)
        var result = [UInt8]()
        result.reserveCapacity(range.count * 2)
        for byte in self[range] {
            result.append(alphabet[Int(byte >> 4)])
            result.append(alphabet[Int(byte & 0x0f)])
        }
        return String(decoding: result, as: UTF8.self)
    }
}
