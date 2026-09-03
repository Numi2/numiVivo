import Foundation
import NumiVivoCore

public enum VivoFidelity: UInt32, Codable, Sendable, CaseIterable {
    case logic = 0
    case deterministic = 1
    case stochastic = 2
    case spatial = 3
    case tissue = 4

    public var label: String {
        switch self {
        case .logic: "F0"
        case .deterministic: "F1"
        case .stochastic: "F2"
        case .spatial: "F3"
        case .tissue: "F4"
        }
    }
}

public struct VivoCompilerLimits: Sendable, Codable, Equatable {
    public var maximumSpecies: UInt32
    public var maximumParameters: UInt32
    public var maximumReactions: UInt32
    public var maximumRules: UInt32
    public var maximumConstraints: UInt32
    public var maximumExpressionInstructions: UInt32
    public var maximumStoichiometryTerms: UInt32
    public var maximumTemporalStates: UInt32

    public init(
        maximumSpecies: UInt32 = 16_384,
        maximumParameters: UInt32 = 65_536,
        maximumReactions: UInt32 = 65_536,
        maximumRules: UInt32 = 16_384,
        maximumConstraints: UInt32 = 16_384,
        maximumExpressionInstructions: UInt32 = 1_048_576,
        maximumStoichiometryTerms: UInt32 = 1_048_576,
        maximumTemporalStates: UInt32 = 262_144
    ) {
        self.maximumSpecies = maximumSpecies
        self.maximumParameters = maximumParameters
        self.maximumReactions = maximumReactions
        self.maximumRules = maximumRules
        self.maximumConstraints = maximumConstraints
        self.maximumExpressionInstructions = maximumExpressionInstructions
        self.maximumStoichiometryTerms = maximumStoichiometryTerms
        self.maximumTemporalStates = maximumTemporalStates
    }
}

public struct VivoCompilerOptions: Sendable, Codable, Equatable {
    public var fidelity: VivoFidelity
    public var strictUnits: Bool
    public var strictSafety: Bool
    public var deterministicPack: Bool
    public var permitHypotheticalParameters: Bool
    public var requireTermination: Bool
    public var limits: VivoCompilerLimits

    public init(
        fidelity: VivoFidelity = .stochastic,
        strictUnits: Bool = true,
        strictSafety: Bool = true,
        deterministicPack: Bool = true,
        permitHypotheticalParameters: Bool = true,
        requireTermination: Bool = true,
        limits: VivoCompilerLimits = .init()
    ) {
        self.fidelity = fidelity
        self.strictUnits = strictUnits
        self.strictSafety = strictSafety
        self.deterministicPack = deterministicPack
        self.permitHypotheticalParameters = permitHypotheticalParameters
        self.requireTermination = requireTermination
        self.limits = limits
    }
}

public enum VivoCompilerStatus: String, Sendable, Codable {
    case accepted
    case invalidArgument
    case parseError
    case validationError
    case safetyRejected
    case compileError
    case invalidPack
    case resourceLimit
    case outOfMemory
    case internalError

    init(_ status: NVivoStatus) {
        switch status {
        case NVIVO_STATUS_OK: self = .accepted
        case NVIVO_STATUS_INVALID_ARGUMENT: self = .invalidArgument
        case NVIVO_STATUS_PARSE_ERROR: self = .parseError
        case NVIVO_STATUS_VALIDATION_ERROR: self = .validationError
        case NVIVO_STATUS_SAFETY_REJECTED: self = .safetyRejected
        case NVIVO_STATUS_COMPILE_ERROR: self = .compileError
        case NVIVO_STATUS_INVALID_PACK: self = .invalidPack
        case NVIVO_STATUS_RESOURCE_LIMIT: self = .resourceLimit
        case NVIVO_STATUS_OUT_OF_MEMORY: self = .outOfMemory
        default: self = .internalError
        }
    }
}

public struct VivoCompilerReport: Sendable, Codable, Equatable {
    public let status: VivoCompilerStatus
    public let json: String

    public var accepted: Bool { status == .accepted }

    public init(status: VivoCompilerStatus, json: String) {
        self.status = status
        self.json = json
    }
}

public struct VivoCompilation: Sendable {
    public let programPack: VivoProgramPack
    public let report: VivoCompilerReport

    public init(programPack: VivoProgramPack, report: VivoCompilerReport) {
        self.programPack = programPack
        self.report = report
    }
}

public struct VivoCompilerFailure: Error, Sendable, CustomStringConvertible {
    public let report: VivoCompilerReport

    public var description: String {
        "NumiVivo compiler returned \(report.status.rawValue): \(report.json)"
    }

    public init(report: VivoCompilerReport) {
        self.report = report
    }
}

public struct VivoProgramCompiler: Sendable {
    private static let strictUnitsFlag: UInt32 = 1 << 0
    private static let strictSafetyFlag: UInt32 = 1 << 1
    private static let deterministicPackFlag: UInt32 = 1 << 2
    private static let permitHypotheticalFlag: UInt32 = 1 << 3
    private static let requireTerminationFlag: UInt32 = 1 << 4

    public init() {}

    public func compile(
        source: Data,
        options: VivoCompilerOptions = .init()
    ) throws -> VivoCompilation {
        guard !source.isEmpty else {
            throw VivoCompilerFailure(report: .init(
                status: .invalidArgument,
                json: #"{"diagnostics":[{"severity":"fatal","code":"NVSW001","message":"VivoProgram source is empty."}]}"#
            ))
        }

        var rawOptions = makeRawOptions(options)
        var packBuffer = NVivoByteBuffer(data: nil, size: 0)
        var diagnosticsBuffer = NVivoByteBuffer(data: nil, size: 0)
        defer {
            nvivo_buffer_release(&packBuffer)
            nvivo_buffer_release(&diagnosticsBuffer)
        }

        let status: NVivoStatus = source.withUnsafeBytes { rawBytes in
            let bytes = rawBytes.bindMemory(to: UInt8.self)
            return nvivo_compile_program_json(
                bytes.baseAddress,
                bytes.count,
                &rawOptions,
                &packBuffer,
                &diagnosticsBuffer
            )
        }
        let report = VivoCompilerReport(
            status: VivoCompilerStatus(status),
            json: string(from: diagnosticsBuffer)
        )
        guard status == NVIVO_STATUS_OK else {
            throw VivoCompilerFailure(report: report)
        }
        guard let pointer = packBuffer.data, packBuffer.size > 0 else {
            throw VivoCompilerFailure(report: .init(
                status: .internalError,
                json: #"{"diagnostics":[{"severity":"fatal","code":"NVSW002","message":"Compiler returned no ProgramPack bytes."}]}"#
            ))
        }
        let data = Data(bytes: pointer, count: packBuffer.size)
        return VivoCompilation(programPack: try VivoProgramPack(data: data), report: report)
    }

    public func compile(
        source: String,
        options: VivoCompilerOptions = .init()
    ) throws -> VivoCompilation {
        guard let data = source.data(using: .utf8) else {
            throw VivoCompilerFailure(report: .init(
                status: .invalidArgument,
                json: #"{"diagnostics":[{"severity":"fatal","code":"NVSW003","message":"VivoProgram source is not valid UTF-8."}]}"#
            ))
        }
        return try compile(source: data, options: options)
    }

    public func compile(
        fileURL: URL,
        options: VivoCompilerOptions = .init()
    ) throws -> VivoCompilation {
        try compile(source: Data(contentsOf: fileURL, options: [.mappedIfSafe]), options: options)
    }

    public func validate(
        source: Data,
        options: VivoCompilerOptions = .init()
    ) -> VivoCompilerReport {
        guard !source.isEmpty else {
            return .init(
                status: .invalidArgument,
                json: #"{"diagnostics":[{"severity":"fatal","code":"NVSW004","message":"VivoProgram source is empty."}]}"#
            )
        }
        var rawOptions = makeRawOptions(options)
        var diagnosticsBuffer = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&diagnosticsBuffer) }
        let status: NVivoStatus = source.withUnsafeBytes { rawBytes in
            let bytes = rawBytes.bindMemory(to: UInt8.self)
            return nvivo_validate_program_json(
                bytes.baseAddress,
                bytes.count,
                &rawOptions,
                &diagnosticsBuffer
            )
        }
        return .init(status: VivoCompilerStatus(status), json: string(from: diagnosticsBuffer))
    }

    public func validate(
        source: String,
        options: VivoCompilerOptions = .init()
    ) -> VivoCompilerReport {
        guard let data = source.data(using: .utf8) else {
            return .init(
                status: .invalidArgument,
                json: #"{"diagnostics":[{"severity":"fatal","code":"NVSW005","message":"VivoProgram source is not valid UTF-8."}]}"#
            )
        }
        return validate(source: data, options: options)
    }

    public func inspect(_ pack: VivoProgramPack, verifySectionHashes: Bool = true) -> VivoCompilerReport {
        var summary = NVivoPackSummary()
        var inspectionBuffer = NVivoByteBuffer(data: nil, size: 0)
        defer { nvivo_buffer_release(&inspectionBuffer) }
        let status: NVivoStatus = pack.data.withUnsafeBytes { rawBytes in
            let bytes = rawBytes.bindMemory(to: UInt8.self)
            return nvivo_inspect_program_pack(
                bytes.baseAddress,
                bytes.count,
                verifySectionHashes ? 1 : 0,
                &summary,
                &inspectionBuffer
            )
        }
        return .init(status: VivoCompilerStatus(status), json: string(from: inspectionBuffer))
    }

    private func makeRawOptions(_ options: VivoCompilerOptions) -> NVivoCompileOptions {
        var result = NVivoCompileOptions()
        nvivo_default_compile_options(&result)
        result.requested_fidelity = options.fidelity.rawValue
        result.flags = 0
        if options.strictUnits { result.flags |= Self.strictUnitsFlag }
        if options.strictSafety { result.flags |= Self.strictSafetyFlag }
        if options.deterministicPack { result.flags |= Self.deterministicPackFlag }
        if options.permitHypotheticalParameters { result.flags |= Self.permitHypotheticalFlag }
        if options.requireTermination { result.flags |= Self.requireTerminationFlag }
        result.limits.maximum_species = options.limits.maximumSpecies
        result.limits.maximum_parameters = options.limits.maximumParameters
        result.limits.maximum_reactions = options.limits.maximumReactions
        result.limits.maximum_rules = options.limits.maximumRules
        result.limits.maximum_constraints = options.limits.maximumConstraints
        result.limits.maximum_expression_instructions = options.limits.maximumExpressionInstructions
        result.limits.maximum_stoichiometry_terms = options.limits.maximumStoichiometryTerms
        result.limits.maximum_temporal_states = options.limits.maximumTemporalStates
        return result
    }

    private func string(from buffer: NVivoByteBuffer) -> String {
        guard let pointer = buffer.data, buffer.size > 0 else { return "{}" }
        return String(decoding: UnsafeBufferPointer(start: pointer, count: buffer.size), as: UTF8.self)
    }
}
