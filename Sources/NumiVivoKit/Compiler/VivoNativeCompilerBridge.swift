import Foundation
import NumiVivoCore

public enum VivoNativeStatus: Int, Codable, Sendable {
    case ok = 0
    case invalidArgument = 1
    case parseError = 2
    case validationError = 3
    case safetyRejected = 4
    case compileError = 5
    case invalidPack = 6
    case resourceLimit = 7
    case outOfMemory = 8
    case internalError = 255
    case unknown = -1

    public var succeeded: Bool { self == .ok }
}

public struct VivoNativeInvocation: Sendable {
    public let status: VivoNativeStatus
    public let primary: Data
    public let diagnostics: Data

    public init(status: VivoNativeStatus, primary: Data, diagnostics: Data) {
        self.status = status
        self.primary = primary
        self.diagnostics = diagnostics
    }

    public var succeeded: Bool { status.succeeded }

    public func diagnosticsString() -> String {
        String(decoding: diagnostics, as: UTF8.self)
    }
}

public struct VivoNativePackInspection: Sendable {
    public let invocation: VivoNativeInvocation
    public let summary: NVivoPackSummary

    public init(invocation: VivoNativeInvocation, summary: NVivoPackSummary) {
        self.invocation = invocation
        self.summary = summary
    }
}

public enum VivoNativeBridgeError: Error, LocalizedError, Sendable {
    case invalidFidelity(UInt32)
    case invalidAnalysisLimits

    public var errorDescription: String? {
        switch self {
        case .invalidFidelity(let value):
            return "native compiler fidelity \(value) is invalid"
        case .invalidAnalysisLimits:
            return "native reaction-network analysis limits are invalid"
        }
    }
}

public enum VivoNativeCompilerBridge {
    public static func validateProgram(
        _ source: Data,
        fidelity: VivoFidelity = .stochastic,
        strictUnits: Bool = true,
        strictSafety: Bool = true,
        requireTermination: Bool = true
    ) throws -> VivoNativeInvocation {
        var options = try compileOptions(
            fidelity: fidelity,
            strictUnits: strictUnits,
            strictSafety: strictSafety,
            requireTermination: requireTermination
        )
        var diagnostics = NVivoByteBuffer(data: nil, size: 0)
        let status = source.withUnsafeBytes { raw -> NVivoStatus in
            nvivo_validate_program_json(
                raw.bindMemory(to: UInt8.self).baseAddress,
                raw.count,
                &options,
                &diagnostics
            )
        }
        return .init(
            status: nativeStatus(status),
            primary: Data(),
            diagnostics: takeBuffer(&diagnostics)
        )
    }

    public static func compileProgram(
        _ source: Data,
        fidelity: VivoFidelity = .stochastic,
        strictUnits: Bool = true,
        strictSafety: Bool = true,
        requireTermination: Bool = true
    ) throws -> VivoNativeInvocation {
        var options = try compileOptions(
            fidelity: fidelity,
            strictUnits: strictUnits,
            strictSafety: strictSafety,
            requireTermination: requireTermination
        )
        var pack = NVivoByteBuffer(data: nil, size: 0)
        var diagnostics = NVivoByteBuffer(data: nil, size: 0)
        let status = source.withUnsafeBytes { raw -> NVivoStatus in
            nvivo_compile_program_json(
                raw.bindMemory(to: UInt8.self).baseAddress,
                raw.count,
                &options,
                &pack,
                &diagnostics
            )
        }
        return .init(
            status: nativeStatus(status),
            primary: takeBuffer(&pack),
            diagnostics: takeBuffer(&diagnostics)
        )
    }

    public static func inspectProgramPack(
        _ pack: Data,
        verifySectionHashes: Bool = true
    ) -> VivoNativePackInspection {
        var summary = NVivoPackSummary()
        var report = NVivoByteBuffer(data: nil, size: 0)
        let status = pack.withUnsafeBytes { raw -> NVivoStatus in
            nvivo_inspect_program_pack(
                raw.bindMemory(to: UInt8.self).baseAddress,
                raw.count,
                verifySectionHashes ? 1 : 0,
                &summary,
                &report
            )
        }
        let invocation = VivoNativeInvocation(
            status: nativeStatus(status),
            primary: takeBuffer(&report),
            diagnostics: Data()
        )
        return .init(invocation: invocation, summary: summary)
    }

    public static func analyzeProgram(
        _ source: Data,
        fidelity: VivoFidelity = .stochastic,
        maximumDenseSpecies: UInt32 = 2_048,
        maximumDenseReactions: UInt32 = 2_048,
        maximumDenseElements: UInt64 = 16_777_216,
        maximumConservationLaws: UInt32 = 256,
        rankTolerance: Double = 1e-10,
        conservationTolerance: Double = 1e-9
    ) throws -> VivoNativeInvocation {
        guard maximumDenseSpecies > 0,
              maximumDenseReactions > 0,
              maximumDenseElements > 0,
              maximumConservationLaws > 0,
              rankTolerance.isFinite,
              rankTolerance > 0,
              conservationTolerance.isFinite,
              conservationTolerance > 0 else {
            throw VivoNativeBridgeError.invalidAnalysisLimits
        }
        var compiler = try compileOptions(
            fidelity: fidelity,
            strictUnits: true,
            strictSafety: false,
            requireTermination: false
        )
        var analysis = NVivoNetworkAnalysisOptions()
        nvivo_default_network_analysis_options(&analysis)
        analysis.maximum_dense_species = maximumDenseSpecies
        analysis.maximum_dense_reactions = maximumDenseReactions
        analysis.maximum_dense_elements = maximumDenseElements
        analysis.maximum_conservation_laws = maximumConservationLaws
        analysis.rank_tolerance = rankTolerance
        analysis.conservation_tolerance = conservationTolerance

        var report = NVivoByteBuffer(data: nil, size: 0)
        var diagnostics = NVivoByteBuffer(data: nil, size: 0)
        let status = source.withUnsafeBytes { raw -> NVivoStatus in
            nvivo_analyze_program_json(
                raw.bindMemory(to: UInt8.self).baseAddress,
                raw.count,
                &compiler,
                &analysis,
                &report,
                &diagnostics
            )
        }
        return .init(
            status: nativeStatus(status),
            primary: takeBuffer(&report),
            diagnostics: takeBuffer(&diagnostics)
        )
    }

    public static func synthesizeMechanisms(
        problem: Data,
        library: Data,
        maximumSolutions: UInt32? = nil,
        maximumVisitedNodes: UInt64? = nil,
        requireIndependentShutdown: Bool? = nil,
        requireMonitor: Bool? = nil,
        requireContextInsulation: Bool? = nil,
        requireResourceBuffering: Bool? = nil,
        requireDistinctOrthogonality: Bool? = nil
    ) -> VivoNativeInvocation {
        var options = NVivoMechanismSynthesisOptions()
        nvivo_default_mechanism_synthesis_options(&options)
        if let maximumSolutions { options.maximum_solutions = maximumSolutions }
        if let maximumVisitedNodes { options.maximum_visited_nodes = maximumVisitedNodes }

        func setFlag(_ flag: UInt32, value: Bool?) {
            guard let value else { return }
            if value { options.flags |= flag }
            else { options.flags &= ~flag }
        }
        setFlag(UInt32(NVIVO_SYNTHESIS_REQUIRE_INDEPENDENT_SHUTDOWN.rawValue), value: requireIndependentShutdown)
        setFlag(UInt32(NVIVO_SYNTHESIS_REQUIRE_MONITOR.rawValue), value: requireMonitor)
        setFlag(UInt32(NVIVO_SYNTHESIS_REQUIRE_CONTEXT_INSULATION.rawValue), value: requireContextInsulation)
        setFlag(UInt32(NVIVO_SYNTHESIS_REQUIRE_RESOURCE_BUFFERING.rawValue), value: requireResourceBuffering)
        setFlag(UInt32(NVIVO_SYNTHESIS_REQUIRE_DISTINCT_ORTHOGONALITY.rawValue), value: requireDistinctOrthogonality)

        var result = NVivoByteBuffer(data: nil, size: 0)
        var diagnostics = NVivoByteBuffer(data: nil, size: 0)
        let status = problem.withUnsafeBytes { problemRaw -> NVivoStatus in
            library.withUnsafeBytes { libraryRaw -> NVivoStatus in
                nvivo_synthesize_mechanisms_json(
                    problemRaw.bindMemory(to: UInt8.self).baseAddress,
                    problemRaw.count,
                    libraryRaw.bindMemory(to: UInt8.self).baseAddress,
                    libraryRaw.count,
                    &options,
                    &result,
                    &diagnostics
                )
            }
        }
        return .init(
            status: nativeStatus(status),
            primary: takeBuffer(&result),
            diagnostics: takeBuffer(&diagnostics)
        )
    }

    private static func compileOptions(
        fidelity: VivoFidelity,
        strictUnits: Bool,
        strictSafety: Bool,
        requireTermination: Bool
    ) throws -> NVivoCompileOptions {
        var options = NVivoCompileOptions()
        nvivo_default_compile_options(&options)
        let fidelityValue = fidelity.rawValue
        guard fidelityValue <= UInt32(NVIVO_FIDELITY_F4_TISSUE.rawValue) else {
            throw VivoNativeBridgeError.invalidFidelity(fidelityValue)
        }
        options.requested_fidelity = fidelityValue
        setCompileFlag(
            &options,
            flag: UInt32(NVIVO_COMPILE_STRICT_UNITS.rawValue),
            enabled: strictUnits
        )
        setCompileFlag(
            &options,
            flag: UInt32(NVIVO_COMPILE_STRICT_SAFETY.rawValue),
            enabled: strictSafety
        )
        setCompileFlag(
            &options,
            flag: UInt32(NVIVO_COMPILE_REQUIRE_TERMINATION.rawValue),
            enabled: requireTermination
        )
        return options
    }

    private static func setCompileFlag(
        _ options: inout NVivoCompileOptions,
        flag: UInt32,
        enabled: Bool
    ) {
        if enabled { options.flags |= flag }
        else { options.flags &= ~flag }
    }

    private static func nativeStatus(_ status: NVivoStatus) -> VivoNativeStatus {
        VivoNativeStatus(rawValue: Int(status.rawValue)) ?? .unknown
    }

    private static func takeBuffer(_ buffer: inout NVivoByteBuffer) -> Data {
        defer { nvivo_buffer_release(&buffer) }
        guard let pointer = buffer.data, buffer.size > 0 else { return Data() }
        return Data(bytes: pointer, count: buffer.size)
    }
}
