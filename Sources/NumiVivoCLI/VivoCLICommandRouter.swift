import Foundation
import NumiVivoKit

struct VivoCLICommandRouter {
    func run(arguments: [String]) -> Int32 {
        do {
            let parsed = try VivoCLIArguments(arguments)
            guard let command = parsed.command else {
                writeStandardOutput(Self.help)
                return VivoCLIExit.success
            }

            switch command {
            case "help", "--help", "-h":
                writeStandardOutput(Self.help)
                return VivoCLIExit.success
            case "version", "--version":
                writeStandardOutput("NumiVivo \(NumiVivo.version) (ProgramPack v\(NumiVivo.programPackVersion))\n")
                return VivoCLIExit.success
            case "validate":
                return try validate(parsed)
            case "compile":
                return try compile(parsed)
            case "inspect":
                return try inspect(parsed)
            case "analyze":
                return try analyze(parsed)
            case "synthesize-mechanisms":
                return try synthesizeMechanisms(parsed)
            case "compile-campaign":
                return try compileCampaign(parsed)
            case "verify-campaign":
                return try verifyCampaign(parsed)
            case "compile-partition":
                return try compilePartition(parsed)
            case "compile-population":
                return try compilePopulation(parsed)
            case "compile-surrogate":
                return try compileSurrogate(parsed)
            case "inspect-checkpoint":
                return try inspectCheckpoint(parsed)
            default:
                throw VivoCLIError.usage("unknown command '\(command)'")
            }
        } catch let error as VivoCLIError {
            writeStandardError("numivivo: \(error.localizedDescription)\n")
            if case .usage = error {
                writeStandardError("Run 'numivivo help' for usage.\n")
            }
            return error.exitCode
        } catch {
            writeStandardError("numivivo: \(error.localizedDescription)\n")
            return VivoCLIExit.software
        }
    }

    private func validate(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "validate <program.json>")
        try arguments.allowOptions([
            "fidelity", "relaxed-units", "relaxed-safety", "allow-no-termination", "diagnostics"
        ])
        let source = try read(arguments.positionals[0])
        let invocation = try VivoNativeCompilerBridge.validateProgram(
            source,
            fidelity: try fidelity(arguments.value("fidelity") ?? "F2"),
            strictUnits: !arguments.flag("relaxed-units"),
            strictSafety: !arguments.flag("relaxed-safety"),
            requireTermination: !arguments.flag("allow-no-termination")
        )
        try publishDiagnostics(invocation.diagnostics, arguments: arguments)
        return exitCode(invocation.status)
    }

    private func compile(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "compile <program.json> --output <program.vivopack>")
        try arguments.requireOption("output")
        try arguments.allowOptions([
            "output", "force", "fidelity", "relaxed-units", "relaxed-safety", "allow-no-termination", "diagnostics"
        ])
        let source = try read(arguments.positionals[0])
        let invocation = try VivoNativeCompilerBridge.compileProgram(
            source,
            fidelity: try fidelity(arguments.value("fidelity") ?? "F2"),
            strictUnits: !arguments.flag("relaxed-units"),
            strictSafety: !arguments.flag("relaxed-safety"),
            requireTermination: !arguments.flag("allow-no-termination")
        )
        try publishDiagnostics(invocation.diagnostics, arguments: arguments)
        guard invocation.succeeded else { return exitCode(invocation.status) }
        try write(
            invocation.primary,
            destination: arguments.requiredValue("output"),
            force: arguments.flag("force"),
            binary: true
        )
        return VivoCLIExit.success
    }

    private func inspect(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "inspect <program.vivopack>")
        try arguments.allowOptions(["skip-section-hashes", "output", "force"])
        let pack = try read(arguments.positionals[0])
        let inspection = VivoNativeCompilerBridge.inspectProgramPack(
            pack,
            verifySectionHashes: !arguments.flag("skip-section-hashes")
        )
        try publishPrimary(
            inspection.invocation.primary,
            arguments: arguments,
            force: arguments.flag("force")
        )
        return exitCode(inspection.invocation.status)
    }

    private func analyze(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "analyze <program.json>")
        try arguments.allowOptions([
            "fidelity", "maximum-dense-species", "maximum-dense-reactions",
            "maximum-dense-elements", "maximum-conservation-laws", "rank-tolerance",
            "conservation-tolerance", "output", "diagnostics", "force"
        ])
        let invocation = try VivoNativeCompilerBridge.analyzeProgram(
            try read(arguments.positionals[0]),
            fidelity: try fidelity(arguments.value("fidelity") ?? "F2"),
            maximumDenseSpecies: try arguments.uint32("maximum-dense-species", default: 2_048),
            maximumDenseReactions: try arguments.uint32("maximum-dense-reactions", default: 2_048),
            maximumDenseElements: try arguments.uint64("maximum-dense-elements", default: 16_777_216),
            maximumConservationLaws: try arguments.uint32("maximum-conservation-laws", default: 256),
            rankTolerance: try arguments.double("rank-tolerance", default: 1e-10),
            conservationTolerance: try arguments.double("conservation-tolerance", default: 1e-9)
        )
        try publishDiagnostics(invocation.diagnostics, arguments: arguments)
        if !invocation.primary.isEmpty {
            try publishPrimary(
                invocation.primary,
                arguments: arguments,
                force: arguments.flag("force")
            )
        }
        return exitCode(invocation.status)
    }

    private func synthesizeMechanisms(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "synthesize-mechanisms <problem.json> --library <library.json>")
        try arguments.requireOption("library")
        try arguments.allowOptions([
            "library", "output", "diagnostics", "force", "maximum-solutions",
            "maximum-visited-nodes", "allow-missing-shutdown", "allow-missing-monitor",
            "require-context-insulation", "require-resource-buffering",
            "require-distinct-orthogonality"
        ])
        let invocation = VivoNativeCompilerBridge.synthesizeMechanisms(
            problem: try read(arguments.positionals[0]),
            library: try read(arguments.requiredValue("library")),
            maximumSolutions: try arguments.optionalUInt32("maximum-solutions"),
            maximumVisitedNodes: try arguments.optionalUInt64("maximum-visited-nodes"),
            requireIndependentShutdown: arguments.flag("allow-missing-shutdown") ? false : nil,
            requireMonitor: arguments.flag("allow-missing-monitor") ? false : nil,
            requireContextInsulation: arguments.flag("require-context-insulation") ? true : nil,
            requireResourceBuffering: arguments.flag("require-resource-buffering") ? true : nil,
            requireDistinctOrthogonality: arguments.flag("require-distinct-orthogonality") ? true : nil
        )
        try publishDiagnostics(invocation.diagnostics, arguments: arguments)
        if !invocation.primary.isEmpty {
            try publishPrimary(
                invocation.primary,
                arguments: arguments,
                force: arguments.flag("force")
            )
        }
        return exitCode(invocation.status)
    }

    private func compileCampaign(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "compile-campaign <definition.json> --output <manifest.json>")
        try arguments.requireOption("output")
        try arguments.allowOptions(["output", "force"])
        let compiled = try VivoSourceArtifactCompiler.compileCampaign(
            from: try read(arguments.positionals[0])
        )
        try write(
            try VivoSourceArtifactCompiler.canonicalJSON(compiled),
            destination: arguments.requiredValue("output"),
            force: arguments.flag("force"),
            binary: false
        )
        writeStandardError("campaign \(compiled.artifactFingerprint)\n")
        return VivoCLIExit.success
    }

    private func verifyCampaign(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "verify-campaign <manifest.json>")
        try arguments.allowOptions(["output", "force"])
        let loaded = try VivoValidatedArtifactLoader.campaignManifest(
            from: try read(arguments.positionals[0])
        )
        let report: [String: Any] = [
            "valid": true,
            "sourceFingerprint": loaded.sourceFingerprint,
            "sourceBytes": loaded.sourceBytes,
            "manifestFingerprint": loaded.value.fingerprint,
            "ledgerHead": loaded.value.ledgerHead,
            "candidateCount": loaded.value.candidates.count,
            "jobCount": loaded.value.jobs.count
        ]
        try publishPrimary(
            try JSONSerialization.data(withJSONObject: report, options: [.prettyPrinted, .sortedKeys]),
            arguments: arguments,
            force: arguments.flag("force")
        )
        return VivoCLIExit.success
    }

    private func compilePartition(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "compile-partition <source.json> --output <model.json>")
        try arguments.requireOption("output")
        try arguments.allowOptions(["output", "force"])
        let compiled = try VivoSourceArtifactCompiler.physiologicalPartitionModel(
            from: try read(arguments.positionals[0])
        )
        try write(
            try VivoSourceArtifactCompiler.canonicalJSON(compiled),
            destination: arguments.requiredValue("output"),
            force: arguments.flag("force"),
            binary: false
        )
        writeStandardError("partition model \(compiled.artifactFingerprint)\n")
        return VivoCLIExit.success
    }

    private func compilePopulation(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "compile-population <source.json> --output <model.json>")
        try arguments.requireOption("output")
        try arguments.allowOptions(["output", "force"])
        let compiled = try VivoSourceArtifactCompiler.populationModel(
            from: try read(arguments.positionals[0])
        )
        try write(
            try VivoSourceArtifactCompiler.canonicalJSON(compiled),
            destination: arguments.requiredValue("output"),
            force: arguments.flag("force"),
            binary: false
        )
        writeStandardError("population model \(compiled.artifactFingerprint)\n")
        return VivoCLIExit.success
    }

    private func compileSurrogate(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "compile-surrogate <source.json> --output <contract.json>")
        try arguments.requireOption("output")
        try arguments.allowOptions(["output", "force"])
        let compiled = try VivoSourceArtifactCompiler.surrogateContract(
            from: try read(arguments.positionals[0])
        )
        try write(
            try VivoSourceArtifactCompiler.canonicalJSON(compiled),
            destination: arguments.requiredValue("output"),
            force: arguments.flag("force"),
            binary: false
        )
        writeStandardError("surrogate contract \(compiled.artifactFingerprint)\n")
        return VivoCLIExit.success
    }

    private func inspectCheckpoint(_ arguments: VivoCLIArguments) throws -> Int32 {
        try arguments.requirePositionals(1, usage: "inspect-checkpoint <checkpoint.vivocheckpoint>")
        try arguments.allowOptions(["require-signature", "output", "force"])
        let checkpoint = try VivoValidatedArtifactLoader.checkpoint(
            from: try read(arguments.positionals[0]),
            requireSignature: arguments.flag("require-signature")
        )
        let report = VivoCheckpointInspectionReport(
            valid: true,
            packageFingerprint: checkpoint.packageFingerprint,
            signed: checkpoint.signature != nil,
            manifest: checkpoint.manifest
        )
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        try publishPrimary(
            try encoder.encode(report),
            arguments: arguments,
            force: arguments.flag("force")
        )
        return VivoCLIExit.success
    }

    private func publishDiagnostics(
        _ data: Data,
        arguments: VivoCLIArguments
    ) throws {
        guard !data.isEmpty else { return }
        if let destination = arguments.value("diagnostics") {
            try write(
                data,
                destination: destination,
                force: arguments.flag("force"),
                binary: false
            )
        } else {
            FileHandle.standardError.write(data)
            if data.last != 0x0a { FileHandle.standardError.write(Data("\n".utf8)) }
        }
    }

    private func publishPrimary(
        _ data: Data,
        arguments: VivoCLIArguments,
        force: Bool
    ) throws {
        if let output = arguments.value("output") {
            try write(data, destination: output, force: force, binary: false)
        } else {
            FileHandle.standardOutput.write(data)
            if data.last != 0x0a { FileHandle.standardOutput.write(Data("\n".utf8)) }
        }
    }

    private func read(_ path: String) throws -> Data {
        if path == "-" {
            let data = FileHandle.standardInput.readDataToEndOfFile()
            guard data.count <= 64 * 1_024 * 1_024 else {
                throw VivoCLIError.data("standard input exceeds 64 MiB")
            }
            return data
        }
        let url = URL(fileURLWithPath: path)
        do {
            return try VivoValidatedArtifactLoader.data(at: url)
        } catch {
            throw VivoCLIError.input(path, error.localizedDescription)
        }
    }

    private func write(
        _ data: Data,
        destination: String,
        force: Bool,
        binary: Bool
    ) throws {
        if destination == "-" {
            FileHandle.standardOutput.write(data)
            if !binary, data.last != 0x0a {
                FileHandle.standardOutput.write(Data("\n".utf8))
            }
            return
        }
        let url = URL(fileURLWithPath: destination)
        let manager = FileManager.default
        if manager.fileExists(atPath: url.path), !force {
            throw VivoCLIError.output(destination, "file exists; pass --force to replace it")
        }
        do {
            try manager.createDirectory(
                at: url.deletingLastPathComponent(),
                withIntermediateDirectories: true
            )
            try data.write(to: url, options: [.atomic])
        } catch {
            throw VivoCLIError.output(destination, error.localizedDescription)
        }
    }

    private func fidelity(_ value: String) throws -> VivoFidelity {
        switch value.uppercased() {
        case "F0", "LOGIC": return .logic
        case "F1", "DETERMINISTIC": return .deterministic
        case "F2", "STOCHASTIC": return .stochastic
        case "F3", "SPATIAL": return .spatial
        case "F4", "TISSUE": return .tissue
        default: throw VivoCLIError.usage("invalid fidelity '\(value)'; use F0 through F4")
        }
    }

    private func exitCode(_ status: VivoNativeStatus) -> Int32 {
        switch status {
        case .ok: return VivoCLIExit.success
        case .invalidArgument: return VivoCLIExit.usage
        case .parseError, .validationError, .safetyRejected, .invalidPack: return VivoCLIExit.data
        case .resourceLimit, .outOfMemory: return VivoCLIExit.temporaryFailure
        case .compileError, .internalError, .unknown: return VivoCLIExit.software
        }
    }

    static let help = """
    NumiVivo \(NumiVivo.version)

    Usage:
      numivivo validate <program.json> [--fidelity F0|F1|F2|F3|F4]
      numivivo compile <program.json> --output <program.vivopack>
      numivivo inspect <program.vivopack>
      numivivo analyze <program.json>
      numivivo synthesize-mechanisms <problem.json> --library <library.json>
      numivivo compile-campaign <definition.json> --output <manifest.json>
      numivivo verify-campaign <manifest.json>
      numivivo compile-partition <source.json> --output <model.json>
      numivivo compile-population <source.json> --output <model.json>
      numivivo compile-surrogate <source.json> --output <contract.json>
      numivivo inspect-checkpoint <checkpoint.vivocheckpoint>
      numivivo version

    Common options:
      --output <path>       Write primary output; '-' writes to standard output.
      --diagnostics <path>  Write diagnostic JSON instead of standard error.
      --force               Replace an existing output file.

    This command line operates on computational model artifacts. It does not
    generate nucleotide sequences, laboratory procedures, or clinical advice.
    """ + "\n"
}

private struct VivoCheckpointInspectionReport: Encodable {
    let valid: Bool
    let packageFingerprint: String
    let signed: Bool
    let manifest: VivoCheckpointManifest
}

private struct VivoCLIArguments {
    let command: String?
    let positionals: [String]
    private let options: [String: String?]

    init(_ raw: [String]) throws {
        command = raw.first
        var positionals: [String] = []
        var options: [String: String?] = [:]
        var index = raw.isEmpty ? 0 : 1
        var optionParsing = true

        while index < raw.count {
            let token = raw[index]
            if optionParsing, token == "--" {
                optionParsing = false
                index += 1
                continue
            }
            if optionParsing, token.hasPrefix("--") {
                let body = String(token.dropFirst(2))
                guard !body.isEmpty else {
                    throw VivoCLIError.usage("empty option")
                }
                let pieces = body.split(separator: "=", maxSplits: 1, omittingEmptySubsequences: false)
                let name = String(pieces[0])
                guard options[name] == nil else {
                    throw VivoCLIError.usage("option --\(name) was provided more than once")
                }
                if pieces.count == 2 {
                    guard !pieces[1].isEmpty else {
                        throw VivoCLIError.usage("option --\(name) has an empty value")
                    }
                    options[name] = String(pieces[1])
                } else if index + 1 < raw.count, !raw[index + 1].hasPrefix("--") {
                    options[name] = raw[index + 1]
                    index += 1
                } else {
                    options[name] = .some(nil)
                }
            } else if optionParsing, token.hasPrefix("-") {
                throw VivoCLIError.usage("short options other than -h are not supported: \(token)")
            } else {
                positionals.append(token)
            }
            index += 1
        }
        self.positionals = positionals
        self.options = options
    }

    func requirePositionals(_ count: Int, usage: String) throws {
        guard positionals.count == count else {
            throw VivoCLIError.usage("expected: numivivo \(usage)")
        }
    }

    func allowOptions(_ allowed: Set<String>) throws {
        if let unknown = options.keys.sorted().first(where: { !allowed.contains($0) }) {
            throw VivoCLIError.usage("unknown option --\(unknown)")
        }
    }

    func requireOption(_ name: String) throws {
        guard options.keys.contains(name), value(name) != nil else {
            throw VivoCLIError.usage("missing required option --\(name)")
        }
    }

    func value(_ name: String) -> String? {
        guard let entry = options[name] else { return nil }
        return entry
    }

    func requiredValue(_ name: String) throws -> String {
        guard let value = value(name) else {
            throw VivoCLIError.usage("option --\(name) requires a value")
        }
        return value
    }

    func flag(_ name: String) -> Bool {
        guard options.keys.contains(name) else { return false }
        return options[name] == nil
    }

    func optionalUInt32(_ name: String) throws -> UInt32? {
        guard let value = value(name) else { return nil }
        guard let parsed = UInt32(value), parsed > 0 else {
            throw VivoCLIError.usage("--\(name) requires a positive UInt32")
        }
        return parsed
    }

    func optionalUInt64(_ name: String) throws -> UInt64? {
        guard let value = value(name) else { return nil }
        guard let parsed = UInt64(value), parsed > 0 else {
            throw VivoCLIError.usage("--\(name) requires a positive UInt64")
        }
        return parsed
    }

    func uint32(_ name: String, default defaultValue: UInt32) throws -> UInt32 {
        try optionalUInt32(name) ?? defaultValue
    }

    func uint64(_ name: String, default defaultValue: UInt64) throws -> UInt64 {
        try optionalUInt64(name) ?? defaultValue
    }

    func double(_ name: String, default defaultValue: Double) throws -> Double {
        guard let value = value(name) else { return defaultValue }
        guard let parsed = Double(value), parsed.isFinite, parsed > 0 else {
            throw VivoCLIError.usage("--\(name) requires a positive finite number")
        }
        return parsed
    }
}

private enum VivoCLIError: Error, LocalizedError {
    case usage(String)
    case data(String)
    case input(String, String)
    case output(String, String)

    var errorDescription: String? {
        switch self {
        case .usage(let message): return message
        case .data(let message): return message
        case .input(let path, let message): return "cannot read \(path): \(message)"
        case .output(let path, let message): return "cannot write \(path): \(message)"
        }
    }

    var exitCode: Int32 {
        switch self {
        case .usage: return VivoCLIExit.usage
        case .data: return VivoCLIExit.data
        case .input: return VivoCLIExit.noInput
        case .output: return VivoCLIExit.cannotCreate
        }
    }
}

private enum VivoCLIExit {
    static let success: Int32 = 0
    static let usage: Int32 = 64
    static let data: Int32 = 65
    static let noInput: Int32 = 66
    static let software: Int32 = 70
    static let cannotCreate: Int32 = 73
    static let temporaryFailure: Int32 = 75
}

private func writeStandardOutput(_ text: String) {
    FileHandle.standardOutput.write(Data(text.utf8))
}

private func writeStandardError(_ text: String) {
    FileHandle.standardError.write(Data(text.utf8))
}
