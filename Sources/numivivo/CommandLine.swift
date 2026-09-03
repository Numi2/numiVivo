import Foundation
import NumiVivoKit
#if canImport(Darwin)
import Darwin
#endif

struct NumiVivoCommandLine {
    static func main() async {
        do {
            let invocation = try Invocation(arguments: Array(CommandLine.arguments.dropFirst()))
            switch invocation.command {
            case "help", "--help", "-h", nil:
                writeStandardOutput(helpText)
            case "version", "--version":
                writeStandardOutput("numivivo 0.1.0\n")
            case "validate":
                try validate(invocation)
            case "compile":
                try compile(invocation)
            case "inspect":
                try inspect(invocation)
            case "run":
                try await run(invocation)
            default:
                throw CLIError.usage("Unknown command '\(invocation.command!)'.\n\n\(helpText)")
            }
        } catch let error as CLIError {
            writeStandardError("numivivo: \(error.description)\n")
            exit(error.exitCode)
        } catch {
            writeStandardError("numivivo: \(error)\n")
            exit(1)
        }
    }

    private static func validate(_ invocation: Invocation) throws {
        let sourceURL = try invocation.requiredInputURL()
        let configuration = try compilerConfiguration(invocation)
        let report = try VivoCompiler.validate(
            file: sourceURL,
            configuration: configuration
        )
        writeStandardOutput(report)
        writeStandardOutput(Data("\n".utf8))
    }

    private static func compile(_ invocation: Invocation) throws {
        let sourceURL = try invocation.requiredInputURL()
        let outputURL = try invocation.outputURL(
            defaultExtension: "nvivo",
            replacingExtensionOf: sourceURL
        )
        let configuration = try compilerConfiguration(invocation)
        let output = try VivoCompiler.compile(
            file: sourceURL,
            configuration: configuration
        )
        try ensureParentDirectory(of: outputURL)
        try output.programPack.data.write(to: outputURL, options: [.atomic])

        if let reportPath = invocation.value("report") {
            let reportURL = URL(fileURLWithPath: reportPath)
            try ensureParentDirectory(of: reportURL)
            try output.reportJSON.write(to: reportURL, options: [.atomic])
        }

        let response: [String: Any] = [
            "status": "compiled",
            "source": sourceURL.path,
            "output": outputURL.path,
            "bytes": output.programPack.data.count,
            "sourceFingerprint": output.programPack.header.sourceFingerprint,
            "contentFingerprint": output.programPack.header.contentFingerprint,
            "fidelity": output.programPack.header.fidelity.rawValue,
            "species": output.programPack.runtimeContract.speciesCount,
            "reactions": output.programPack.runtimeContract.reactionCount,
            "rules": output.programPack.runtimeContract.ruleCount,
            "monitors": output.programPack.runtimeContract.monitorCount
        ]
        try writeJSON(response)
    }

    private static func inspect(_ invocation: Invocation) throws {
        let inputURL = try invocation.requiredInputURL()
        let pack = try VivoProgramPack(contentsOf: inputURL)
        let sectionValues: [[String: Any]] = pack.sections.values
            .sorted { $0.type.rawValue < $1.type.rawValue }
            .map { section in
                [
                    "type": section.type.rawValue,
                    "name": String(describing: section.type),
                    "offset": section.offset,
                    "size": section.size,
                    "stride": section.stride,
                    "count": section.count,
                    "alignment": section.alignment,
                    "fingerprint": section.fingerprint
                ]
            }
        let speciesValues: [[String: Any]] = pack.species.map { species in
            [
                "index": species.index,
                "id": species.id,
                "compartment": species.compartment,
                "unit": species.unit,
                "flags": species.flags,
                "initialValue": species.initialValue,
                "minimum": species.minimum,
                "maximum": species.maximum
            ]
        }
        let cohortValues: [[String: Any]] = pack.cohorts.map { cohort in
            [
                "index": cohort.index,
                "reactionOffset": cohort.reactionOffset,
                "reactionCount": cohort.reactionCount,
                "rateLaw": cohort.rateLaw,
                "flags": cohort.flags,
                "maximumStableStep": cohort.maximumStableStep,
                "stiffnessEstimate": cohort.stiffnessEstimate,
                "preferredThreads": cohort.preferredThreads
            ]
        }
        var response: [String: Any] = [
            "header": [
                "major": pack.header.major,
                "minor": pack.header.minor,
                "compilerABI": pack.header.compilerABI,
                "flags": pack.header.flags,
                "fidelity": pack.header.fidelity.rawValue,
                "sectionCount": pack.header.sectionCount,
                "totalBytes": pack.header.totalBytes,
                "sourceFingerprint": pack.header.sourceFingerprint,
                "contentFingerprint": pack.header.contentFingerprint
            ],
            "contract": [
                "speciesCount": pack.runtimeContract.speciesCount,
                "parameterCount": pack.runtimeContract.parameterCount,
                "reactionCount": pack.runtimeContract.reactionCount,
                "ruleCount": pack.runtimeContract.ruleCount,
                "monitorCount": pack.runtimeContract.monitorCount,
                "cohortCount": pack.runtimeContract.cohortCount,
                "temporalStateCount": pack.runtimeContract.temporalStateCount,
                "maximumExpressionStack": pack.runtimeContract.maximumExpressionStack,
                "featureFlags": pack.runtimeContract.featureFlags,
                "authoritativeScalarBytes": pack.runtimeContract.authoritativeScalarBytes,
                "randomStreamVersion": pack.runtimeContract.randomStreamVersion
            ],
            "sections": sectionValues,
            "species": speciesValues,
            "cohorts": cohortValues
        ]
        if let manifest = pack.manifestJSON,
           let object = try? JSONSerialization.jsonObject(with: manifest) {
            response["manifest"] = object
        }
        try writeJSON(response)
    }

    private static func run(_ invocation: Invocation) async throws {
        let inputURL = try invocation.requiredInputURL()
        let mode = try stepMode(invocation.value("mode") ?? "stochastic")
        let fidelity: VivoFidelity = mode == .stochastic ? .f2Stochastic : .f1Deterministic
        let pack: VivoProgramPack
        if inputURL.pathExtension.lowercased() == "nvivo" {
            pack = try VivoProgramPack(contentsOf: inputURL)
        } else {
            pack = try VivoCompiler.compile(
                file: inputURL,
                configuration: VivoCompilerConfiguration(
                    fidelity: fidelity,
                    strictUnits: !invocation.flag("relaxed-units"),
                    strictSafety: !invocation.flag("relaxed-safety"),
                    deterministicPack: true,
                    permitHypotheticalParameters: invocation.flag("permit-hypothetical"),
                    requireTermination: !invocation.flag("allow-unterminated")
                )
            ).programPack
        }

        let cells = try invocation.integer("cells", default: 1, minimum: 1)
        let duration = try invocation.double("duration", default: 1.0, minimumExclusive: 0)
        let deltaTime = try invocation.double("dt", default: 0.01, minimumExclusive: 0)
        let minimumSubstep = try invocation.double(
            "minimum-substep",
            default: min(deltaTime, 0.001),
            minimumExclusive: 0
        )
        let maximumSubstep = try invocation.double(
            "maximum-substep",
            default: deltaTime,
            minimumExclusive: 0
        )
        let seed = try invocation.unsignedInteger("seed", default: 0)
        let sampleEvery = try invocation.integer("sample-every", default: 0, minimum: 0)
        let eventCapacity = try invocation.integer("event-capacity", default: 65_536, minimum: 1)
        let runtime = try NumiVivoRuntime(
            programPack: pack,
            configuration: VivoRuntimeConfiguration(
                cellCapacity: cells,
                activeCellCount: cells,
                mode: mode,
                seed: seed,
                maximumSubstep: Float(maximumSubstep),
                minimumSubstep: Float(minimumSubstep),
                maximumInternalSubsteps: try invocation.integer(
                    "maximum-internal-substeps",
                    default: 4_096,
                    minimum: 1
                ),
                eventCapacity: eventCapacity,
                maximumPrivateBytes: try invocation.integer(
                    "maximum-private-bytes",
                    default: 8 * 1_024 * 1_024 * 1_024,
                    minimum: 1
                ),
                maximumDelayBytes: try invocation.integer(
                    "maximum-delay-bytes",
                    default: 2 * 1_024 * 1_024 * 1_024,
                    minimum: 0
                )
            )
        )
        let observations = observationIDs(invocation, pack: pack)
        let staticInputs = try parseInputs(invocation.values("input"), pack: pack, cellCount: cells)
        let writer = try OutputWriter(path: invocation.value("output"))
        defer { writer.close() }

        let memory: [String: Any] = [
            "type": "runtime",
            "device": runtime.deviceInfo.name,
            "registryID": runtime.deviceInfo.registryID,
            "programFingerprint": pack.header.contentFingerprint,
            "mode": mode.rawValue,
            "cellCount": cells,
            "privateBytes": runtime.memoryReport.privateBytes,
            "sharedBytes": runtime.memoryReport.sharedBytes,
            "delaySlotCount": runtime.memoryReport.delaySlotCount
        ]
        try writer.writeJSONLine(memory)

        let totalSteps = Int(ceil(duration / deltaTime))
        var elapsed = 0.0
        for stepIndex in 0..<totalSteps {
            let remaining = duration - elapsed
            let stepDelta = min(deltaTime, remaining)
            if stepDelta <= 0 { break }
            let result = try await runtime.step(
                deltaTime: Float(stepDelta),
                inputs: staticInputs
            )
            elapsed += stepDelta

            let shouldSample = sampleEvery > 0
                ? ((stepIndex + 1) % sampleEvery == 0 || stepIndex + 1 == totalSteps)
                : stepIndex + 1 == totalSteps
            var stateObject: [String: Any] = [:]
            if shouldSample, !observations.isEmpty {
                let slices = try await runtime.readState(species: observations)
                for slice in slices { stateObject[slice.species] = slice.values }
            }

            let stepObject: [String: Any] = [
                "type": "step",
                "status": result.status.rawValue,
                "logicalStep": result.committedLogicalStep,
                "time": result.committedTime,
                "requestedDeltaTime": result.requestedDeltaTime,
                "executedSubsteps": result.executedSubsteps,
                "stateVersion": result.stateVersion,
                "diagnosticFlags": result.diagnostics.flags,
                "nonFiniteCount": result.diagnostics.nonFiniteCount,
                "boundViolationCount": result.diagnostics.boundViolationCount,
                "monitorViolationCount": result.diagnostics.monitorViolationCount,
                "stochasticFallbackCount": result.diagnostics.stochasticFallbackCount,
                "fluxTruncationCount": result.diagnostics.fluxTruncationCount,
                "eventOverflowCount": result.diagnostics.eventOverflowCount,
                "events": result.events.map { event in
                    [
                        "cell": event.cellIndex,
                        "kind": event.kind,
                        "subject": event.subject,
                        "logicalStep": event.logicalStep,
                        "value0": event.value0,
                        "value1": event.value1,
                        "flags": event.flags
                    ]
                },
                "state": stateObject
            ]
            if shouldSample || result.status != .committed || !result.events.isEmpty {
                try writer.writeJSONLine(stepObject)
            }
            if result.status != .committed { break }
        }
    }

    private static func compilerConfiguration(_ invocation: Invocation) throws -> VivoCompilerConfiguration {
        let fidelity = try parseFidelity(invocation.value("fidelity") ?? "F2")
        return VivoCompilerConfiguration(
            fidelity: fidelity,
            strictUnits: !invocation.flag("relaxed-units"),
            strictSafety: !invocation.flag("relaxed-safety"),
            deterministicPack: true,
            permitHypotheticalParameters: invocation.flag("permit-hypothetical"),
            requireTermination: !invocation.flag("allow-unterminated")
        )
    }

    private static func parseFidelity(_ value: String) throws -> VivoFidelity {
        switch value.uppercased() {
        case "F0", "LOGIC": return .f0Logic
        case "F1", "DETERMINISTIC": return .f1Deterministic
        case "F2", "STOCHASTIC": return .f2Stochastic
        case "F3", "SPATIAL": return .f3Spatial
        case "F4", "TISSUE": return .f4Tissue
        default: throw CLIError.usage("Unknown fidelity '\(value)'. Expected F0, F1, F2, F3, or F4.")
        }
    }

    private static func stepMode(_ value: String) throws -> VivoStepMode {
        switch value.lowercased() {
        case "deterministic", "f1": return .deterministic
        case "stochastic", "f2": return .stochastic
        default: throw CLIError.usage("Unknown mode '\(value)'. Expected deterministic or stochastic.")
        }
    }

    private static func observationIDs(
        _ invocation: Invocation,
        pack: VivoProgramPack
    ) -> [String] {
        if let value = invocation.value("observe") {
            return value.split(separator: ",").map(String.init).filter { !$0.isEmpty }
        }
        let outputs = pack.species.filter(\.isOutput).map(\.id)
        return outputs.isEmpty ? Array(pack.species.prefix(8).map(\.id)) : outputs
    }

    private static func parseInputs(
        _ values: [String],
        pack: VivoProgramPack,
        cellCount: Int
    ) throws -> [VivoInputUpdate] {
        let speciesByID = pack.speciesByID
        return try values.map { value in
            guard let equals = value.lastIndex(of: "=") else {
                throw CLIError.usage("Input '\(value)' must use SIGNAL[:CELL]=VALUE.")
            }
            let target = String(value[..<equals])
            let numeric = String(value[value.index(after: equals)...])
            guard let parsedValue = Float(numeric), parsedValue.isFinite else {
                throw CLIError.usage("Input '\(value)' has an invalid numeric value.")
            }
            let components = target.split(separator: ":", omittingEmptySubsequences: false)
            let signal = String(components[0])
            guard let species = speciesByID[signal], species.isExternallyOwned else {
                throw CLIError.usage("Input signal '\(signal)' is not externally owned by the program.")
            }
            let cell: UInt32
            if components.count == 1 {
                cell = 0
            } else if components.count == 2,
                      let parsedCell = UInt32(components[1]),
                      parsedCell < UInt32(cellCount) {
                cell = parsedCell
            } else {
                throw CLIError.usage("Input '\(value)' has an invalid cell index.")
            }
            return VivoInputUpdate(signal: signal, cellIndex: cell, value: parsedValue)
        }
    }

    private static func ensureParentDirectory(of url: URL) throws {
        let parent = url.deletingLastPathComponent()
        try FileManager.default.createDirectory(
            at: parent,
            withIntermediateDirectories: true
        )
    }

    private static func writeJSON(_ object: Any) throws {
        let data = try JSONSerialization.data(
            withJSONObject: object,
            options: [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        )
        writeStandardOutput(data)
        writeStandardOutput(Data("\n".utf8))
    }

    private static func writeStandardOutput(_ string: String) {
        writeStandardOutput(Data(string.utf8))
    }

    private static func writeStandardOutput(_ data: Data) {
        FileHandle.standardOutput.write(data)
    }

    private static func writeStandardError(_ string: String) {
        FileHandle.standardError.write(Data(string.utf8))
    }

    private static let helpText = """
    NumiVivo molecular-program compiler and Apple-silicon runtime

    USAGE
      numivivo validate PROGRAM.json [--fidelity F0|F1|F2|F3|F4]
      numivivo compile PROGRAM.json [-o PROGRAM.nvivo] [--report REPORT.json]
      numivivo inspect PROGRAM.nvivo
      numivivo run PROGRAM.json|PROGRAM.nvivo --cells N --duration SECONDS --dt SECONDS

    COMMON COMPILER OPTIONS
      --fidelity VALUE             Minimum execution fidelity; default F2
      --relaxed-units              Permit implicit unitless literal interpretation
      --relaxed-safety             Downgrade selected safety findings
      --permit-hypothetical        Permit explicitly hypothetical parameters
      --allow-unterminated         Do not require a declared termination path

    RUN OPTIONS
      --mode deterministic|stochastic
      --seed UINT64
      --observe ID,ID,...
      --input SIGNAL[:CELL]=VALUE  Repeatable static input staged every logical step
      --sample-every N             Emit sampled state every N steps; default final only
      --output FILE.ndjson         Write newline-delimited JSON instead of stdout
      --minimum-substep SECONDS
      --maximum-substep SECONDS
      --maximum-internal-substeps N
      --event-capacity N
      --maximum-private-bytes N
      --maximum-delay-bytes N
    """
}

private struct Invocation {
    let command: String?
    let positionals: [String]
    private let options: [String: [String]]
    private let flags: Set<String>

    init(arguments: [String]) throws {
        self.command = arguments.first
        var positionals: [String] = []
        var options: [String: [String]] = [:]
        var flags: Set<String> = []
        var index = command == nil ? 0 : 1

        while index < arguments.count {
            let argument = arguments[index]
            if argument == "--" {
                positionals.append(contentsOf: arguments[(index + 1)...])
                break
            }
            if argument == "-o" {
                guard index + 1 < arguments.count else {
                    throw CLIError.usage("-o requires a path.")
                }
                options["output", default: []].append(arguments[index + 1])
                index += 2
                continue
            }
            if argument.hasPrefix("--") {
                let body = String(argument.dropFirst(2))
                if let equals = body.firstIndex(of: "=") {
                    let name = String(body[..<equals])
                    let value = String(body[body.index(after: equals)...])
                    options[name, default: []].append(value)
                    index += 1
                    continue
                }
                if index + 1 < arguments.count, !arguments[index + 1].hasPrefix("-") {
                    options[body, default: []].append(arguments[index + 1])
                    index += 2
                } else {
                    flags.insert(body)
                    index += 1
                }
                continue
            }
            positionals.append(argument)
            index += 1
        }
        self.positionals = positionals
        self.options = options
        self.flags = flags
    }

    func value(_ name: String) -> String? { options[name]?.last }
    func values(_ name: String) -> [String] { options[name] ?? [] }
    func flag(_ name: String) -> Bool { flags.contains(name) }

    func requiredInputURL() throws -> URL {
        guard let path = positionals.first else {
            throw CLIError.usage("A program path is required.")
        }
        return URL(fileURLWithPath: path)
    }

    func outputURL(defaultExtension: String, replacingExtensionOf input: URL) throws -> URL {
        if let path = value("output") { return URL(fileURLWithPath: path) }
        return input.deletingPathExtension().appendingPathExtension(defaultExtension)
    }

    func integer(_ name: String, default fallback: Int, minimum: Int) throws -> Int {
        guard let value = value(name) else { return fallback }
        guard let parsed = Int(value), parsed >= minimum else {
            throw CLIError.usage("--\(name) must be an integer greater than or equal to \(minimum).")
        }
        return parsed
    }

    func unsignedInteger(_ name: String, default fallback: UInt64) throws -> UInt64 {
        guard let value = value(name) else { return fallback }
        guard let parsed = UInt64(value) else {
            throw CLIError.usage("--\(name) must be an unsigned integer.")
        }
        return parsed
    }

    func double(_ name: String, default fallback: Double, minimumExclusive minimum: Double) throws -> Double {
        guard let value = value(name) else { return fallback }
        guard let parsed = Double(value), parsed.isFinite, parsed > minimum else {
            throw CLIError.usage("--\(name) must be a finite number greater than \(minimum).")
        }
        return parsed
    }
}

private enum CLIError: Error, CustomStringConvertible {
    case usage(String)
    case data(String)

    var exitCode: Int32 {
        switch self {
        case .usage: return 64
        case .data: return 65
        }
    }

    var description: String {
        switch self {
        case .usage(let message), .data(let message): return message
        }
    }
}

private final class OutputWriter {
    private let handle: FileHandle
    private let shouldClose: Bool

    init(path: String?) throws {
        if let path {
            let url = URL(fileURLWithPath: path)
            try FileManager.default.createDirectory(
                at: url.deletingLastPathComponent(),
                withIntermediateDirectories: true
            )
            FileManager.default.createFile(atPath: url.path, contents: nil)
            guard let handle = FileHandle(forWritingAtPath: url.path) else {
                throw CLIError.data("Unable to open output file '\(path)'.")
            }
            self.handle = handle
            self.shouldClose = true
        } else {
            self.handle = .standardOutput
            self.shouldClose = false
        }
    }

    func writeJSONLine(_ object: Any) throws {
        var data = try JSONSerialization.data(
            withJSONObject: object,
            options: [.sortedKeys, .withoutEscapingSlashes]
        )
        data.append(0x0a)
        handle.write(data)
    }

    func close() {
        if shouldClose { try? handle.close() }
    }
}
