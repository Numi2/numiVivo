import Foundation
import NumiVivoKit

struct VivoForceFieldCLICommands {
    static func handles(_ name: String?) -> Bool {
        guard let name else { return false }
        return ["forcefield-validate", "forcefield-assign", "forcefield-compile", "forcefield-help"].contains(name)
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw CLIError.usage("missing force-field command") }
            if command == "forcefield-help" { FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0 }
            let args = try Arguments(Array(arguments.dropFirst()))
            switch command {
            case "forcefield-validate": try validate(args)
            case "forcefield-assign": try assign(args)
            case "forcefield-compile": try compile(args)
            default: throw CLIError.usage("unknown force-field command")
            }
            return 0
        } catch let error as CLIError {
            diagnostic(error.description); return 64
        } catch let error as VivoArtifactValidationError {
            diagnostic(error.description); return 65
        } catch {
            diagnostic(error.localizedDescription); return 70
        }
    }

    private func validate(_ args: Arguments) throws {
        try args.requirePositionals(1); try args.allow(["output", "force"])
        let library = try loadLibrary(args.positionals[0])
        try library.validate()
        let report = ValidationReport(schemaVersion: 1, fingerprint: try library.fingerprint(),
                                      identifier: library.identifier, version: library.version,
                                      atomTypeCount: UInt32(library.atomTypes.count),
                                      residueTemplateCount: UInt32(library.residueTemplates.count),
                                      bondParameterCount: UInt32(library.bondParameters.count),
                                      angleParameterCount: UInt32(library.angleParameters.count),
                                      torsionParameterCount: UInt32(library.torsionParameters.count))
        try write(try VivoCanonicalJSON.encode(report), to: args.options["output"] ?? "-", force: args.force)
    }

    private func assign(_ args: Arguments) throws {
        try args.requirePositionals(1); try args.allow(["library", "output", "force"])
        let structure = try loadStructure(args.positionals[0]).structure
        let library = try loadLibrary(try args.required("library"))
        let assignment = try VivoResidueTemplateAssigner.assign(structure: structure, library: library)
        try write(try VivoCanonicalJSON.encode(assignment), to: try args.required("output"), force: args.force)
        diagnostic("assigned=\(structure.atoms.count - assignment.unresolvedAtoms.count) unresolved=\(assignment.unresolvedAtoms.count)")
    }

    private func compile(_ args: Arguments) throws {
        try args.requirePositionals(1)
        try args.allow(["library", "assignment", "allow-missing-bonds", "allow-missing-angles", "allow-missing-torsions", "output", "report", "force"])
        let document = try loadStructure(args.positionals[0])
        let library = try loadLibrary(try args.required("library"))
        let assignment: VivoForceFieldAssignment?
        if let path = args.options["assignment"] {
            assignment = try VivoCanonicalJSON.decode(VivoForceFieldAssignment.self,
                                                       from: Data(contentsOf: URL(fileURLWithPath: path), options: [.mappedIfSafe]))
        } else { assignment = nil }
        let compiled = try VivoForceFieldCompiler.compile(
            structure: document.structure, library: library, assignment: assignment,
            options: .init(requireAllBondParameters: !args.flags.contains("allow-missing-bonds"),
                           requireAllAngleParameters: !args.flags.contains("allow-missing-angles"),
                           requireAllProperTorsions: !args.flags.contains("allow-missing-torsions"))
        )
        try write(try VivoCanonicalJSON.encode(compiled.system), to: try args.required("output"), force: args.force)
        if let report = args.options["report"] {
            try write(try VivoCanonicalJSON.encode(compiled.report), to: report, force: args.force)
        }
        diagnostic("classical-system \(try compiled.system.fingerprint()) particles=\(compiled.report.particleCount) bonds=\(compiled.report.bondCount) angles=\(compiled.report.angleCount) torsions=\(compiled.report.torsionTermCount)")
    }

    private func loadLibrary(_ path: String) throws -> VivoForceFieldLibrary {
        let data = try Data(contentsOf: URL(fileURLWithPath: path), options: [.mappedIfSafe])
        return try VivoCanonicalJSON.decode(VivoForceFieldLibrary.self, from: data)
    }
    private func loadStructure(_ path: String) throws -> VivoMolecularStructureDocument {
        try VivoMolecularStructureDocument.decode(Data(contentsOf: URL(fileURLWithPath: path), options: [.mappedIfSafe]))
    }
    private func write(_ data: Data, to path: String, force: Bool) throws {
        if path == "-" {
            FileHandle.standardOutput.write(data); if data.last != 0x0a { FileHandle.standardOutput.write(Data([0x0a])) }; return
        }
        if FileManager.default.fileExists(atPath: path), !force { throw CLIError.usage("output exists: \(path); use --force") }
        let url = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        try data.write(to: url, options: .atomic)
    }
    private func diagnostic(_ message: String) {
        FileHandle.standardError.write(Data(("numivivo forcefield: " + message + "\n").utf8))
    }

    private struct ValidationReport: Encodable {
        let schemaVersion: UInt32
        let fingerprint: VivoFingerprint
        let identifier: String
        let version: String
        let atomTypeCount: UInt32
        let residueTemplateCount: UInt32
        let bondParameterCount: UInt32
        let angleParameterCount: UInt32
        let torsionParameterCount: UInt32
    }
    private enum CLIError: Error, CustomStringConvertible {
        case usage(String)
        var description: String { switch self { case .usage(let value): value } }
    }
    private struct Arguments {
        var positionals: [String] = []
        var options: [String: String] = [:]
        var flags = Set<String>()
        var force: Bool { flags.contains("force") }
        init(_ raw: [String]) throws {
            let flagNames: Set<String> = ["force", "allow-missing-bonds", "allow-missing-angles", "allow-missing-torsions"]
            var i = 0
            while i < raw.count {
                let token = raw[i]
                if token.hasPrefix("--") {
                    let key = String(token.dropFirst(2))
                    if flagNames.contains(key) {
                        guard flags.insert(key).inserted else { throw CLIError.usage("duplicate --\(key)") }
                    } else {
                        guard options[key] == nil, i + 1 < raw.count, !raw[i + 1].hasPrefix("--") else {
                            throw CLIError.usage("duplicate option or missing value: --\(key)")
                        }
                        options[key] = raw[i + 1]; i += 1
                    }
                } else { positionals.append(token) }
                i += 1
            }
        }
        func requirePositionals(_ count: Int) throws {
            guard positionals.count == count else { throw CLIError.usage("expected exactly \(count) positional argument(s)") }
        }
        func allow(_ allowed: Set<String>) throws {
            let supplied = Set(options.keys).union(flags)
            guard supplied.isSubset(of: allowed) else { throw CLIError.usage("unknown option(s): \(supplied.subtracting(allowed).sorted().joined(separator: ", "))") }
        }
        func required(_ key: String) throws -> String {
            guard let value = options[key], !value.isEmpty else { throw CLIError.usage("missing --\(key)") }
            return value
        }
    }

    static let help = """
    Force-field foundation

      numivivo forcefield-validate <library.json> [--output <report.json>]
      numivivo forcefield-assign <structure.json> --library <library.json> --output <assignment.json>
      numivivo forcefield-compile <structure.json> --library <library.json> --output <system.json>
          [--assignment <assignment.json>] [--report <report.json>]
          [--allow-missing-bonds] [--allow-missing-angles] [--allow-missing-torsions] [--force]

    All parameters use NumiVivo canonical MD units: nm, ps, dalton, elementary
    charge and kJ/mol. Strict compilation rejects missing interaction parameters by default.
    """ + "\n"
}
