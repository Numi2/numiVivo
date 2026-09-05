import Foundation
import NumiVivoKit

struct VivoStructureCLICommands {
    static func handles(_ name: String?) -> Bool {
        guard let name else { return false }
        return ["structure-import", "structure-export", "structure-inspect", "structure-select", "structure-help"].contains(name)
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw CLIError.usage("missing structure command") }
            if command == "structure-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let parsed = try Arguments(Array(arguments.dropFirst()))
            switch command {
            case "structure-import": return try importStructure(parsed)
            case "structure-export": return try exportStructure(parsed)
            case "structure-inspect": return try inspectStructure(parsed)
            case "structure-select": return try selectAtoms(parsed)
            default: throw CLIError.usage("unknown structure command")
            }
        } catch let error as CLIError {
            diagnostic(error.description); return 64
        } catch let error as VivoArtifactValidationError {
            diagnostic(error.description); return 65
        } catch {
            diagnostic(error.localizedDescription); return 70
        }
    }

    private func importStructure(_ args: Arguments) throws -> Int32 {
        try args.requirePositionals(1)
        try args.allow(["format", "identifier", "output", "force"])
        let path = args.positionals[0]
        let format = try resolvedFormat(args.options["format"], path: path)
        let data = try Data(contentsOf: URL(fileURLWithPath: path), options: [.mappedIfSafe])
        let identifier = args.options["identifier"] ?? URL(fileURLWithPath: path).deletingPathExtension().lastPathComponent
        let imported = try VivoStructureCodec.decode(data, format: format, identifier: identifier)
        let document = try VivoMolecularStructureDocument(structure: imported.structure,
                                                          sourceFingerprint: imported.sourceFingerprint,
                                                          sourceFormat: imported.format)
        try write(try document.canonicalData(), to: args.options["output"] ?? "-", force: args.force)
        diagnostic("structure \(document.structureFingerprint) atoms=\(document.structure.atoms.count) bonds=\(document.structure.bonds.count)")
        return 0
    }

    private func exportStructure(_ args: Arguments) throws -> Int32 {
        try args.requirePositionals(1)
        try args.allow(["format", "output", "conformer", "force"])
        let output = try args.required("output")
        let format = try resolvedFormat(args.options["format"], path: output)
        let document = try VivoMolecularStructureDocument.decode(Data(contentsOf: URL(fileURLWithPath: args.positionals[0]), options: [.mappedIfSafe]))
        let conformer = try args.nonnegativeInt("conformer", default: 0)
        let bytes = try VivoStructureCodec.encode(document.structure, format: format, conformerIndex: conformer)
        try write(bytes, to: output, force: args.force)
        return 0
    }

    private func inspectStructure(_ args: Arguments) throws -> Int32 {
        try args.requirePositionals(1)
        try args.allow(["output", "force"])
        let document = try VivoMolecularStructureDocument.decode(Data(contentsOf: URL(fileURLWithPath: args.positionals[0]), options: [.mappedIfSafe]))
        let report = try VivoStructureValidator.validate(document.structure)
        let payload = InspectReport(schemaVersion: 1,
                                    structureFingerprint: document.structureFingerprint,
                                    documentFingerprint: try document.fingerprint(),
                                    sourceFingerprint: document.sourceFingerprint,
                                    sourceFormat: document.sourceFormat,
                                    identifier: document.structure.identifier,
                                    atomCount: report.atomCount, bondCount: report.bondCount,
                                    residueCount: report.residueCount, chainCount: report.chainCount,
                                    conformerCount: report.conformerCount,
                                    connectedComponents: report.connectedComponents,
                                    hasPeriodicCell: report.hasPeriodicCell,
                                    elements: Dictionary(grouping: document.structure.atoms, by: { $0.element.symbol })
                                        .mapValues { UInt32($0.count) })
        try write(try VivoCanonicalJSON.encode(payload), to: args.options["output"] ?? "-", force: args.force)
        return 0
    }

    private func selectAtoms(_ args: Arguments) throws -> Int32 {
        try args.requirePositionals(1)
        try args.allow(["query", "conformer", "no-periodic", "output", "force"])
        let query = try args.required("query")
        let document = try VivoMolecularStructureDocument.decode(Data(contentsOf: URL(fileURLWithPath: args.positionals[0]), options: [.mappedIfSafe]))
        let selection = try VivoSelectionQuery.parse(query)
        let conformer = try args.nonnegativeInt("conformer", default: 0)
        let result = try VivoSelectionEvaluator.evaluate(selection, in: document.structure,
                                                         conformerIndex: conformer,
                                                         periodic: !args.flags.contains("no-periodic"))
        let payload = SelectionReport(schemaVersion: 1, structureFingerprint: document.structureFingerprint,
                                      query: query, count: UInt32(result.atomIndices.count),
                                      atomIndices: result.atomIndices)
        try write(try VivoCanonicalJSON.encode(payload), to: args.options["output"] ?? "-", force: args.force)
        return 0
    }

    private func resolvedFormat(_ explicit: String?, path: String) throws -> VivoStructureFormat {
        if let explicit {
            switch explicit.lowercased() {
            case "json", "vivojson", "vivo-json": return .vivoJSON
            case "pdb": return .pdb
            case "cif", "mmcif": return .mmcif
            case "sdf", "mol": return .sdf
            case "mol2": return .mol2
            case "smi", "smiles": return .smiles
            default: throw CLIError.usage("unknown structure format '\(explicit)'")
            }
        }
        guard let inferred = VivoStructureFormat.infer(path: path) else {
            throw CLIError.usage("cannot infer structure format from '\(path)'; provide --format")
        }
        return inferred
    }

    private func write(_ data: Data, to path: String, force: Bool) throws {
        if path == "-" {
            FileHandle.standardOutput.write(data)
            if data.last != 0x0a { FileHandle.standardOutput.write(Data([0x0a])) }
            return
        }
        if FileManager.default.fileExists(atPath: path), !force {
            throw CLIError.usage("output exists: \(path); use --force to replace it")
        }
        let url = URL(fileURLWithPath: path)
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        try data.write(to: url, options: .atomic)
    }

    private func diagnostic(_ message: String) {
        FileHandle.standardError.write(Data(("numivivo structure: " + message + "\n").utf8))
    }

    static let help = """
    Native molecular structure commands

      numivivo structure-import <input> [--format pdb|mmcif|sdf|mol2|smiles|json]
          [--identifier <id>] [--output <structure.json>] [--force]
      numivivo structure-export <structure.json> --output <path>
          [--format pdb|sdf|mol2|json] [--conformer 0] [--force]
      numivivo structure-inspect <structure.json> [--output <report.json>]
      numivivo structure-select <structure.json> --query <selection>
          [--conformer 0] [--no-periodic] [--output <selection.json>]

    Selection examples:
      element(C,N,O,S)
      and(protein,resseq(450,500),not(hydrogen(true)))
      within(0.35,resname(LIG),false)

    Coordinates are canonicalized to nanometres. Import never silently drops
    unsupported stereochemistry or substitutes a lossy file format.
    """ + "\n"

    private struct InspectReport: Encodable {
        let schemaVersion: UInt32
        let structureFingerprint: VivoFingerprint
        let documentFingerprint: VivoFingerprint
        let sourceFingerprint: VivoFingerprint?
        let sourceFormat: VivoStructureFormat?
        let identifier: String
        let atomCount: UInt32
        let bondCount: UInt32
        let residueCount: UInt32
        let chainCount: UInt32
        let conformerCount: UInt32
        let connectedComponents: UInt32
        let hasPeriodicCell: Bool
        let elements: [String: UInt32]
    }

    private struct SelectionReport: Encodable {
        let schemaVersion: UInt32
        let structureFingerprint: VivoFingerprint
        let query: String
        let count: UInt32
        let atomIndices: [UInt32]
    }

    private enum CLIError: Error, CustomStringConvertible {
        case usage(String)
        var description: String { switch self { case .usage(let message): message } }
    }

    private struct Arguments {
        var positionals: [String] = []
        var options: [String: String] = [:]
        var flags = Set<String>()
        var force: Bool { flags.contains("force") }

        init(_ raw: [String]) throws {
            var i = 0
            while i < raw.count {
                let token = raw[i]
                if token.hasPrefix("--") {
                    let key = String(token.dropFirst(2))
                    if key == "force" || key == "no-periodic" {
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
        func required(_ key: String) throws -> String {
            guard let value = options[key], !value.isEmpty else { throw CLIError.usage("missing --\(key)") }
            return value
        }
        func allow(_ allowed: Set<String>) throws {
            let supplied = Set(options.keys).union(flags)
            guard supplied.isSubset(of: allowed) else {
                throw CLIError.usage("unknown option(s): \(supplied.subtracting(allowed).sorted().joined(separator: ", "))")
            }
        }
        func nonnegativeInt(_ key: String, default fallback: Int) throws -> Int {
            guard let raw = options[key] else { return fallback }
            guard let value = Int(raw), value >= 0 else { throw CLIError.usage("--\(key) requires a nonnegative integer") }
            return value
        }
    }
}
