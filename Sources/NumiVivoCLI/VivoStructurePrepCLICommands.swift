import Foundation
import NumiVivoKit

struct VivoStructurePrepCLICommands {
    static func handles(_ name: String?) -> Bool {
        guard let name else { return false }
        return ["structure-resolve-altloc", "structure-build-topology", "structure-slice", "structure-prep-help"].contains(name)
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw CLIError.usage("missing structure preparation command") }
            if command == "structure-prep-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let args = try Arguments(Array(arguments.dropFirst()))
            switch command {
            case "structure-resolve-altloc": try resolveAltLoc(args)
            case "structure-build-topology": try buildTopology(args)
            case "structure-slice": try slice(args)
            default: throw CLIError.usage("unknown structure preparation command")
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

    private func resolveAltLoc(_ args: Arguments) throws {
        try args.requirePositionals(1)
        try args.allow(["preferred", "output", "mapping", "force"])
        let document = try load(args.positionals[0])
        let policy: VivoAlternateLocationPolicy = args.options["preferred"].map { .preferred($0) } ?? .highestOccupancy
        let rewrite = try VivoAlternateLocationResolver.resolve(document.structure, policy: policy)
        let output = try VivoMolecularStructureDocument(structure: rewrite.structure,
                                                        sourceFingerprint: document.sourceFingerprint,
                                                        sourceFormat: document.sourceFormat)
        try write(try output.canonicalData(), to: try args.required("output"), force: args.force)
        try writeMappingIfRequested(rewrite, path: args.options["mapping"], force: args.force)
        diagnostic("resolved alternate locations: \(document.structure.atoms.count) -> \(rewrite.structure.atoms.count) atoms")
    }

    private func buildTopology(_ args: Arguments) throws {
        try args.requirePositionals(1)
        try args.allow(["conformer", "no-hydrogens", "no-disulfides", "perceive-covalent",
                        "allow-inter-residue", "include-metals", "output", "force"])
        let document = try load(args.positionals[0])
        var structure = document.structure
        let conformer = try args.nonnegativeInt("conformer", default: 0)
        let report = try VivoBiopolymerTopologyBuilder.apply(
            to: &structure, conformerIndex: conformer,
            inferHydrogenAttachment: !args.flags.contains("no-hydrogens"),
            inferDisulfides: !args.flags.contains("no-disulfides")
        )
        var perceived: UInt32 = 0
        if args.flags.contains("perceive-covalent") {
            perceived = try VivoCovalentBondPerception.apply(
                to: &structure, conformerIndex: conformer,
                options: .init(allowInterResidue: args.flags.contains("allow-inter-residue"),
                               includeMetals: args.flags.contains("include-metals"))
            )
        }
        let output = try VivoMolecularStructureDocument(structure: structure,
                                                        sourceFingerprint: document.sourceFingerprint,
                                                        sourceFormat: document.sourceFormat)
        try write(try output.canonicalData(), to: try args.required("output"), force: args.force)
        diagnostic("topology bonds=\(structure.bonds.count) template=\(report.templateBondsAdded) peptide=\(report.peptideBondsAdded) disulfide=\(report.disulfideBondsAdded) hydrogen=\(report.hydrogenBondsAdded) perceived=\(perceived) unresolvedResidues=\(report.unresolvedResidues.count)")
    }

    private func slice(_ args: Arguments) throws {
        try args.requirePositionals(1)
        try args.allow(["query", "conformer", "no-periodic", "identifier", "output", "mapping", "force"])
        let document = try load(args.positionals[0])
        let query = try args.required("query")
        let selection = try VivoSelectionQuery.parse(query)
        let selected = try VivoSelectionEvaluator.evaluate(selection, in: document.structure,
                                                            conformerIndex: try args.nonnegativeInt("conformer", default: 0),
                                                            periodic: !args.flags.contains("no-periodic"))
        let rewrite = try VivoStructureSlicer.slice(document.structure, atomIndices: selected.atomIndices,
                                                    identifier: args.options["identifier"])
        let output = try VivoMolecularStructureDocument(structure: rewrite.structure,
                                                        sourceFingerprint: document.sourceFingerprint,
                                                        sourceFormat: document.sourceFormat)
        try write(try output.canonicalData(), to: try args.required("output"), force: args.force)
        try writeMappingIfRequested(rewrite, path: args.options["mapping"], force: args.force)
        diagnostic("slice atoms=\(rewrite.structure.atoms.count) bonds=\(rewrite.structure.bonds.count)")
    }

    private func load(_ path: String) throws -> VivoMolecularStructureDocument {
        try VivoMolecularStructureDocument.decode(Data(contentsOf: URL(fileURLWithPath: path), options: [.mappedIfSafe]))
    }

    private func writeMappingIfRequested(_ rewrite: VivoStructureRewrite, path: String?, force: Bool) throws {
        guard let path else { return }
        let mapping = MappingReport(schemaVersion: 1, oldToNew: rewrite.oldToNew, newToOld: rewrite.newToOld)
        try write(try VivoCanonicalJSON.encode(mapping), to: path, force: force)
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
        FileHandle.standardError.write(Data(("numivivo structure-prep: " + message + "\n").utf8))
    }

    private struct MappingReport: Encodable {
        let schemaVersion: UInt32
        let oldToNew: [UInt32?]
        let newToOld: [UInt32]
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
            let flagNames: Set<String> = ["force", "no-hydrogens", "no-disulfides", "perceive-covalent", "allow-inter-residue", "include-metals", "no-periodic"]
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
        func nonnegativeInt(_ key: String, default fallback: Int) throws -> Int {
            guard let raw = options[key] else { return fallback }
            guard let value = Int(raw), value >= 0 else { throw CLIError.usage("--\(key) requires a nonnegative integer") }
            return value
        }
    }

    static let help = """
    Native molecular structure preparation

      numivivo structure-resolve-altloc <structure.json> --output <resolved.json>
          [--preferred A] [--mapping <map.json>] [--force]
      numivivo structure-build-topology <structure.json> --output <topology.json>
          [--conformer 0] [--no-hydrogens] [--no-disulfides]
          [--perceive-covalent] [--allow-inter-residue] [--include-metals] [--force]
      numivivo structure-slice <structure.json> --query <selection> --output <slice.json>
          [--identifier <id>] [--mapping <map.json>] [--conformer 0] [--no-periodic] [--force]

    Standard amino-acid bonds and peptide links are deterministic. Geometry-based
    covalent perception is opt-in and conservative. Mapping files make every atom
    deletion/reindexing explicit for downstream force-field and QM state.
    """ + "\n"
}
