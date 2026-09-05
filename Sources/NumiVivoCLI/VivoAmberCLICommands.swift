import Foundation
import NumiVivoKit

struct VivoAmberCLICommands {
    static func handles(_ name: String?) -> Bool {
        name == "amber-import" || name == "amber-help"
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw CLIError.usage("missing AMBER command") }
            if command == "amber-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8)); return 0
            }
            let args = try Arguments(Array(arguments.dropFirst()))
            try args.requirePositionals(1)
            try args.allow(["restart", "identifier", "structure", "system", "state", "mapping", "force"])
            let prmtopPath = args.positionals[0]
            let prmtop = try Data(contentsOf: URL(fileURLWithPath: prmtopPath), options: [.mappedIfSafe])
            let restart = try args.options["restart"].map {
                try Data(contentsOf: URL(fileURLWithPath: $0), options: [.mappedIfSafe])
            }
            let identifier = args.options["identifier"] ?? URL(fileURLWithPath: prmtopPath).deletingPathExtension().lastPathComponent
            let imported = try VivoAmberImporter.importSystem(prmtopData: prmtop, restartData: restart, identifier: identifier)

            try write(try imported.structure.canonicalData(), to: try args.required("structure"), force: args.force)
            try write(try VivoCanonicalJSON.encode(imported.system), to: try args.required("system"), force: args.force)
            if let statePath = args.options["state"] {
                guard let state = imported.initialState else {
                    throw CLIError.usage("--state requested but no --restart was supplied")
                }
                try write(try VivoCanonicalJSON.encode(state), to: statePath, force: args.force)
            }
            if let mappingPath = args.options["mapping"] {
                let mapping = Mapping(schemaVersion: 1,
                                      particleToStructureAtom: imported.particleToStructureAtom,
                                      amberAtomTypes: imported.amberAtomTypes)
                try write(try VivoCanonicalJSON.encode(mapping), to: mappingPath, force: args.force)
            }
            let virtualSites = imported.particleToStructureAtom.filter { $0 == nil }.count
            diagnostic("imported particles=\(imported.system.particles.count) physicalAtoms=\(imported.structure.structure.atoms.count) virtualSites=\(virtualSites) bonds=\(imported.system.bonds.count) angles=\(imported.system.angles.count) torsions=\(imported.system.torsions.count)")
            return 0
        } catch let error as CLIError {
            diagnostic(error.description); return 64
        } catch let error as VivoArtifactValidationError {
            diagnostic(error.description); return 65
        } catch {
            diagnostic(error.localizedDescription); return 70
        }
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
        FileHandle.standardError.write(Data(("numivivo amber: " + message + "\n").utf8))
    }

    private struct Mapping: Encodable {
        let schemaVersion: UInt32
        let particleToStructureAtom: [UInt32?]
        let amberAtomTypes: [String]
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
            var i = 0
            while i < raw.count {
                let token = raw[i]
                if token.hasPrefix("--") {
                    let key = String(token.dropFirst(2))
                    if key == "force" {
                        guard flags.insert(key).inserted else { throw CLIError.usage("duplicate --force") }
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
            guard supplied.isSubset(of: allowed) else {
                throw CLIError.usage("unknown option(s): \(supplied.subtracting(allowed).sorted().joined(separator: ", "))")
            }
        }
        func required(_ key: String) throws -> String {
            guard let value = options[key], !value.isEmpty else { throw CLIError.usage("missing --\(key)") }
            return value
        }
    }

    static let help = """
    AMBER compatibility bridge

      numivivo amber-import <system.prmtop> --structure <structure.json> --system <classical-system.json>
          [--restart <system.rst7>] [--state <initial-state.json>]
          [--mapping <particle-map.json>] [--identifier <id>] [--force]

    The importer preserves AMBER bonded parameters, exclusions, per-dihedral
    1-4 SCEE/SCNB scaling, explicit Lennard-Jones type-pair coefficients, charges,
    and massless extra points. AMBER restart coordinates are converted Å -> nm.
    Restart velocities are detected and preserved as an explicit unconverted Wave-B
    boundary rather than being silently interpreted with the wrong units.
    """ + "\n"
}
