import Foundation
import NumiVivoKit

struct VivoAmberExecutableCLICommands {
    static func handles(_ name: String?) -> Bool {
        name == "amber-import-md" || name == "amber-md-help"
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw Error.usage("missing AMBER MD command") }
            if command == "amber-md-help" {
                FileHandle.standardOutput.write(Data(Self.help.utf8))
                return 0
            }
            let a = try Arguments(Array(arguments.dropFirst()))
            try a.requirePositionals(1, "amber-import-md <system.prmtop> --restart <system.rst7> --structure <structure.json> --system <system.json> --state <state.json>")
            try a.allow(["restart", "structure", "system", "state", "mapping", "identifier", "force"])
            let restartPath = try a.required("restart")
            let structurePath = try a.required("structure")
            let systemPath = try a.required("system")
            let statePath = try a.required("state")
            try distinctDestinations([structurePath, systemPath, statePath] + [a.options["mapping"]].compactMap { $0 })
            let result = try VivoAmberExecutableImporter.importSystem(
                prmtopData: try Data(contentsOf: URL(fileURLWithPath: a.positionals[0])),
                restartData: try Data(contentsOf: URL(fileURLWithPath: restartPath)),
                identifier: a.options["identifier"] ?? "amber-md"
            )
            guard let state = result.initialState else {
                throw Error.data("executable AMBER conversion did not produce an initial state")
            }
            try write(VivoCanonicalJSON.encode(result.structure), to: structurePath, force: a.force)
            try write(VivoCanonicalJSON.encode(result.system), to: systemPath, force: a.force)
            try write(VivoCanonicalJSON.encode(state), to: statePath, force: a.force)
            if let mapping = a.options["mapping"] {
                let value = ParticleMap(particleToStructureAtom: result.particleToStructureAtom,
                                        amberAtomTypes: result.amberAtomTypes,
                                        virtualSites: result.system.linearVirtualSites ?? [])
                try write(VivoCanonicalJSON.encode(value), to: mapping, force: a.force)
            }
            let virtualCount = result.system.linearVirtualSites?.count ?? 0
            FileHandle.standardError.write(Data(
                "AMBER MD: \(result.system.particles.count) particles, \(virtualCount) resolved virtual sites, system \(try result.system.fingerprint())\n".utf8
            ))
            return 0
        } catch let error as Error {
            FileHandle.standardError.write(Data("numivivo amber-md: \(error.description)\n".utf8))
            return error.exitCode
        } catch {
            FileHandle.standardError.write(Data("numivivo amber-md: \(error.localizedDescription)\n".utf8))
            return 65
        }
    }

    private struct ParticleMap: Codable {
        let particleToStructureAtom: [UInt32?]
        let amberAtomTypes: [String]
        let virtualSites: [VivoLinearVirtualSite]
    }

    private func distinctDestinations(_ paths: [String]) throws {
        let normalized = paths.map { URL(fileURLWithPath: $0).standardizedFileURL.path }
        guard Set(normalized).count == normalized.count else {
            throw Error.usage("AMBER MD output destinations must be distinct")
        }
    }

    private func write(_ data: Data, to path: String, force: Bool) throws {
        let url = URL(fileURLWithPath: path)
        if FileManager.default.fileExists(atPath: url.path), !force {
            throw Error.usage("output exists: \(path); use --force")
        }
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(),
                                                withIntermediateDirectories: true)
        try data.write(to: url, options: .atomic)
    }

    static let help = """
    NumiVivo executable AMBER conversion

      numivivo amber-import-md system.prmtop --restart system.rst7 \
          --structure structure.json --system classical-system.json \
          --state initial-state.json [--mapping particle-map.json]

    Unlike raw `amber-import`, this command requires restart geometry and resolves
    every massless site into an authoritative Wave B construction rule. Current
    automatic resolution is deliberately limited to one-extra-point O-H-H water
    residues that are linear on the molecular bisector within 1e-5 nm.
    """ + "\n"
}

private enum Error: Swift.Error, CustomStringConvertible {
    case usage(String)
    case data(String)
    var description: String { switch self { case .usage(let v), .data(let v): return v } }
    var exitCode: Int32 { switch self { case .usage: return 64; case .data: return 65 } }
}

private struct Arguments {
    let positionals: [String]
    let options: [String: String]
    let force: Bool
    init(_ raw: [String]) throws {
        var p: [String] = [], o: [String: String] = [:], f = false, i = 0
        while i < raw.count {
            let token = raw[i]
            if token == "--force" {
                guard !f else { throw Error.usage("duplicate --force") }
                f = true; i += 1; continue
            }
            if token.hasPrefix("--") {
                let key = String(token.dropFirst(2))
                guard o[key] == nil, i + 1 < raw.count, !raw[i + 1].hasPrefix("--") else {
                    throw Error.usage("missing value or duplicate option \(token)")
                }
                o[key] = raw[i + 1]; i += 2
            } else { p.append(token); i += 1 }
        }
        positionals = p; options = o; force = f
    }
    func requirePositionals(_ count: Int, _ usage: String) throws {
        guard positionals.count == count else { throw Error.usage(usage) }
    }
    func allow(_ keys: Set<String>) throws {
        guard options.keys.allSatisfy(keys.contains) else { throw Error.usage("unknown option") }
    }
    func required(_ key: String) throws -> String {
        guard let value = options[key], !value.isEmpty else { throw Error.usage("missing --\(key)") }
        return value
    }
}
