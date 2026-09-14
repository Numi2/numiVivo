import Foundation
import NumiVivoKit

/// Read-only assessment of a cross-scale evidence graph. The command reports
/// the contract state for every declared boundary; it never invents missing
/// links or upgrades hypotheses to biological outcomes.
struct VivoCrossScaleCLICommands {
    private static let maximumGraphBytes = 64 * 1_024 * 1_024

    static func handles(_ name: String?) -> Bool {
        ["cross-scale-assess", "cross-scale-help"].contains(name ?? "")
    }

    func run(arguments: [String]) -> Int32 {
        do {
            guard let command = arguments.first else { throw VivoOmicsError.invalid(Self.help) }
            switch command {
            case "cross-scale-help":
                guard arguments.count == 1 else { throw VivoOmicsError.invalid("cross-scale-help takes no arguments") }
                FileHandle.standardOutput.write(Data(Self.help.utf8))
                return 0
            case "cross-scale-assess":
                var input: String?
                var requireQualified = false
                for token in arguments.dropFirst() {
                    if token == "--require-qualified" {
                        guard !requireQualified else { throw VivoOmicsError.invalid("duplicate --require-qualified") }
                        requireQualified = true
                    } else {
                        guard input == nil, token != "-", !token.hasPrefix("--"), !token.isEmpty else {
                            throw VivoOmicsError.invalid("expected one graph JSON file and optional --require-qualified")
                        }
                        input = token
                    }
                }
                guard let input else {
                    throw VivoOmicsError.invalid("cross-scale-assess <graph.json> [--require-qualified]")
                }
                let rawURL = URL(fileURLWithPath: input).standardizedFileURL
                guard rawURL.isFileURL, rawURL.path.utf8.count <= 8_192,
                      !rawURL.path.contains("\0"), !rawURL.lastPathComponent.isEmpty else {
                    throw VivoOmicsError.invalid("cross-scale graph path must be a bounded local file")
                }
                // Resolve symlinks in existing parent directories (for example
                // macOS /tmp) while preserving the leaf for O_NOFOLLOW input
                // admission in VivoSingleCellCampaignIO.
                let url = rawURL.deletingLastPathComponent().resolvingSymlinksInPath()
                    .appendingPathComponent(rawURL.lastPathComponent).standardizedFileURL
                let bytes = try VivoSingleCellCampaignIO.readDocument(url, maximumBytes: Self.maximumGraphBytes)
                let graph = try VivoCanonicalJSON.decode(VivoCrossScaleEvidenceGraph.self, from: bytes)
                let assessment = try graph.assess()
                FileHandle.standardOutput.write(try VivoCanonicalJSON.encode(assessment))
                FileHandle.standardOutput.write(Data("\n".utf8))
                if requireQualified {
                    do {
                        try graph.requireQualifiedResearchPath()
                    } catch {
                        FileHandle.standardError.write(Data("numivivo cross-scale: \(Self.message(error))\n".utf8))
                        return 2
                    }
                }
                return 0
            default:
                throw VivoOmicsError.invalid(Self.help)
            }
        } catch {
            FileHandle.standardError.write(Data("numivivo cross-scale: \(Self.message(error))\n".utf8))
            return 65
        }
    }

    private static func message(_ error: Error) -> String {
        (error as? LocalizedError)?.errorDescription ?? String(describing: error)
    }

    static let help = """
    Cross-scale evidence assessment
      cross-scale-assess <graph.json> [--require-qualified]
      cross-scale-help

    The input is a bounded VivoCrossScaleEvidenceGraph JSON document. Output is
    a machine-readable assessment of the seven adjacent boundaries. Assessment
    succeeds for incomplete and hypothesis-only graphs so their limitations
    remain visible. --require-qualified returns exit code 2 unless every link
    has source, model, validation and held-out observation evidence. A
    qualified research path is not clinical, treatment or patient-specific
    authorization.
    """
}
