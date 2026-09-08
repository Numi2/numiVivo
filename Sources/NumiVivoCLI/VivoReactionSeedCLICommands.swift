import Foundation
import NumiVivoKit

/// Geometry-only preparation of a fresh reaction request from a verified
/// accepted sampling prefix. Numerical qualification is a separate workflow
/// operation whose input can be this publication's `request` output artifact.
struct VivoReactionSeedCLICommands {
    static func handles(_ name: String?) -> Bool {
        ["reaction-seed", "reaction-seed-verify", "reaction-encounter-request", "reaction-seed-help"].contains(name ?? "")
    }

    private struct Conditions: Decodable {
        struct Reservoir: Decodable {
            let componentIndex: Int
            let value: Double
            let unit: VivoConditionalEncounterReservoirUnit
        }
        let schema: String
        let identifier: String
        let reactantEndpointIdentifier: String
        let temperatureK: Double
        let taggedComponentIndex: Int
        let reservoirs: [Reservoir]
        let observationTimesSeconds: [Double]
        let includeReservoirSensitivities: Bool?
        let numerics: VivoKineticSensitivityConfiguration?
    }

    func run(arguments raw: [String]) async -> Int32 {
        do {
            try Task.checkCancellation()
            if raw == ["reaction-seed-help"] {
                try FileHandle.standardOutput.write(contentsOf: Data(Self.help.utf8)); return 0
            }
            guard let command = raw.first, Self.handles(command), command != "reaction-seed-help",
                  raw.count >= 2, !raw[1].hasPrefix("--") else {
                throw VivoChemistryError.invalid("reaction seed requires one request or receipt JSON file")
            }
            var options: [String: String] = [:], allPayloads = false, index = 2
            let allowed: Set<String> = command == "reaction-encounter-request" ? ["--conditions", "--output"] :
                ["--store", "--output", "--read-limits"]
            while index < raw.count {
                let key = raw[index]
                if key == "--verify-all-payloads" {
                    guard !allPayloads, command != "reaction-encounter-request" else { throw VivoChemistryError.invalid("duplicate or inapplicable reaction seed switch") }
                    allPayloads = true; index += 1
                } else {
                    guard allowed.contains(key), options[key] == nil, index + 1 < raw.count,
                          !raw[index + 1].isEmpty, !raw[index + 1].hasPrefix("--") else {
                        throw VivoChemistryError.invalid("unknown, duplicate or incomplete reaction seed option")
                    }
                    options[key] = raw[index + 1]; index += 2
                }
            }
            if command != "reaction-encounter-request", options["--store"] == nil {
                throw VivoChemistryError.invalid("reaction seed requires --store")
            }
            let root = try options["--store"].map { try VivoWorkflowCLIDocumentPaths.canonicalURL(URL(fileURLWithPath: $0)) }
            let input = URL(fileURLWithPath: raw[1])
            let limits = try VivoMolecularSamplingCLILimits.load(options["--read-limits"])
            var destination: VivoKineticsDocumentIO.PreparedOutput?
            if let path = options["--output"] {
                let output = try VivoWorkflowCLIDocumentPaths.canonicalURL(URL(fileURLWithPath: path))
                if let root {
                    let prefix = root.path == "/" ? "/" : root.path + "/"
                    guard output != root, !output.path.hasPrefix(prefix) else {
                        throw VivoChemistryError.invalid("reaction seed output aliases artifact storage")
                    }
                }
                let inputs = [input] + [options["--read-limits"], options["--conditions"]].compactMap { $0 }.map { URL(fileURLWithPath: $0) }
                for source in inputs {
                    guard try !VivoWorkflowCLIDocumentPaths.aliases(output, source) else {
                        throw VivoChemistryError.invalid("reaction seed output aliases an input")
                    }
                }
                guard !FileManager.default.fileExists(atPath: output.path) else {
                    throw VivoChemistryError.invalid("reaction seed output exists; outputs are no-clobber")
                }
                try Task.checkCancellation()
                destination = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: output)
            }
            try Task.checkCancellation()
            let validation: VivoMolecularSamplingSelectionValidation? = allPayloads ? .allPayloads : nil
            let data: Data
            if command == "reaction-encounter-request" {
                guard let conditionsPath = options["--conditions"] else {
                    throw VivoChemistryError.invalid("encounter request requires --conditions")
                }
                let conditions = try VivoKineticsDocumentIO.read(Conditions.self, from: URL(fileURLWithPath: conditionsPath))
                guard conditions.schema == "numivivo.org/conditional-encounter-conditions/v1" else {
                    throw VivoChemistryError.invalid("encounter conditions schema")
                }
                let wrapper = try VivoKineticsDocumentIO.read(VivoReactionCalculationResult.self, from: input)
                let source: VivoTransitionStateTheoryResult
                switch wrapper {
                case .connectedReaction(let value), .transitionStateTheory(let value): source = value
                default: throw VivoChemistryError.invalid("encounter request needs a complete TST reaction result")
                }
                let components = try VivoConditionalEncounter.componentBindings(in: source)
                func component(_ index: Int) throws -> VivoConditionalEncounterComponent {
                    guard components.indices.contains(index) else { throw VivoChemistryError.invalid("encounter component index") }
                    return components[index]
                }
                let request = try VivoConditionalEncounterRequest(identifier: conditions.identifier,
                    reactantEndpointIdentifier: conditions.reactantEndpointIdentifier, temperatureK: conditions.temperatureK,
                    taggedComponent: component(conditions.taggedComponentIndex), reservoirs: conditions.reservoirs.map {
                        try .init(component: component($0.componentIndex), value: $0.value, unit: $0.unit)
                    }, observationTimesSeconds: conditions.observationTimesSeconds,
                    includeReservoirSensitivities: conditions.includeReservoirSensitivities ?? false,
                    numerics: conditions.numerics ?? .init())
                try VivoConditionalEncounter.validateConditions(request, source: source)
                data = try VivoCanonicalJSON.encode(request)
            } else if command == "reaction-seed" {
                let store = try VivoArtifactStore(rootURL: root!, createIfNeeded: false)
                let implementation = try VivoWorkflowCLIImplementation.fingerprint()
                let request = try VivoKineticsDocumentIO.read(VivoSamplingReactionSeedRequest.self, from: input)
                let verified = try await VivoMolecularSamplingExporter.verifiedExport(request.sourceExport, store: store,
                    implementationFingerprint: implementation, limits: limits, minimumValidation: validation)
                let publication = try await VivoSamplingReactionSeed.publish(request, source: verified,
                    implementationFingerprint: implementation)
                data = try VivoCanonicalJSON.encode(publication.receipt)
            } else {
                let store = try VivoArtifactStore(rootURL: root!, createIfNeeded: false)
                let implementation = try VivoWorkflowCLIImplementation.fingerprint()
                let receipt = try VivoKineticsDocumentIO.read(VivoSamplingReactionSeedReceipt.self, from: input)
                let verification = try await VivoSamplingReactionSeed.verify(receipt, store: store,
                    implementationFingerprint: implementation, limits: limits, minimumValidation: validation)
                data = try VivoCanonicalJSON.encode(verification)
            }
            try Task.checkCancellation()
            if let destination { try destination.write(data) }
            else {
                try FileHandle.standardOutput.write(contentsOf: data)
                try FileHandle.standardOutput.write(contentsOf: Data("\n".utf8))
            }
            return 0
        } catch is CancellationError {
            try? FileHandle.standardError.write(contentsOf: Data("Reaction seed cancelled.\n".utf8)); return 130
        } catch {
            try? FileHandle.standardError.write(contentsOf: Data("Reaction seed failed: \(error)\n".utf8)); return 65
        }
    }

    static let help = """
    Create initial reaction geometry from a verified accepted sampling export:
      numivivo reaction-seed request.json --store STORE [--output receipt.json]
      numivivo reaction-seed-verify receipt.json --store STORE [--output verified.json]
      Both commands accept [--read-limits limits.json] [--verify-all-payloads].
      numivivo reaction-encounter-request reaction-result.json --conditions conditions.json
          [--output encounter-request.json]
    A request declares the source export, destination reaction calculation,
    atom mappings and model-transfer statement. Only initial nuclear positions
    change. Masses, temperature, electronic model and reaction assumptions remain
    explicit destination inputs. Periodic geometry requires an explicit policy
    outside this interface. Partial sampling remains a geometry seed only.
    Feed the receipt's typed request output artifact to the existing
    vivo.native.reaction-qualification workflow operation for fresh qualification.
    The vivo.native.conditional-encounter operation then consumes a qualified
    reaction result and an explicit tagged-reactant/reservoir/time request.
    reaction-encounter-request binds explicitly selected component indices to
    their atom identities and qualified-point hashes. It checks conditions and
    creates a request; numerical reaction validation occurs when it is executed.
    Verification reads the immutable source chain without running chemistry or
    regenerating artifacts. CLI source receipts require this executable identity.
    File outputs are no-clobber and outside inputs and artifact storage.

    """
}
