import Foundation

public struct VivoChemicalStateVariant: Codable, Sendable, Equatable {
    public var identifier: String
    public var boundProtonOffset: Int
    public var removeHydrogens: [UInt32]
    public var addHydrogens: [VivoHydrogenPlacement]
    public var chargeEdits: [VivoPreparationChargeEdit]
    public var residueEdits: [VivoPreparationResidueEdit]
    public var origin: VivoKineticOrigin
    public var evidence: VivoKineticEvidence

    public init(identifier: String, boundProtonOffset: Int,
                removeHydrogens: [UInt32] = [], addHydrogens: [VivoHydrogenPlacement] = [],
                chargeEdits: [VivoPreparationChargeEdit] = [], residueEdits: [VivoPreparationResidueEdit] = [],
                origin: VivoKineticOrigin, evidence: VivoKineticEvidence) {
        self.identifier = identifier; self.boundProtonOffset = boundProtonOffset
        self.removeHydrogens = removeHydrogens; self.addHydrogens = addHydrogens
        self.chargeEdits = chargeEdits; self.residueEdits = residueEdits
        self.origin = origin; self.evidence = evidence
    }
    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 256, (-32...32).contains(boundProtonOffset),
              Set(removeHydrogens).count == removeHydrogens.count,
              Set(addHydrogens.map(\.identifier)).count == addHydrogens.count,
              Set(chargeEdits.map(\.sourceAtom)).count == chargeEdits.count,
              Set(residueEdits.map(\.sourceResidue)).count == residueEdits.count else {
            throw VivoChemistryError.invalid("chemical-state enumeration variant")
        }
        try evidence.validate(origin: origin)
    }
}

public struct VivoChemicalStateSite: Codable, Sendable, Equatable {
    public var identifier: String
    public var variants: [VivoChemicalStateVariant]
    public init(identifier: String, variants: [VivoChemicalStateVariant]) {
        self.identifier = identifier; self.variants = variants
    }
    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 256, variants.count >= 2, variants.count <= 64,
              Set(variants.map(\.identifier)).count == variants.count else {
            throw VivoChemistryError.invalid("chemical-state enumeration site")
        }
        for variant in variants { try variant.validate() }
    }
}

public struct VivoChemicalStateEnumerationRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-enumeration/v1"
    public var schema: String
    public var identifier: String
    public var structure: VivoMolecularStructure
    public var protonationSourceIdentifier: String
    public var referencePH: Double
    public var sites: [VivoChemicalStateSite]
    public var maximumStates: Int
    public var forceField: VivoForceFieldLibrary?
    public var compilationOptions: VivoForceFieldCompilationOptions?

    public init(identifier: String, structure: VivoMolecularStructure,
                protonationSourceIdentifier: String, referencePH: Double,
                sites: [VivoChemicalStateSite], maximumStates: Int = 4096,
                forceField: VivoForceFieldLibrary? = nil,
                compilationOptions: VivoForceFieldCompilationOptions? = nil) {
        schema = Self.schema; self.identifier = identifier; self.structure = structure
        self.protonationSourceIdentifier = protonationSourceIdentifier; self.referencePH = referencePH
        self.sites = sites; self.maximumStates = maximumStates
        self.forceField = forceField; self.compilationOptions = compilationOptions
    }
}

public struct VivoChemicalStateEnumerationResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-enumeration-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let catalogRequest: VivoChemicalStateCatalogRequest
    public let catalog: VivoChemicalStateCatalogResult
    public let combinationCount: Int
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Cartesian enumeration of explicitly supplied local chemical alternatives.
/// Adjacency changes exactly one local site. This prevents hand-written catalogs
/// from accidentally omitting combinations, but does not discover the local rules.
public enum VivoChemicalStateEnumeration {
    public static let interpretation = "Bounded combinatorial enumeration of explicit local protonation/tautomer state rules followed by authoritative molecular preparation and catalog validation. Generated states differ by supplied edits only. Rule discovery, pKa prediction, energetic ranking and chemical completeness are not inferred."

    public static func calculate(_ request: VivoChemicalStateEnumerationRequest) throws -> VivoChemicalStateEnumerationResult {
        guard request.schema == VivoChemicalStateEnumerationRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              !request.protonationSourceIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.referencePH.isFinite, (-10...30).contains(request.referencePH),
              !request.sites.isEmpty, request.sites.count <= 64,
              Set(request.sites.map(\.identifier)).count == request.sites.count,
              (1...100_000).contains(request.maximumStates) else {
            throw VivoChemistryError.invalid("chemical-state enumeration request")
        }
        _ = try VivoStructureValidator.validate(request.structure)
        for site in request.sites { try site.validate() }
        var count = 1
        for site in request.sites {
            let product = count.multipliedReportingOverflow(by: site.variants.count)
            guard !product.overflow, product.partialValue <= request.maximumStates else {
                throw VivoChemistryError.resourceLimit("chemical-state Cartesian product exceeds declared maximumStates")
            }
            count = product.partialValue
        }
        struct Choice {
            var variantIndices: [Int]
        }
        var choices: [Choice] = [.init(variantIndices: [])]
        for site in request.sites {
            var expanded: [Choice] = []
            expanded.reserveCapacity(choices.count * site.variants.count)
            for choice in choices {
                for index in site.variants.indices {
                    var next = choice.variantIndices; next.append(index)
                    expanded.append(.init(variantIndices: next))
                }
            }
            choices = expanded
        }
        let baseCharge = request.structure.atoms.reduce(0) { $0 + Int($1.formalCharge) }
        func stateIdentifier(_ choice: Choice) -> String {
            zip(request.sites, choice.variantIndices).map { "\($0.0.identifier)=\($0.0.variants[$0.1].identifier)" }.joined(separator: ";")
        }
        var entries: [VivoChemicalStateCatalogEntry] = []
        entries.reserveCapacity(choices.count)
        for choice in choices {
            let selected = zip(request.sites, choice.variantIndices).map { $0.0.variants[$0.1] }
            let removed = selected.flatMap(\.removeHydrogens)
            let added = selected.flatMap(\.addHydrogens)
            let charges = selected.flatMap(\.chargeEdits)
            let residues = selected.flatMap(\.residueEdits)
            guard Set(removed).count == removed.count,
                  Set(added.map(\.identifier)).count == added.count,
                  Set(charges.map(\.sourceAtom)).count == charges.count,
                  Set(residues.map(\.sourceResidue)).count == residues.count else {
                throw VivoChemistryError.invalid("chemical-state local rules conflict when combined")
            }
            let protonOffset = selected.reduce(0) { $0 + $1.boundProtonOffset }
            guard (-128...128).contains(protonOffset) else { throw VivoChemistryError.invalid("combined proton offset") }
            let identifier = stateIdentifier(choice)
            let preparation = VivoMolecularPreparationRequest(structure: request.structure,
                microstateIdentifier: identifier,
                protonationSourceIdentifier: request.protonationSourceIdentifier,
                pH: request.referencePH, expectedFormalCharge: baseCharge + protonOffset,
                removeHydrogens: removed, addHydrogens: added,
                chargeEdits: charges, residueEdits: residues,
                forceField: request.forceField, compilationOptions: request.compilationOptions)
            var neighbors: [String] = []
            for siteIndex in request.sites.indices {
                for alternate in request.sites[siteIndex].variants.indices where alternate != choice.variantIndices[siteIndex] {
                    var other = choice; other.variantIndices[siteIndex] = alternate
                    neighbors.append(stateIdentifier(other))
                }
            }
            struct SourceEvidence: Codable {
                let schema: String
                let enumerationIdentifier: String
                let stateIdentifier: String
                let selectedSiteVariants: [String]
                let sourceEvidence: [VivoKineticEvidence]
            }
            let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(SourceEvidence(
                schema: "numivivo.org/chemical-state-enumeration-state-evidence/v1",
                enumerationIdentifier: request.identifier, stateIdentifier: identifier,
                selectedSiteVariants: zip(request.sites, choice.variantIndices).map { "\($0.0.identifier)=\($0.0.variants[$0.1].identifier)" },
                sourceEvidence: selected.map(\.evidence))))
            let origin: VivoKineticOrigin = selected.contains(where: { $0.origin == .assumed }) ? .assumed : .calculated
            let evidence = VivoKineticEvidence(source: "NumiVivo explicit chemical-state rule enumeration",
                locator: identifier, sourceFingerprint: evidenceID.hex)
            entries.append(.init(identifier: identifier, boundProtonOffset: protonOffset,
                preparation: preparation, neighbors: neighbors.sorted(),
                stateOrigin: origin, stateEvidence: evidence))
        }
        entries.sort { $0.identifier < $1.identifier }
        let catalogRequest = VivoChemicalStateCatalogRequest(identifier: request.identifier, entries: entries)
        let catalog = try VivoChemicalStateCatalog.calculate(catalogRequest)
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let requestFingerprint: VivoFingerprint
            let catalogEvidenceFingerprint: VivoFingerprint
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/chemical-state-enumeration-evidence/v1",
            requestFingerprint: requestID, catalogEvidenceFingerprint: catalog.evidenceFingerprint)))
        return .init(schema: VivoChemicalStateEnumerationResult.schema,
            requestFingerprint: requestID, catalogRequest: catalogRequest, catalog: catalog,
            combinationCount: count, interpretation: interpretation, evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoChemicalStateEnumerationResult,
                                request: VivoChemicalStateEnumerationRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoChemistryError.invalid("chemical-state enumeration does not reconstruct")
        }
    }
}