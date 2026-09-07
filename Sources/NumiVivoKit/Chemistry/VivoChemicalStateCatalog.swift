import Foundation

public struct VivoChemicalStateCatalogEntry: Codable, Sendable, Equatable {
    public var identifier: String
    public var boundProtonOffset: Int
    public var preparation: VivoMolecularPreparationRequest
    public var neighbors: [String]
    public var stateOrigin: VivoKineticOrigin
    public var stateEvidence: VivoKineticEvidence

    public init(identifier: String, boundProtonOffset: Int,
                preparation: VivoMolecularPreparationRequest,
                neighbors: [String], stateOrigin: VivoKineticOrigin,
                stateEvidence: VivoKineticEvidence) {
        self.identifier = identifier
        self.boundProtonOffset = boundProtonOffset
        self.preparation = preparation
        self.neighbors = neighbors
        self.stateOrigin = stateOrigin
        self.stateEvidence = stateEvidence
    }
}

public struct VivoChemicalStateCatalogRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-catalog/v1"
    public var schema: String
    public var identifier: String
    public var entries: [VivoChemicalStateCatalogEntry]

    public init(identifier: String, entries: [VivoChemicalStateCatalogEntry]) {
        schema = Self.schema
        self.identifier = identifier
        self.entries = entries
    }
}

public struct VivoPreparedChemicalState: Codable, Sendable, Equatable {
    public let identifier: String
    public let boundProtonOffset: Int
    public let expectedFormalCharge: Int
    public let preparationRequestFingerprint: VivoFingerprint
    public let preparationResult: VivoMolecularPreparationResult
    public let preparedStructureFingerprint: VivoFingerprint
    public let hydrogenCount: Int
    public let neighbors: [String]
    public let stateOrigin: VivoKineticOrigin
    public let stateEvidence: VivoKineticEvidence
}

public struct VivoChemicalStateCatalogResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/chemical-state-catalog-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let sourceStructureFingerprint: VivoFingerprint
    public let states: [VivoPreparedChemicalState]
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Validates and materializes an explicitly supplied chemical-state graph using
/// the existing authoritative preparation transaction. It proves mapping and
/// proton-count consistency; it does not claim the catalog is chemically complete.
public enum VivoChemicalStateCatalog {
    public static let interpretation = "Preparation-bound explicit chemical-state catalog with shared source structure, reciprocal state graph, persistent heavy-atom source mapping and proton-count consistency. The catalog is supplied by the caller; state enumeration, pKa prediction and completeness are not inferred."

    public static func calculate(_ request: VivoChemicalStateCatalogRequest) throws -> VivoChemicalStateCatalogResult {
        guard request.schema == VivoChemicalStateCatalogRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              !request.entries.isEmpty, request.entries.count <= 4096,
              Set(request.entries.map(\.identifier)).count == request.entries.count else {
            throw VivoChemistryError.invalid("chemical-state catalog identity or capacity")
        }
        let identifiers = Set(request.entries.map(\.identifier))
        for entry in request.entries {
            guard !entry.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
                  entry.identifier == entry.preparation.microstateIdentifier,
                  (-128...128).contains(entry.boundProtonOffset),
                  !entry.neighbors.isEmpty, entry.neighbors.count <= 4096,
                  Set(entry.neighbors).count == entry.neighbors.count,
                  !entry.neighbors.contains(entry.identifier),
                  entry.neighbors.allSatisfy(identifiers.contains) else {
                throw VivoChemistryError.invalid("chemical-state catalog entry identity or graph")
            }
            try entry.stateEvidence.validate(origin: entry.stateOrigin)
        }
        for entry in request.entries {
            for neighbor in entry.neighbors {
                guard request.entries.first(where: { $0.identifier == neighbor })?.neighbors.contains(entry.identifier) == true else {
                    throw VivoChemistryError.invalid("chemical-state catalog graph must be reciprocal")
                }
            }
        }
        guard let first = request.entries.first else { throw VivoChemistryError.invalid("empty chemical-state catalog") }
        let sourceFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(first.preparation.structure))
        let sourceAtoms = first.preparation.structure.atoms
        let heavySource = sourceAtoms.filter { $0.element.atomicNumber > 1 }.map(\.index)
        let sourceHydrogens = sourceAtoms.filter { $0.element.atomicNumber == 1 }.count
        var prepared: [VivoPreparedChemicalState] = []
        prepared.reserveCapacity(request.entries.count)
        for entry in request.entries {
            guard try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(entry.preparation.structure)) == sourceFingerprint else {
                throw VivoChemistryError.invalid("chemical-state entries do not share one source molecular structure")
            }
            let result = try VivoMolecularPreparation.prepare(entry.preparation)
            let heavyMapped = heavySource.compactMap { sourceIndex -> UInt32? in
                guard Int(sourceIndex) < result.sourceToPrepared.count else { return nil }
                return result.sourceToPrepared[Int(sourceIndex)]
            }
            guard heavyMapped.count == heavySource.count,
                  heavyMapped.allSatisfy({ result.structure.atoms[Int($0)].element.atomicNumber > 1 }) else {
                throw VivoChemistryError.invalid("chemical-state preparation removed or remapped a heavy source atom")
            }
            let hydrogens = result.structure.atoms.filter { $0.element.atomicNumber == 1 }.count
            let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(entry.preparation))
            let structureID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result.structure))
            prepared.append(.init(identifier: entry.identifier,
                boundProtonOffset: entry.boundProtonOffset,
                expectedFormalCharge: entry.preparation.expectedFormalCharge,
                preparationRequestFingerprint: requestID,
                preparationResult: result,
                preparedStructureFingerprint: structureID,
                hydrogenCount: hydrogens,
                neighbors: entry.neighbors,
                stateOrigin: entry.stateOrigin,
                stateEvidence: entry.stateEvidence))
        }
        guard let baseline = prepared.first else { throw VivoChemistryError.invalid("empty prepared chemical-state catalog") }
        for state in prepared {
            let declaredDelta = state.boundProtonOffset - baseline.boundProtonOffset
            let actualDelta = state.hydrogenCount - baseline.hydrogenCount
            guard declaredDelta == actualDelta else {
                throw VivoChemistryError.invalid("declared bound-proton offset differs from prepared hydrogen-count difference for \(state.identifier)")
            }
            let chargeDelta = state.expectedFormalCharge - baseline.expectedFormalCharge
            guard chargeDelta == declaredDelta else {
                throw VivoChemistryError.invalid("formal-charge difference is inconsistent with proton-only state change for \(state.identifier)")
            }
        }
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoChemicalStateCatalogRequest
            let states: [VivoPreparedChemicalState]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/chemical-state-catalog-evidence/v1",
            request: request, states: prepared)))
        _ = sourceHydrogens // retained intentionally in source identity, not used as a hidden reference.
        return .init(schema: VivoChemicalStateCatalogResult.schema,
            requestFingerprint: requestID,
            sourceStructureFingerprint: sourceFingerprint,
            states: prepared,
            interpretation: interpretation,
            evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoChemicalStateCatalogResult,
                                request: VivoChemicalStateCatalogRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoChemistryError.invalid("chemical-state catalog does not reconstruct")
        }
    }
}