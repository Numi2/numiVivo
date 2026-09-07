import Foundation

public enum VivoMechanismBondState: String, Codable, Sendable, CaseIterable {
    case absent, single, double, triple, aromatic, dative
    public init(_ order: VivoBondOrder?) {
        switch order {
        case nil: self = .absent
        case .single?: self = .single
        case .double?: self = .double
        case .triple?: self = .triple
        case .aromatic?: self = .aromatic
        case .dative?: self = .dative
        case .unknown?: self = .absent
        }
    }
    public var order: VivoBondOrder? {
        switch self {
        case .absent: return nil
        case .single: return .single
        case .double: return .double
        case .triple: return .triple
        case .aromatic: return .aromatic
        case .dative: return .dative
        }
    }
    public var orderValue: Double {
        switch self {
        case .absent: return 0
        case .single, .dative: return 1
        case .double: return 2
        case .triple: return 3
        case .aromatic: return 1.5
        }
    }
}

public struct VivoMechanismEditableBond: Codable, Sendable, Equatable {
    public var atomA: UInt32
    public var atomB: UInt32
    /// Allowed target states, excluding an implicit no-change option.
    public var allowedTargets: [VivoMechanismBondState]
    public init(atomA: UInt32, atomB: UInt32, allowedTargets: [VivoMechanismBondState]) {
        self.atomA = atomA; self.atomB = atomB; self.allowedTargets = allowedTargets
    }
}

public struct VivoMechanismAtomValenceLimit: Codable, Sendable, Equatable {
    public var atomIndex: UInt32
    public var maximumBondOrderSum: Double
    public init(atomIndex: UInt32, maximumBondOrderSum: Double) {
        self.atomIndex = atomIndex; self.maximumBondOrderSum = maximumBondOrderSum
    }
}

public struct VivoMechanismCandidateRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mechanism-candidate-request/v1"
    public var schema: String
    public var identifier: String
    public var structure: VivoMolecularStructure
    public var editableBonds: [VivoMechanismEditableBond]
    public var valenceLimits: [VivoMechanismAtomValenceLimit]
    public var maximumEditsPerCandidate: Int
    public var maximumCandidates: Int
    public init(identifier: String, structure: VivoMolecularStructure,
                editableBonds: [VivoMechanismEditableBond],
                valenceLimits: [VivoMechanismAtomValenceLimit] = [],
                maximumEditsPerCandidate: Int = 2,
                maximumCandidates: Int = 10_000) {
        schema = Self.schema; self.identifier = identifier; self.structure = structure
        self.editableBonds = editableBonds; self.valenceLimits = valenceLimits
        self.maximumEditsPerCandidate = maximumEditsPerCandidate; self.maximumCandidates = maximumCandidates
    }
}

public struct VivoMechanismBondChange: Codable, Sendable, Equatable, Hashable {
    public let atomA: UInt32
    public let atomB: UInt32
    public let from: VivoMechanismBondState
    public let to: VivoMechanismBondState
}

public struct VivoMechanismCandidate: Codable, Sendable, Equatable {
    public let identifier: String
    public let changes: [VivoMechanismBondChange]
    public let productStructureFingerprint: VivoFingerprint
    public let editMagnitude: Double
    public let suggestedCoordinate: VivoQMMMReactionCoordinate?
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

public struct VivoMechanismCandidateResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mechanism-candidate-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let sourceStructureFingerprint: VivoFingerprint
    public let candidates: [VivoMechanismCandidate]
    public let rejectedByValence: Int
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Enumerates graph hypotheses only. It does not infer electron flow, proton
/// reservoirs, barriers, transition states, or chemical feasibility. Every
/// candidate must still pass electronic/path/saddle/connectivity qualification.
public enum VivoMechanismCandidateEnumeration {
    public static let interpretation = "Bounded mapped bond-change hypothesis enumeration from caller-declared editable atom pairs and allowed target bond states. Atom identities are preserved and explicit valence caps are enforced. Ranking uses only edit magnitude; emitted reaction coordinates are heuristic seeds, not mechanism evidence."

    private struct Pair: Hashable { let a: UInt32; let b: UInt32 }
    private static func pair(_ a: UInt32, _ b: UInt32) -> Pair { .init(a:min(a,b), b:max(a,b)) }

    public static func enumerate(_ request: VivoMechanismCandidateRequest) throws -> VivoMechanismCandidateResult {
        guard request.schema == VivoMechanismCandidateRequest.schema,
              !request.identifier.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              !request.editableBonds.isEmpty, request.editableBonds.count <= 64,
              (1...8).contains(request.maximumEditsPerCandidate),
              (1...100_000).contains(request.maximumCandidates) else {
            throw VivoChemistryError.invalid("mechanism-candidate request identity or capacity")
        }
        _ = try VivoStructureValidator.validate(request.structure)
        let atomCount = UInt32(request.structure.atoms.count)
        var editablePairs = Set<Pair>()
        for rule in request.editableBonds {
            guard rule.atomA < atomCount, rule.atomB < atomCount, rule.atomA != rule.atomB,
                  editablePairs.insert(pair(rule.atomA,rule.atomB)).inserted,
                  !rule.allowedTargets.isEmpty, Set(rule.allowedTargets).count == rule.allowedTargets.count else {
                throw VivoChemistryError.invalid("mechanism editable bond")
            }
        }
        var limits:[UInt32:Double]=[:]
        for limit in request.valenceLimits {
            guard limit.atomIndex < atomCount, limits[limit.atomIndex] == nil,
                  limit.maximumBondOrderSum.isFinite, limit.maximumBondOrderSum > 0, limit.maximumBondOrderSum <= 12 else {
                throw VivoChemistryError.invalid("mechanism atom valence limit")
            }
            limits[limit.atomIndex]=limit.maximumBondOrderSum
        }
        var current:[Pair:VivoMechanismBondState]=[:]
        for bond in request.structure.bonds {
            guard bond.order != .unknown else {
                if editablePairs.contains(pair(bond.atomA,bond.atomB)) {
                    throw VivoChemistryError.unsupported("editable mechanism bonds require explicit bond order")
                }
                continue
            }
            current[pair(bond.atomA,bond.atomB)] = VivoMechanismBondState(bond.order)
        }
        struct Option { let pair:Pair; let from:VivoMechanismBondState; let targets:[VivoMechanismBondState] }
        let options = request.editableBonds.map { rule -> Option in
            let p=pair(rule.atomA,rule.atomB), from=current[p] ?? .absent
            return .init(pair:p,from:from,targets:rule.allowedTargets.filter{$0 != from})
        }
        guard options.allSatisfy({ !$0.targets.isEmpty }) else {
            throw VivoChemistryError.invalid("editable bond has no actual target state different from the source")
        }
        var raw:[[VivoMechanismBondChange]]=[]
        var chosen:[VivoMechanismBondChange]=[]
        func walk(_ index:Int) throws {
            if raw.count > request.maximumCandidates { throw VivoChemistryError.resourceLimit("mechanism candidate enumeration capacity") }
            if index == options.count {
                if !chosen.isEmpty && chosen.count <= request.maximumEditsPerCandidate { raw.append(chosen) }
                return
            }
            try walk(index+1)
            guard chosen.count < request.maximumEditsPerCandidate else { return }
            let o=options[index]
            for target in o.targets {
                chosen.append(.init(atomA:o.pair.a,atomB:o.pair.b,from:o.from,to:target))
                try walk(index+1)
                chosen.removeLast()
            }
        }
        try walk(0)
        guard raw.count <= request.maximumCandidates else { throw VivoChemistryError.resourceLimit("mechanism candidate capacity") }

        func valence(_ bonds:[VivoMolecularBond]) throws -> [Double] {
            var sums=[Double](repeating:0,count:request.structure.atoms.count)
            for bond in bonds {
                let state=VivoMechanismBondState(bond.order)
                guard bond.order != .unknown else { continue }
                sums[Int(bond.atomA)]+=state.orderValue; sums[Int(bond.atomB)]+=state.orderValue
            }
            return sums
        }
        var accepted:[VivoMechanismCandidate]=[], rejected=0
        for changes in raw {
            var product=request.structure
            var byPair=Dictionary(uniqueKeysWithValues: product.bonds.map{(pair($0.atomA,$0.atomB),$0)})
            for change in changes {
                let p=pair(change.atomA,change.atomB)
                guard (byPair[p].map{VivoMechanismBondState($0.order)} ?? .absent) == change.from else {
                    throw VivoChemistryError.invalid("mechanism source bond state changed during enumeration")
                }
                if let order=change.to.order {
                    byPair[p]=VivoMolecularBond(atomA:p.a,atomB:p.b,order:order)
                } else { byPair.removeValue(forKey:p) }
            }
            product.bonds=byPair.values.sorted {
                let a=pair($0.atomA,$0.atomB),b=pair($1.atomA,$1.atomB)
                return a.a == b.a ? a.b < b.b : a.a < b.a
            }
            _ = try VivoStructureValidator.validate(product)
            let sums=try valence(product.bonds)
            if limits.contains(where:{ sums[Int($0.key)] > $0.value + 1e-12 }) { rejected += 1; continue }
            let magnitude=changes.reduce(0.0){$0+abs($1.to.orderValue-$1.from.orderValue)}
            let forming=changes.filter{$0.from == .absent && $0.to != .absent}
            let breaking=changes.filter{$0.from != .absent && $0.to == .absent}
            let coordinate:VivoQMMMReactionCoordinate?
            if changes.count == 1 {
                let c=changes[0]
                coordinate=.init(identifier:"mechanism-\(accepted.count)-seed",kind:.distance,atomIndices:[c.atomA,c.atomB])
            } else if forming.count == 1 && breaking.count == 1 {
                coordinate=.init(identifier:"mechanism-\(accepted.count)-seed",kind:.distanceDifference,
                    atomIndices:[forming[0].atomA,forming[0].atomB,breaking[0].atomA,breaking[0].atomB])
            } else { coordinate=nil }
            let productID=try VivoStructureCodec.fingerprint(product)
            struct Evidence:Codable { let schema:String;let source:VivoFingerprint;let changes:[VivoMechanismBondChange];let product:VivoFingerprint }
            let sourceID=try VivoStructureCodec.fingerprint(request.structure)
            let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(schema:"numivivo.org/mechanism-candidate-evidence/v1",source:sourceID,changes:changes,product:productID)))
            accepted.append(.init(identifier:"candidate-\(accepted.count)",changes:changes,
                productStructureFingerprint:productID,editMagnitude:magnitude,suggestedCoordinate:coordinate,
                interpretation:"graph hypothesis only; requires independent electronic/path/saddle/connectivity qualification",
                evidenceFingerprint:evidenceID))
        }
        accepted.sort {
            if $0.editMagnitude != $1.editMagnitude { return $0.editMagnitude < $1.editMagnitude }
            let a=$0.changes.map{"\($0.atomA)-\($0.atomB):\($0.to.rawValue)"}.joined(separator:";")
            let b=$1.changes.map{"\($0.atomA)-\($0.atomB):\($0.to.rawValue)"}.joined(separator:";")
            return a < b
        }
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        let sourceID=try VivoStructureCodec.fingerprint(request.structure)
        struct ResultEvidence:Codable { let schema:String;let request:VivoFingerprint;let candidates:[VivoMechanismCandidate];let rejected:Int }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(ResultEvidence(schema:"numivivo.org/mechanism-candidate-result-evidence/v1",request:requestID,candidates:accepted,rejected:rejected)))
        return .init(schema:VivoMechanismCandidateResult.schema,requestFingerprint:requestID,
            sourceStructureFingerprint:sourceID,candidates:accepted,rejectedByValence:rejected,
            interpretation:interpretation,evidenceFingerprint:evidenceID)
    }
}
