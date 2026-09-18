import Foundation

/// Bounded, nonperiodic heavy-atom geometry on the existing molecular container.
/// Contact cutoffs are declared analysis choices, not energies or binding evidence.
public enum VivoMolecularInterface {
    public struct Plan: Codable, Sendable, Equatable {
        public let schemaVersion: Int
        public let binderChains: [String]
        public let targetChains: [String]
        public let conformerID: String
        public let contactDistanceNM: Double
        public let shortDistanceNM: Double
        public let maximumPairEvaluations: Int
        public init(binderChains: [String], targetChains: [String], conformerID: String,
                    contactDistanceNM: Double = 0.45, shortDistanceNM: Double = 0.20,
                    maximumPairEvaluations: Int = 25_000_000) {
            self.schemaVersion = 1; self.binderChains = binderChains; self.targetChains = targetChains
            self.conformerID = conformerID; self.contactDistanceNM = contactDistanceNM
            self.shortDistanceNM = shortDistanceNM; self.maximumPairEvaluations = maximumPairEvaluations
        }
    }
    public struct ResidueContact: Codable, Sendable, Equatable {
        public let binderResidueIndex: UInt32
        public let targetResidueIndex: UInt32
        public let minimumDistanceNM: Double
    }
    public struct Report: Codable, Sendable, Equatable {
        public let schemaVersion: Int
        public let method: String
        public let structureSHA256: String
        public let plan: Plan
        public let binderHeavyAtomIndices: [UInt32]
        public let targetHeavyAtomIndices: [UInt32]
        public let binderContactAtomIndices: [UInt32]
        public let targetContactAtomIndices: [UInt32]
        public let ignoredAtomCount: Int
        public let evaluatedAtomPairs: Int
        public let contactAtomPairs: Int
        public let shortDistanceAtomPairs: Int
        public let minimumDistanceNM: Double
        public let geometricCentroidDistanceNM: Double
        public let binderGeometricRadiusOfGyrationNM: Double
        public let residueContacts: [ResidueContact]
        public let limitations: [String]
    }

    public static func analyze(_ structure: VivoMolecularStructure, plan: Plan) throws -> Report {
        // Bound allocation before the general validator constructs topology indices.
        try require(structure.atoms.count <= 100_000 && structure.bonds.count <= 400_000
                    && structure.residues.count <= 100_000 && structure.chains.count <= 10_000
                    && structure.conformers.count <= 32, "interface input capacity exceeded")
        try VivoStructureValidator.validate(structure)
        try require(structure.periodicCell == nil, "interface geometry requires nonperiodic, explicitly unwrapped input")
        let binder = Set(plan.binderChains), target = Set(plan.targetChains)
        try require(plan.schemaVersion == 1 && !binder.isEmpty && !target.isEmpty
                    && binder.count == plan.binderChains.count && target.count == plan.targetChains.count
                    && binder.isDisjoint(with: target), "invalid or overlapping interface chain groups")
        try require(binder.union(target).isSubset(of: Set(structure.chains.map(\.identifier))), "unknown interface chain")
        try require(plan.contactDistanceNM.isFinite && plan.shortDistanceNM.isFinite
                    && plan.shortDistanceNM > 0 && plan.shortDistanceNM < plan.contactDistanceNM
                    && plan.contactDistanceNM <= 2, "require 0 < shortDistanceNM < contactDistanceNM <= 2")
        try require((1...100_000_000).contains(plan.maximumPairEvaluations), "invalid pair-evaluation budget")
        guard let conformer = structure.conformers.first(where: { $0.identifier == plan.conformerID }) else {
            throw VivoArtifactValidationError.invalid("unknown interface conformer")
        }
        let positions = conformer.positionsNM
        func group(_ atom: VivoMolecularAtom) -> String? {
            guard let r = atom.residueIndex, let c = structure.residues[Int(r)].chainIndex else { return nil }
            return structure.chains[Int(c)].identifier
        }
        var left: [UInt32] = [], right: [UInt32] = [], atomNames = Set<String>()
        var representedChains = Set<String>()
        for atom in structure.atoms {
            guard let chain = group(atom), binder.contains(chain) || target.contains(chain),
                  atom.element.atomicNumber != 1 else { continue }
            try require(atom.alternateLocation == nil || atom.alternateLocation == "",
                        "resolve alternate locations before interface analysis")
            try require(atom.occupancy == nil || atom.occupancy == 1,
                        "partial/zero occupancy is not supported by interface geometry")
            try require(atomNames.insert("\(atom.residueIndex!):\(atom.name)").inserted,
                        "duplicate residue/atom identity in selected interface")
            let p = positions[Int(atom.index)]
            try require(max(abs(p.x), abs(p.y), abs(p.z)) <= 1e6, "interface coordinate magnitude exceeds bound")
            representedChains.insert(chain)
            if binder.contains(chain) { left.append(atom.index) } else { right.append(atom.index) }
        }
        try require(representedChains == binder.union(target), "each selected chain requires heavy atoms")
        let work = left.count.multipliedReportingOverflow(by: right.count)
        try require(!work.overflow && work.partialValue <= plan.maximumPairEvaluations,
                    "interface pair-evaluation budget exceeded")
        let leftSet = Set(left), rightSet = Set(right)
        for bond in structure.bonds {
            try require(!((leftSet.contains(bond.atomA) && rightSet.contains(bond.atomB))
                        || (leftSet.contains(bond.atomB) && rightSet.contains(bond.atomA))),
                        "covalently linked partners require a different interface model")
        }
        let contact2 = plan.contactDistanceNM * plan.contactDistanceNM
        let short2 = plan.shortDistanceNM * plan.shortDistanceNM
        var minimum2 = Double.infinity, count = 0, shortCount = 0
        var leftContact = Set<UInt32>(), rightContact = Set<UInt32>()
        var residueMin: [UInt64: Double] = [:]
        for i in left {
            let pi = positions[Int(i)], ri = structure.atoms[Int(i)].residueIndex!
            for j in right {
                let d2 = (pi - positions[Int(j)]).squaredNorm
                minimum2 = min(minimum2, d2)
                if d2 < short2 { shortCount += 1 }
                if d2 <= contact2 {
                    count += 1; leftContact.insert(i); rightContact.insert(j)
                    let rj = structure.atoms[Int(j)].residueIndex!
                    let key = UInt64(ri) << 32 | UInt64(rj)
                    if let old = residueMin[key] { residueMin[key] = min(old, d2) }
                    else {
                        try require(residueMin.count < 100_000, "contact-residue-pair capacity exceeded")
                        residueMin[key] = d2
                    }
                }
            }
        }
        func centroid(_ indices: [UInt32]) -> VivoVector3D {
            // Difference from an anchor reduces cancellation under common translation.
            let anchor = positions[Int(indices[0])]
            return anchor + indices.reduce(.zero) { $0 + (positions[Int($1)] - anchor) } / Double(indices.count)
        }
        let center = centroid(left)
        let radius = sqrt(left.reduce(0.0) { $0 + (positions[Int($1)] - center).squaredNorm } / Double(left.count))
        let contacts = residueMin.keys.sorted().map { key in
            ResidueContact(binderResidueIndex: UInt32(key >> 32), targetResidueIndex: UInt32(key & 0xffffffff),
                           minimumDistanceNM: sqrt(residueMin[key]!))
        }
        return Report(schemaVersion: 1, method: "numivivo-heavy-atom-interface-v1",
            structureSHA256: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(structure)).hex,
            plan: plan, binderHeavyAtomIndices: left, targetHeavyAtomIndices: right,
            binderContactAtomIndices: leftContact.sorted(), targetContactAtomIndices: rightContact.sorted(),
            ignoredAtomCount: structure.atoms.count - left.count - right.count,
            evaluatedAtomPairs: work.partialValue, contactAtomPairs: count, shortDistanceAtomPairs: shortCount,
            minimumDistanceNM: sqrt(minimum2), geometricCentroidDistanceNM: (center - centroid(right)).norm,
            binderGeometricRadiusOfGyrationNM: radius, residueContacts: contacts,
            limitations: ["Geometry only; not binding affinity, hydrogen-bond assignment or a simulation.",
                "Contacts use <= contactDistanceNM; short-distance pairs use < shortDistanceNM.",
                "Short-distance pairs are not force-field or van-der-Waals clash classifications.",
                "Hydrogens and unselected atoms are excluded; geometric centroids are not mass weighted.",
                "No periodic images, alternative-location averaging or occupancy weighting.",
                "Declared cutoffs are analysis choices, not experimentally calibrated predictors."])
    }

    private static func require(_ condition: Bool, _ message: String) throws {
        if !condition { throw VivoArtifactValidationError.invalid(message) }
    }
}
