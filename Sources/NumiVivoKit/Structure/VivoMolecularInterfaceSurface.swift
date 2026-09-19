import Foundation

extension VivoMolecularInterface {
    /// Explicit heavy-atom Shrake–Rupley geometry. Radii are analysis inputs,
    /// not assigned force-field parameters; probe and all radii are in nm.
    public struct SurfacePlan: Codable, Sendable, Equatable {
        public let schemaVersion: Int
        public let radiusProfile: String
        public let radiiNM: [String: Double]
        public let probeRadiusNM: Double
        public let pointsPerAtom: Int
        public let overlapThresholdNM: Double
        public let maximumNeighborPairs: Int
        public let maximumPointTests: Int
        public let maximumNeighborSearchTests: Int
        public init(radiusProfile: String, radiiNM: [String: Double],
                    probeRadiusNM: Double = 0.14, pointsPerAtom: Int = 960,
                    overlapThresholdNM: Double = 0.04, maximumNeighborPairs: Int = 2_000_000,
                    maximumPointTests: Int = 100_000_000, maximumNeighborSearchTests: Int = 10_000_000) {
            schemaVersion = 1; self.radiusProfile = radiusProfile; self.radiiNM = radiiNM
            self.probeRadiusNM = probeRadiusNM; self.pointsPerAtom = pointsPerAtom
            self.overlapThresholdNM = overlapThresholdNM; self.maximumNeighborPairs = maximumNeighborPairs
            self.maximumPointTests = maximumPointTests; self.maximumNeighborSearchTests = maximumNeighborSearchTests
        }
    }
    public struct AtomSurface: Codable, Sendable, Equatable {
        public let atomIndex: UInt32
        public let isolatedAreaNM2: Double
        public let complexAreaNM2: Double
        public var buriedAreaNM2: Double { isolatedAreaNM2 - complexAreaNM2 }
    }
    public struct SurfaceReport: Codable, Sendable, Equatable {
        public let schemaVersion: Int
        public let method: String
        public let structureSHA256: String
        public let interfacePlan: Plan
        public let plan: SurfacePlan
        public let binderIsolatedAreaNM2: Double
        public let targetIsolatedAreaNM2: Double
        public let complexAreaNM2: Double
        public let binderBuriedAreaNM2: Double
        public let targetBuriedAreaNM2: Double
        /// SASA(A) + SASA(B) - SASA(A+B), NOT divided by two.
        public var buriedAreaSumNM2: Double { binderBuriedAreaNM2 + targetBuriedAreaNM2 }
        /// Half the summed loss, an explicitly named geometric area convention.
        public var halfBuriedAreaSumNM2: Double { buriedAreaSumNM2 / 2 }
        public let crossPartnerOverlapPairs: Int
        public let maximumCrossPartnerOverlapNM: Double
        public let neighborPairs: Int
        public let pointTests: Int
        public let neighborSearchTests: Int
        public let atoms: [AtomSurface]
        public let limitations: [String]
    }
    private struct SurfaceCell: Hashable {
        let x: Int; let y: Int; let z: Int
    }

    public static func analyzeSurface(_ structure: VivoMolecularStructure, interfacePlan: Plan,
                                       plan: SurfacePlan) throws -> SurfaceReport {
        // Reuse all structural, chain, occupancy, coordinate and topology guards.
        let interface = try analyze(structure, plan: interfacePlan)
        try surfaceRequire(plan.schemaVersion == 1 && !plan.radiusProfile.isEmpty
            && plan.radiusProfile.utf8.count <= 1024, "invalid surface schema/radius profile")
        try surfaceRequire((32...4096).contains(plan.pointsPerAtom)
            && plan.probeRadiusNM.isFinite && plan.probeRadiusNM > 0 && plan.probeRadiusNM <= 0.5
            && plan.overlapThresholdNM.isFinite && (0...0.5).contains(plan.overlapThresholdNM)
            && (1...5_000_000).contains(plan.maximumNeighborPairs)
            && (1...500_000_000).contains(plan.maximumPointTests)
            && (1...100_000_000).contains(plan.maximumNeighborSearchTests), "invalid surface resolution or budget")
        try surfaceRequire(!plan.radiiNM.isEmpty && plan.radiiNM.count <= 118, "invalid radius table")
        for (symbol, radius) in plan.radiiNM {
            try surfaceRequire(VivoElement.from(symbol: symbol)?.symbol == symbol
                && radius.isFinite && radius > 0 && radius <= 0.5, "invalid element/radius: \(symbol)")
        }
        let ids = (interface.binderHeavyAtomIndices + interface.targetHeavyAtomIndices).sorted()
        let left = Set(interface.binderHeavyAtomIndices)
        let conformer = structure.conformers.first { $0.identifier == interfacePlan.conformerID }!
        let positions = ids.map { conformer.positionsNM[Int($0)] }
        let isLeft = ids.map { left.contains($0) }
        let radii = try ids.map { id -> Double in
            let atom = structure.atoms[Int(id)]
            guard VivoElement.from(symbol: atom.element.symbol)?.atomicNumber == atom.element.atomicNumber,
                  let radius = plan.radiiNM[atom.element.symbol] else {
                throw VivoArtifactValidationError.invalid("surface requires an explicit radius for every selected element")
            }
            return radius
        }
        let expanded = radii.map { $0 + plan.probeRadiusNM }
        let width = 2 * expanded.max()!
        var cells: [SurfaceCell: [Int]] = [:]
        func cell(_ p: VivoVector3D) -> SurfaceCell {
            .init(x: Int(floor(p.x / width)), y: Int(floor(p.y / width)), z: Int(floor(p.z / width)))
        }
        for i in ids.indices { cells[cell(positions[i]), default: []].append(i) }
        var neighbors = Array(repeating: [Int](), count: ids.count)
        var pairs = 0, searchTests = 0, overlapPairs = 0, maximumOverlap = 0.0
        for i in ids.indices {
            let c = cell(positions[i])
            for dx in -1...1 { for dy in -1...1 { for dz in -1...1 {
                for j in cells[.init(x: c.x + dx, y: c.y + dy, z: c.z + dz)] ?? [] where j > i {
                    try surfaceRequire(searchTests < plan.maximumNeighborSearchTests, "surface neighbor-search budget exceeded")
                    searchTests += 1
                    let d2 = (positions[i] - positions[j]).squaredNorm
                    let limit = expanded[i] + expanded[j]
                    if d2 < limit * limit {
                        try surfaceRequire(d2 > 1e-24, "coincident heavy-atom centers make surface ownership ambiguous")
                        try surfaceRequire(pairs < plan.maximumNeighborPairs, "surface neighbor-pair budget exceeded")
                        neighbors[i].append(j); neighbors[j].append(i); pairs += 1
                        if isLeft[i] != isLeft[j] {
                            let overlap = radii[i] + radii[j] - sqrt(d2)
                            maximumOverlap = max(maximumOverlap, overlap)
                            if overlap > plan.overlapThresholdNM { overlapPairs += 1 }
                        }
                    }
                }
            } } }
        }
        // A deterministic golden-angle quadrature, like the reference family of
        // Shrake–Rupley implementations. Finite quadrature is not rotation exact.
        let n = plan.pointsPerAtom, angle = Double.pi * (3 - sqrt(5.0))
        let directions = (0..<n).map { k -> VivoVector3D in
            let z = 1 - (2 * Double(k) + 1) / Double(n)
            let r = sqrt(max(0, 1 - z * z)), phi = Double(k) * angle
            return .init(r * cos(phi), r * sin(phi), z)
        }
        var result: [AtomSurface] = [], tests = 0
        var binderIsolated = 0.0, targetIsolated = 0.0, complex = 0.0, binderBuried = 0.0, targetBuried = 0.0
        for i in ids.indices {
            // Same-partner occlusion first: a buried point cannot contribute to
            // either isolated or complex area. Keep order deterministic.
            let same = neighbors[i].filter { isLeft[$0] == isLeft[i] }.sorted()
            let other = neighbors[i].filter { isLeft[$0] != isLeft[i] }.sorted()
            var isolatedCount = 0, complexCount = 0
            func occluded(_ delta: VivoVector3D, _ js: [Int]) throws -> Bool {
                for j in js {
                    try surfaceRequire(tests < plan.maximumPointTests, "surface point-test budget exceeded")
                    tests += 1
                    // Relative displacements avoid adding/subtracting large offsets.
                    if (delta - (positions[j] - positions[i])).squaredNorm < expanded[j] * expanded[j] { return true }
                }
                return false
            }
            for direction in directions {
                let delta = direction * expanded[i]
                if try occluded(delta, same) { continue }
                isolatedCount += 1
                if try !occluded(delta, other) { complexCount += 1 }
            }
            let unit = 4 * Double.pi * expanded[i] * expanded[i] / Double(n)
            let isolated = Double(isolatedCount) * unit, bound = Double(complexCount) * unit
            result.append(.init(atomIndex: ids[i], isolatedAreaNM2: isolated, complexAreaNM2: bound))
            complex += bound
            if isLeft[i] { binderIsolated += isolated; binderBuried += isolated - bound }
            else { targetIsolated += isolated; targetBuried += isolated - bound }
        }
        return SurfaceReport(schemaVersion: 1, method: "numivivo-heavy-atom-shrake-rupley-v1",
            structureSHA256: interface.structureSHA256, interfacePlan: interfacePlan, plan: plan,
            binderIsolatedAreaNM2: binderIsolated, targetIsolatedAreaNM2: targetIsolated,
            complexAreaNM2: complex, binderBuriedAreaNM2: binderBuried, targetBuriedAreaNM2: targetBuried,
            crossPartnerOverlapPairs: overlapPairs, maximumCrossPartnerOverlapNM: maximumOverlap,
            neighborPairs: pairs, pointTests: tests, neighborSearchTests: searchTests, atoms: result,
            limitations: ["Static heavy-atom SASA and radius overlap, not energy, affinity, polarity or hydrogen-bond assignment.",
                "Only the selected partners occlude surface; solvent, hydrogens and unselected atoms are excluded.",
                "Radii and probe are explicit analysis choices, not an assigned force field.",
                "Finite sphere quadrature has orientation/resolution error; repeat at higher resolution for convergence.",
                "Buried area is separated into both partner losses; half their sum is only an area convention."])
    }
    private static func surfaceRequire(_ condition: Bool, _ reason: String) throws {
        if !condition { throw VivoArtifactValidationError.invalid(reason) }
    }
}
