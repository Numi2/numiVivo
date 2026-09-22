import Foundation

/// Reuses the MD owner's ordered connectivity rather than choosing independent
/// periodic images for two pull groups. A closure mismatch is a winding molecule,
/// not permission to pick a different COM image.
public struct VivoProteinWholeMoleculeLayout: Codable, Sendable, Equatable {
    public let orderedMembers: [UInt32]
    public let parents: [UInt32]
    public let closureEdges: [[UInt32]]
    public init(orderedMembers: [UInt32], parents: [UInt32], closureEdges: [[UInt32]]) {
        self.orderedMembers = orderedMembers; self.parents = parents; self.closureEdges = closureEdges
    }
    public func validate(particleCount: Int) throws {
        guard particleCount > 0, orderedMembers.count == parents.count, !orderedMembers.isEmpty,
              orderedMembers.count <= particleCount, parents[0] == .max else {
            throw VivoProteinStressError.invalid("whole-molecule traversal shape")
        }
        var visited = Set<UInt32>()
        for i in orderedMembers.indices {
            let atom = orderedMembers[i]
            guard Int(atom) < particleCount, !visited.contains(atom),
                  i == 0 || visited.contains(parents[i]) else {
                throw VivoProteinStressError.invalid("whole-molecule traversal is not a rooted tree")
            }
            visited.insert(atom)
        }
        guard closureEdges.count <= 4_000_000,
              closureEdges.allSatisfy({ $0.count == 2 && $0[0] != $0[1] && visited.contains($0[0]) && visited.contains($0[1]) }) else {
            throw VivoProteinStressError.invalid("whole-molecule closure edges")
        }
    }
    public func reconstruct(positionsNM: [VivoVector3D], toleranceNM: Double = 1e-6,
                            minimumImage: (VivoVector3D) throws -> VivoVector3D) throws -> [VivoVector3D] {
        try validate(particleCount: positionsNM.count)
        guard positionsNM.allSatisfy(\.isFinite), toleranceNM.isFinite, toleranceNM > 0 else {
            throw VivoProteinStressError.invalid("whole-molecule coordinates or tolerance")
        }
        var whole = positionsNM
        for i in orderedMembers.indices.dropFirst() {
            let child = Int(orderedMembers[i]), parent = Int(parents[i])
            let delta = try minimumImage(positionsNM[child] - positionsNM[parent])
            whole[child] = whole[parent] + delta
            guard delta.isFinite, whole[child].isFinite else { throw VivoProteinStressError.invalid("nonfinite molecular image") }
        }
        for edge in closureEdges {
            let a = Int(edge[0]), b = Int(edge[1])
            let displacement = try minimumImage(positionsNM[b] - positionsNM[a])
            let residual = ((whole[b] - whole[a]) - displacement).norm
            guard residual.isFinite, residual <= toleranceNM else {
                throw VivoProteinStressError.invalid("periodically winding or inconsistent molecular connectivity")
            }
        }
        return whole
    }
}
