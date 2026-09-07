import Foundation

enum VivoNuclearMolecularGroups {
    static func build(system: VivoClassicalSystem, document: VivoMolecularStructureDocument? = nil) throws -> [[Int]] {
        try VivoClassicalSystemValidator.validate(system)
        let atoms = system.particles.filter { $0.role == .atom }, n = atoms.count
        let slot = Dictionary(uniqueKeysWithValues: atoms.enumerated().map { ($0.element.index,$0.offset) })
        var graph = [Set<Int>](repeating: [],count: n)
        func join(_ a: UInt32,_ b: UInt32) throws {
            guard let i = slot[a], let j = slot[b] else { throw VivoChemistryError.invalid("nuclear molecular connectivity must join physical atoms") }
            graph[i].insert(j); graph[j].insert(i)
        }
        for bond in system.bonds { try join(bond.a,bond.b) }
        for constraint in system.constraints { try join(constraint.a,constraint.b) }
        let sites = try system.resolvedVirtualSiteGraph()
        for parents in sites.physicalAncestors {
            guard let first = parents.first else { continue }
            for parent in parents.dropFirst() { try join(first,parent) }
        }
        if let document {
            try document.validate()
            guard document.structureFingerprint == system.structureFingerprint else { throw VivoChemistryError.invalid("nuclear source structure identity") }
            let atomMap = Dictionary(uniqueKeysWithValues: atoms.compactMap { p in p.atomIndex.map { ($0,p.index) } })
            for bond in document.structure.bonds {
                guard let a = atomMap[bond.atomA], let b = atomMap[bond.atomB] else { throw VivoChemistryError.invalid("nuclear source atom mapping") }
                try join(a,b)
            }
        }
        var visited = Set<Int>(), groups: [[Int]] = []
        for root in 0..<n where !visited.contains(root) {
            var queue = [root], cursor = 0; visited.insert(root)
            while cursor < queue.count {
                let i = queue[cursor]; cursor += 1
                for j in graph[i].sorted() where visited.insert(j).inserted { queue.append(j) }
            }
            groups.append(queue.sorted())
        }
        return groups
    }
}
