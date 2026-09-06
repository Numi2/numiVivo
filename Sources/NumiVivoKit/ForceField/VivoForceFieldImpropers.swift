import Foundation

/// Periodic improper angle evaluated in exactly this atom order. The center is
/// separately declared because force fields do not share one ordering convention.
public struct VivoForceFieldTemplateImproper: Codable, Sendable, Equatable, Hashable {
    public var orderedAtoms: [String]
    public var center: String
    public init(orderedAtoms: [String], center: String) { self.orderedAtoms=orderedAtoms;self.center=center }
}
public struct VivoExplicitImproperTorsion: Codable, Sendable, Equatable, Hashable {
    public var orderedAtoms: [UInt32]
    public var center: UInt32
    public init(orderedAtoms: [UInt32], center: UInt32) { self.orderedAtoms=orderedAtoms;self.center=center }
}

extension VivoForceFieldCompiler {
    /// One matching authority for existing periodic torsions. Prefer the most
    /// specific family, retain all its Fourier terms, reject equally specific
    /// competing families instead of summing unrelated wildcard patterns.
    static func selectedTorsions(_ types: [String], improper: Bool,
                                 library: VivoForceFieldLibrary) throws -> [VivoForceFieldTorsionParameter] {
        var families: [[String]: [VivoForceFieldTorsionParameter]] = [:]
        var best = -1
        for parameter in library.torsionParameters where parameter.improper == improper {
            let pattern=[parameter.typeA,parameter.typeB,parameter.typeC,parameter.typeD]
            func matches(_ p:[String]) -> Bool {
                zip(p,types).allSatisfy { $0.0 == "*" || $0.0 == $0.1 }
            }
            // Full reversal preserves this dihedral; no other permutations are
            // allowed. In particular, peripheral atoms are never sorted.
            let reverse=Array(pattern.reversed())
            guard matches(pattern) || matches(reverse) else { continue }
            let specificity=pattern.filter{$0 != "*"}.count
            guard specificity>=best else { continue }
            if specificity>best { families.removeAll();best=specificity }
            let canonical=pattern.lexicographicallyPrecedes(reverse) ? pattern:reverse
            families[canonical,default:[]].append(parameter)
        }
        guard families.count<=1 else {
            throw VivoArtifactValidationError.incompatible("ambiguous equal-specificity torsion families for \(types)")
        }
        let terms=families.values.first ?? []
        // A repeated phase/periodicity in the SAME family is not another
        // physical torsion. Reject rather than silently double its coefficient.
        var seen=Set<String>()
        for term in terms {
            guard seen.insert("\(term.periodicity):\(term.phaseRadians.bitPattern)").inserted else {
                throw VivoArtifactValidationError.invalid("duplicate Fourier term in a matched torsion family")
            }
        }
        return terms
    }

    static func compileImpropers(structure: VivoMolecularStructure, topology: VivoStructureIndex,
                                 library: VivoForceFieldLibrary, particles: [VivoClassicalParticle],
                                 additional: [VivoExplicitImproperTorsion]) throws -> [VivoPeriodicTorsion] {
        let hasParameters=library.torsionParameters.contains(where: \.improper)
        var templates:[String:VivoForceFieldResidueTemplate]=[:]
        for template in library.residueTemplates { for alias in template.residueNames { templates[alias.uppercased()]=template } }
        var assignments=additional
        for residue in structure.residues {
            guard let template=templates[residue.name.uppercased()] else {
                if hasParameters { throw VivoArtifactValidationError.unresolved("improper topology requires a native template for residue \(residue.index)") }
                continue
            }
            guard let impropers=template.impropers else {
                if hasParameters {
                    throw VivoArtifactValidationError.unresolved("improper topology is unspecified for residue \(residue.index); provide ordered recipes or explicit []")
                }
                continue
            }
            let names=residue.atomIndices.map{structure.atoms[Int($0)].name}
            guard Set(names).count==names.count else { throw VivoArtifactValidationError.invalid("improper topology has duplicate residue atom names") }
            let byName=Dictionary(uniqueKeysWithValues:zip(names,residue.atomIndices))
            for recipe in impropers {
                let indices=recipe.orderedAtoms.compactMap{byName[$0]}
                guard indices.count==4,let center=byName[recipe.center] else {
                    throw VivoArtifactValidationError.unresolved("native improper recipe references a missing prepared atom")
                }
                assignments.append(.init(orderedAtoms:indices,center:center))
            }
        }
        var output:[VivoPeriodicTorsion]=[],seen=Set<[UInt32]>()
        var cache:[[String]:[VivoForceFieldTorsionParameter]]=[:]
        for recipe in assignments {
            let atoms=recipe.orderedAtoms
            guard atoms.count==4,Set(atoms).count==4,atoms.allSatisfy({Int($0)<particles.count}),
                  atoms.contains(recipe.center),atoms.filter({$0 != recipe.center}).allSatisfy({topology.bondedAtoms[Int(recipe.center)].contains($0)}) else {
                throw VivoArtifactValidationError.invalid("improper recipe is not four distinct atoms bonded to its declared center")
            }
            let reverse=Array(atoms.reversed()),canonical=atoms.lexicographicallyPrecedes(reverse) ? atoms:reverse
            guard seen.insert(canonical).inserted else { throw VivoArtifactValidationError.invalid("duplicate ordered improper recipe") }
            let types=atoms.map{particles[Int($0)].typeIdentifier}
            let terms:[VivoForceFieldTorsionParameter]
            if let stored=cache[types] { terms=stored }
            else { terms=try selectedTorsions(types,improper:true,library:library);cache[types]=terms }
            guard !terms.isEmpty else { throw VivoArtifactValidationError.unresolved("no parameter for declared improper \(types)") }
            for term in terms {
                output.append(.init(a:atoms[0],b:atoms[1],c:atoms[2],d:atoms[3],periodicity:term.periodicity,
                                    phaseRadians:term.phaseRadians,barrierKJPerMol:term.barrierKJPerMol,improper:true))
            }
        }
        return output
    }
}
