import Foundation

/// Converts a raw AMBER import into the stricter execution artifact required by
/// Wave B. Only virtual-site constructions that can be established from source
/// topology + restart geometry are admitted. Nothing is guessed from atom names.
public enum VivoAmberExecutableImporter {
    public static func importSystem(prmtopData: Data,
                                    restartData: Data? = nil,
                                    identifier: String = "amber") throws -> VivoAmberImportResult {
        var result = try VivoAmberImporter.importSystem(prmtopData: prmtopData,
                                                        restartData: restartData,
                                                        identifier: identifier)
        let virtualParticles = result.system.particles.filter { $0.role == .virtualSite }.map(\.index)
        let initialBeforeResolution = result.initialState

        if !virtualParticles.isEmpty {
            guard let initial = initialBeforeResolution else {
                throw VivoArtifactValidationError.unresolved(
                    "AMBER virtual sites require restart coordinates to establish their construction geometry"
                )
            }
            let top = try VivoAmberPrmtop(data: prmtopData)
            let atomicNumbers = try top.integers("ATOMIC_NUMBER")
            let residuePointers = try top.integers("RESIDUE_POINTER")
            let residueLabels = try top.strings("RESIDUE_LABEL")
            guard atomicNumbers.count == result.system.particles.count,
                  residuePointers.count == residueLabels.count,
                  residuePointers.first == 1 else {
                throw VivoArtifactValidationError.invalid(
                    "AMBER virtual-site resolver received inconsistent atom/residue arrays"
                )
            }

            let virtualSet = Set(virtualParticles)
            var resolved: [VivoLinearVirtualSite] = []
            var accounted = Set<UInt32>()
            for residue in residuePointers.indices {
                let begin = residuePointers[residue] - 1
                let end = residue + 1 < residuePointers.count
                    ? residuePointers[residue + 1] - 1
                    : atomicNumbers.count
                guard begin >= 0, end > begin, end <= atomicNumbers.count else {
                    throw VivoArtifactValidationError.invalid("AMBER residue particle range is invalid")
                }
                let indices = Array(begin..<end)
                let virtual = indices.filter { virtualSet.contains(UInt32($0)) }
                guard !virtual.isEmpty else { continue }

                // Wave B v1 recognizes the common four-site water family: one
                // massless point and physical O/H/H parents on the molecular bisector.
                let physical = indices.filter { atomicNumbers[$0] > 0 }
                let oxygen = physical.filter { atomicNumbers[$0] == 8 }
                let hydrogens = physical.filter { atomicNumbers[$0] == 1 }
                guard virtual.count == 1, physical.count == 3,
                      oxygen.count == 1, hydrogens.count == 2 else {
                    throw VivoArtifactValidationError.incompatible(
                        "AMBER residue \(residueLabels[residue]) contains virtual-site geometry outside the supported linear O-H-H water contract"
                    )
                }
                let site = virtual[0], o = oxygen[0], h1 = hydrogens[0], h2 = hydrogens[1]
                let positions = initial.positionsNM
                let oh1 = minimumImage(positions[h1] - positions[o], cell: initial.periodicCell)
                let oh2 = minimumImage(positions[h2] - positions[o], cell: initial.periodicCell)
                let os = minimumImage(positions[site] - positions[o], cell: initial.periodicCell)
                let bisector = (oh1 + oh2) * 0.5
                let denominator = bisector.squaredNorm
                guard denominator.isFinite, denominator > 1e-12 else {
                    throw VivoArtifactValidationError.invalid("AMBER water bisector is degenerate")
                }
                let alpha = os.dot(bisector) / denominator
                let residual = (bisector * alpha - os).norm
                guard alpha.isFinite, alpha > -0.25, alpha < 1.25,
                      residual.isFinite, residual <= 1e-5 else {
                    throw VivoArtifactValidationError.incompatible(
                        "AMBER water extra point is not a linear O-H-H bisector site within 1e-5 nm (residual \(residual))"
                    )
                }
                resolved.append(.init(
                    siteParticle: UInt32(site),
                    parentParticles: [UInt32(o), UInt32(h1), UInt32(h2)],
                    weights: [1 - alpha, alpha * 0.5, alpha * 0.5],
                    provenance: "derived from AMBER residue/restart geometry; linear O-H-H water extra point"
                ))
                accounted.insert(UInt32(site))
            }
            guard accounted == virtualSet else {
                let missing = virtualSet.subtracting(accounted).sorted()
                throw VivoArtifactValidationError.incompatible("AMBER virtual particles remain unresolved: \(missing)")
            }
            result.system.linearVirtualSites = resolved.sorted { $0.siteParticle < $1.siteParticle }
            result.system.metadata["virtualSiteResolution"] = "linear-water-from-restart"
            result.system.metadata["virtualSiteCount"] = String(resolved.count)
        }

        // The executable AMBER route targets the conventional 2 fs biomolecular
        // contract: constrain physical H bonds and close three-site water geometry.
        let constrained = try VivoClassicalConstraintCompiler.compile(
            system: result.system,
            structure: result.structure.structure,
            policy: .hydrogenBondsAndRigidWater
        )
        result.system = constrained.system
        result.system.metadata["addedHydrogenBondConstraints"] = String(constrained.report.addedHydrogenBondCount)
        result.system.metadata["addedWaterHHConstraints"] = String(constrained.report.addedWaterHHCount)

        try VivoClassicalSystemValidator.validate(result.system,
                                                   atomCount: UInt32(result.structure.structure.atoms.count))
        let reboundFingerprint = try result.system.fingerprint()
        if let initial = initialBeforeResolution {
            result.initialState = VivoClassicalInitialState(
                systemFingerprint: reboundFingerprint,
                positionsNM: initial.positionsNM,
                periodicCell: initial.periodicCell,
                sourceTimePS: initial.sourceTimePS
            )
            try result.initialState?.validate(particleCount: result.system.particles.count)
        }
        return result
    }

    private static func minimumImage(_ delta: VivoVector3D,
                                     cell: VivoPeriodicCell?) -> VivoVector3D {
        guard let cell else { return delta }
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else { return delta }
        let ra = cell.b.cross(cell.c) / determinant
        let rb = cell.c.cross(cell.a) / determinant
        let rc = cell.a.cross(cell.b) / determinant
        let fx = delta.dot(ra).rounded()
        let fy = delta.dot(rb).rounded()
        let fz = delta.dot(rc).rounded()
        return delta - cell.a * fx - cell.b * fy - cell.c * fz
    }
}
