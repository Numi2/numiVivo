import Foundation

/// Scientific/backend compatibility checks shared by runtime creation, CLI
/// capability reporting and staged-workflow validation. No GPU is required.
public enum VivoMDExecutionPreflight {
    public static func validateCell(_ cell: VivoPeriodicCell, configuration: VivoMDConfiguration) throws {
        guard cell.isValid else { throw VivoMDRuntimeError.unsupported(["periodic cell is invalid"]) }
        let vectors=[cell.a,cell.b,cell.c]
        for vector in vectors {
            for x in [vector.x,vector.y,vector.z] {
                guard Float(x).isFinite else { throw VivoMDRuntimeError.unsupported(["periodic cell is outside FP32"])}
            }
        }
        // The current force, constraint and virtual-site kernels independently
        // round fractional coordinates. That is an exact nearest-image rule for
        // orthogonal cells, NOT arbitrary skew lattices. Do not run the latter
        // merely because their archival representation accepts triclinic cells.
        for (i,j) in [(0,1),(0,2),(1,2)] {
            let cosine=vectors[i].dot(vectors[j])/(vectors[i].norm*vectors[j].norm)
            guard cosine.isFinite,abs(cosine)<=1e-8 else {
                throw VivoMDRuntimeError.unsupported(["skew triclinic MD needs a nearest-lattice-vector backend; this numerical profile supports orthogonal cells"])
            }
        }
        try VivoMDMetalABI.validateMinimumImageCutoff(configuration.cutoffNM,cell:cell)
        if configuration.resolvedNeighborListEnabled {
            try VivoMDMetalABI.validateNeighborRadius(configuration.cutoffNM+configuration.neighborSkinNM,cell:cell)
        }
    }

    public static func blockers(system:VivoClassicalSystem,initial:VivoClassicalInitialState,
                                configuration:VivoMDConfiguration) -> [String] {
        var result:[String]=[]
        if configuration.thermostat == .velocityRescale {result.append("velocity-rescale thermostat is not implemented")}
        if configuration.electrostatics == .pme && !configuration.resolvedNeighborListEnabled {
            result.append("PME real-space evaluation requires neighbor-list execution")
        }
        if configuration.electrostatics == .pme && initial.periodicCell == nil {result.append("PME requires a periodic cell")}
        if configuration.ensemble == .npt && initial.periodicCell == nil {result.append("NPT requires a periodic cell")}
        if let cell=initial.periodicCell {
            do {try validateCell(cell,configuration:configuration)} catch {result.append(String(describing:error))}
        }
        if Set(system.particles.map(\.typeIdentifier)).count>65_535 {result.append("dense LJ type-pair address exceeds the UInt32 shader index")}
        let massive=system.particles.filter{$0.massDa>0}.count
        if UInt64(system.constraints.count)>=3*UInt64(massive) {result.append("constraint count leaves no positive thermal degrees of freedom")}
        for constraint in system.constraints where system.particles[Int(constraint.a)].massDa==0 || system.particles[Int(constraint.b)].massDa==0 {
            result.append("holonomic constraints cannot independently move massless particles")
            break
        }
        return result
    }
}
