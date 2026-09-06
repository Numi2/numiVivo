import Foundation

/// Transparent input fixtures, never a claim to reproduce author-supplied data.
public enum VivoBarrierBenchmarks {
    public static func hydrogenExchange631G(ensemble: Bool = false) -> VivoBarrierConvergenceRequest {
        let coordinates: [Double] = [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
        let snapshots = coordinates.enumerated().map { i,u in
            // A declared collinear H + H2 bond-exchange scan. These geometries
            // are not optimized by the campaign and the center is not certified
            // as a saddle simply by being selected for a barrier difference.
            let left = 1.8 + 1.3*u + u*u, right = 1.8 - 1.3*u + u*u
            return VivoMolecularPathSnapshot(identifier: "h3-\(i)", coordinate: u,
                system: .init(nuclei: [
                    .init(atomicNumber: 1, positionBohr: .init(0,0,-left), structureAtomIndex: 0),
                    .init(atomicNumber: 1, positionBohr: .zero, structureAtomIndex: 1),
                    .init(atomicNumber: 1, positionBohr: .init(0,0,right), structureAtomIndex: 2)
                ], alphaElectrons: 2, betaElectrons: 1))
        }
        let inner = [VivoGaussianPrimitive(exponent: 18.731137, coefficient: 0.0334946),
                     .init(exponent: 2.8253937, coefficient: 0.23472695),
                     .init(exponent: 0.6401217, coefficient: 0.81375733)]
        let basis = VivoGaussianBasis(identifier: "H-6-31G-explicit-v1", shells: (0..<3).flatMap { i in
            [VivoGaussianShell(nucleusIndex: i, angularMomentum: 0, primitives: inner),
             .init(nucleusIndex: i, angularMomentum: 0, primitives: [.init(exponent: 0.1612778, coefficient: 1)])]
        }, source: "Explicit hydrogen 6-31G coefficients; independently checked against PySCF 2.8.0 6-31G basis data")
        return .init(identifier: ensemble ? "mapped-H3-6-31G-ensemble-CAS-ladder" : "mapped-H3-6-31G-CAS-ladder", atomIdentifiers: ["H-left","H-center","H-right"],
            coordinateUnit: "dimensionless declared scan parameter", snapshots: snapshots, basis: basis,
            anchorPointIdentifier: "h3-4", barrierPointIdentifier: "h3-4",
            transportGroups: [Array(0..<6)],
            levels: (3...6).map { count in .init(identifier: "CAS-3e-\(count)o",
                method: .casci(partition: .init(active: Array(0..<count)))) },
            acceptance: .init(), maximumPointEvaluations: 45,
            ensembleOrbitals: ensemble ? .init(pointWeights: [Double](repeating: 1.0/9, count: 9)) : nil)
    }
}
