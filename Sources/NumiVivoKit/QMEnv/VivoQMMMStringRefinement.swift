import Foundation

public struct VivoQMMMStringNodeForce: Codable, Sendable, Equatable {
    public var nodeIdentifier: String
    /// Mean negative free-energy gradient for each component q_d, kJ/mol/nm.
    public var meanForceKJPerMolNM: [Double]
    public var standardErrorKJPerMolNM: [Double]
    public var effectiveSamples: Double
    public init(nodeIdentifier: String,
                meanForceKJPerMolNM: [Double],
                standardErrorKJPerMolNM: [Double],
                effectiveSamples: Double) {
        self.nodeIdentifier = nodeIdentifier
        self.meanForceKJPerMolNM = meanForceKJPerMolNM
        self.standardErrorKJPerMolNM = standardErrorKJPerMolNM
        self.effectiveSamples = effectiveSamples
    }
}

public struct VivoQMMMStringRefinementConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-string-refinement/v2"
    public var schema: String
    /// Step in scaled-coordinate space per kJ/mol of generalized scaled force.
    public var mobilityPerKJPerMol: Double
    /// Optional dimensionless discrete-Laplacian smoothing strength in [0,0.5].
    public var smoothing: Double
    /// Maximum norm of one provisional scaled-node displacement.
    public var maximumScaledNodeDisplacement: Double
    public var minimumEffectiveSamples: Double
    /// Convergence threshold on scaled perpendicular force norm, kJ/mol.
    public var perpendicularForceToleranceKJPerMol: Double
    /// Conservative standard-deviation guard; not an asserted confidence interval.
    public var standardErrorMultiplier: Double
    /// Includes smoothing and arc-length reparameterization, not just force motion.
    public var convergenceScaledDisplacementTolerance: Double
    public init(mobilityPerKJPerMol: Double = 0.002,
                smoothing: Double = 0.05,
                maximumScaledNodeDisplacement: Double = 0.25,
                minimumEffectiveSamples: Double = 50,
                perpendicularForceToleranceKJPerMol: Double = 0.5,
                standardErrorMultiplier: Double = 2,
                convergenceScaledDisplacementTolerance: Double = 0.001) {
        schema = Self.schema
        self.mobilityPerKJPerMol = mobilityPerKJPerMol
        self.smoothing = smoothing
        self.maximumScaledNodeDisplacement = maximumScaledNodeDisplacement
        self.minimumEffectiveSamples = minimumEffectiveSamples
        self.perpendicularForceToleranceKJPerMol = perpendicularForceToleranceKJPerMol
        self.standardErrorMultiplier = standardErrorMultiplier
        self.convergenceScaledDisplacementTolerance = convergenceScaledDisplacementTolerance
    }
    public func validate() throws {
        guard schema == Self.schema,
              mobilityPerKJPerMol.isFinite, mobilityPerKJPerMol > 0, mobilityPerKJPerMol <= 1,
              smoothing.isFinite, smoothing >= 0, smoothing <= 0.5,
              maximumScaledNodeDisplacement.isFinite, maximumScaledNodeDisplacement > 0, maximumScaledNodeDisplacement <= 10,
              minimumEffectiveSamples.isFinite, minimumEffectiveSamples >= 1,
              perpendicularForceToleranceKJPerMol.isFinite, perpendicularForceToleranceKJPerMol > 0,
              standardErrorMultiplier.isFinite, standardErrorMultiplier >= 1, standardErrorMultiplier <= 10,
              convergenceScaledDisplacementTolerance.isFinite, convergenceScaledDisplacementTolerance > 0 else {
            throw VivoChemistryError.invalid("QM/MM string-refinement configuration")
        }
    }
}

public struct VivoQMMMStringNodeDiagnostic: Codable, Sendable, Equatable {
    public let nodeIdentifier: String
    public let scaledPerpendicularForceNormKJPerMol: Double
    public let perpendicularForceStandardDeviationUpperBoundKJPerMol: Double
    public let uncertaintyGuardedPerpendicularForceNormKJPerMol: Double
    public let scaledDisplacementNorm: Double
    public let effectiveSamples: Double
    public let maximumComponentStandardErrorKJPerMolNM: Double
}

public struct VivoQMMMStringRefinementResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-string-refinement-result/v2"
    public let schema: String
    public let sourcePathFingerprint: VivoFingerprint
    public let refinedPath: VivoQMMMPathCoordinate
    public let diagnostics: [VivoQMMMStringNodeDiagnostic]
    public let maximumPerpendicularForceNormKJPerMol: Double
    public let maximumUncertaintyGuardedPerpendicularForceNormKJPerMol: Double
    public let maximumScaledDisplacement: Double
    public let converged: Bool
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// One deterministic finite-temperature string update from caller-supplied mean
/// generalized forces. Sampling/force estimation is separate: this function
/// refuses undersampled nodes and does not invent forces from path geometry.
public enum VivoQMMMStringRefinement {
    public static let interpretation = "Finite-temperature string refinement in dimension-scaled collective-variable space. Endpoints are fixed. Interior mean forces are projected perpendicular to the local path tangent, optionally Laplacian-smoothed, displacement-limited, and reparameterized to equal scaled arc length. Convergence requires a covariance-agnostic standard-error guard on perpendicular force and small final displacement, including reparameterization. This is refinement in the declared Euclidean scaled metric; it does not infer a molecular CV mobility tensor or an unrestrained minimum-free-energy path."

    public static func refine(path: VivoQMMMPathCoordinate,
                              nodeForces: [VivoQMMMStringNodeForce],
                              configuration: VivoQMMMStringRefinementConfiguration = .init()) throws -> VivoQMMMStringRefinementResult {
        try path.validate(); try configuration.validate()
        let dimensions = path.components.count, count = path.nodes.count
        guard count >= 3, nodeForces.count == count - 2,
              Set(nodeForces.map(\.nodeIdentifier)).count == nodeForces.count else {
            throw VivoChemistryError.invalid("string refinement requires exactly one force record per interior node")
        }
        let byID = Dictionary(uniqueKeysWithValues: nodeForces.map { ($0.nodeIdentifier, $0) })
        for index in 1..<(count - 1) {
            guard let force = byID[path.nodes[index].identifier],
                  force.meanForceKJPerMolNM.count == dimensions,
                  force.standardErrorKJPerMolNM.count == dimensions,
                  force.meanForceKJPerMolNM.allSatisfy(\.isFinite),
                  force.standardErrorKJPerMolNM.allSatisfy({ $0.isFinite && $0 >= 0 }),
                  force.effectiveSamples.isFinite,
                  force.effectiveSamples >= configuration.minimumEffectiveSamples else {
                throw VivoChemistryError.convergence("string refinement node is missing, malformed or undersampled")
            }
        }

        let y = path.nodes.map { node in
            zip(node.valuesNM, path.scalesNM).map { $0.0 / $0.1 }
        }
        let originalY = y
        var provisional = y
        var diagnostics: [VivoQMMMStringNodeDiagnostic] = []
        var maximumForce = 0.0, maximumGuardedForce = 0.0

        func norm(_ vector: [Double]) -> Double { sqrt(vector.reduce(0) { $0 + $1*$1 }) }
        func dot(_ a: [Double], _ b: [Double]) -> Double { zip(a,b).reduce(0) { $0 + $1.0*$1.1 } }

        for i in 1..<(count - 1) {
            guard let observation = byID[path.nodes[i].identifier] else {
                throw VivoChemistryError.invalid("missing interior string force")
            }
            let tangent0 = (0..<dimensions).map { y[i+1][$0] - y[i-1][$0] }
            let tangentNorm = norm(tangent0)
            guard tangentNorm.isFinite, tangentNorm > 1e-12 else {
                throw VivoChemistryError.invalid("string path has a degenerate interior tangent")
            }
            let tangent = tangent0.map { $0 / tangentNorm }
            // q_d = scale_d*y_d, therefore -dG/dy_d = F_qd*scale_d.
            let scaledForce = (0..<dimensions).map { observation.meanForceKJPerMolNM[$0] * path.scalesNM[$0] }
            let parallel = dot(scaledForce, tangent)
            let perpendicular = (0..<dimensions).map { scaledForce[$0] - parallel*tangent[$0] }
            let forceNorm = norm(perpendicular)
            maximumForce = max(maximumForce, forceNorm)
            // ||P e_j|| = sqrt(1-t_j^2). Minkowski gives this upper bound
            // on RMS projected error for ANY covariance compatible with the
            // supplied component standard errors. Independence is not assumed.
            var errorBound = 0.0
            for d in 0..<dimensions {
                let scaledError = observation.standardErrorKJPerMolNM[d] * path.scalesNM[d]
                let columnNorm = sqrt(max(0.0, 1.0 - tangent[d]*tangent[d]))
                errorBound += scaledError * columnNorm
            }
            let guardedForce = forceNorm + configuration.standardErrorMultiplier * errorBound
            guard forceNorm.isFinite, errorBound.isFinite, guardedForce.isFinite else {
                throw VivoChemistryError.convergence("string force or uncertainty overflow")
            }
            maximumGuardedForce = max(maximumGuardedForce, guardedForce)
            var displacement = perpendicular.map { configuration.mobilityPerKJPerMol * $0 }
            if configuration.smoothing > 0 {
                for d in 0..<dimensions {
                    displacement[d] += configuration.smoothing * (y[i-1][d] - 2*y[i][d] + y[i+1][d])
                }
            }
            let displacementNorm = norm(displacement)
            if displacementNorm > configuration.maximumScaledNodeDisplacement {
                let factor = configuration.maximumScaledNodeDisplacement / displacementNorm
                displacement = displacement.map { $0 * factor }
            }
            for d in 0..<dimensions { provisional[i][d] += displacement[d] }
            diagnostics.append(.init(nodeIdentifier: path.nodes[i].identifier,
                scaledPerpendicularForceNormKJPerMol: forceNorm,
                perpendicularForceStandardDeviationUpperBoundKJPerMol: errorBound,
                uncertaintyGuardedPerpendicularForceNormKJPerMol: guardedForce,
                scaledDisplacementNorm: norm(displacement),
                effectiveSamples: observation.effectiveSamples,
                maximumComponentStandardErrorKJPerMolNM: observation.standardErrorKJPerMolNM.max() ?? 0))
        }

        // Equal-arc-length reparameterization of the provisional polyline.
        var cumulative = [Double](repeating: 0, count: count)
        for i in 1..<count {
            let segment = (0..<dimensions).map { provisional[i][$0] - provisional[i-1][$0] }
            let length = norm(segment)
            guard length.isFinite, length > 1e-12 else {
                throw VivoChemistryError.invalid("string refinement collapsed adjacent nodes")
            }
            cumulative[i] = cumulative[i-1] + length
        }
        guard let total = cumulative.last, total.isFinite, total > 0 else {
            throw VivoChemistryError.invalid("string refinement path length")
        }
        var reparameterized = provisional
        reparameterized[0] = provisional[0]
        reparameterized[count-1] = provisional[count-1]
        for targetIndex in 1..<(count-1) {
            let target = total * Double(targetIndex) / Double(count-1)
            var upper = 1
            while upper < count && cumulative[upper] < target { upper += 1 }
            guard upper < count else { throw VivoChemistryError.invalid("string arc-length interpolation") }
            let lower = upper - 1
            let span = cumulative[upper] - cumulative[lower]
            guard span > 0 else { throw VivoChemistryError.invalid("string zero arc-length segment") }
            let t = (target - cumulative[lower]) / span
            reparameterized[targetIndex] = (0..<dimensions).map {
                provisional[lower][$0] * (1-t) + provisional[upper][$0] * t
            }
        }

        var refined = path
        for i in 0..<count {
            refined.nodes[i].valuesNM = (0..<dimensions).map { reparameterized[i][$0] * path.scalesNM[$0] }
        }
        try refined.validate()
        let maximumDisplacement = (0..<count).map { i in
            norm((0..<dimensions).map { reparameterized[i][$0] - originalY[i][$0] })
        }.max() ?? 0
        guard maximumDisplacement.isFinite else { throw VivoChemistryError.convergence("string displacement overflow") }

        let sourceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(path))
        struct Evidence: Codable {
            let schema: String
            let path: VivoQMMMPathCoordinate
            let observations: [VivoQMMMStringNodeForce]
            let configuration: VivoQMMMStringRefinementConfiguration
            let refined: VivoQMMMPathCoordinate
            let diagnostics: [VivoQMMMStringNodeDiagnostic]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/qmmm-string-refinement-evidence/v2", path: path,
            observations: nodeForces, configuration: configuration,
            refined: refined, diagnostics: diagnostics)))
        return .init(schema: VivoQMMMStringRefinementResult.schema,
            sourcePathFingerprint: sourceID,
            refinedPath: refined,
            diagnostics: diagnostics,
            maximumPerpendicularForceNormKJPerMol: maximumForce,
            maximumUncertaintyGuardedPerpendicularForceNormKJPerMol: maximumGuardedForce,
            maximumScaledDisplacement: maximumDisplacement,
            converged: maximumGuardedForce <= configuration.perpendicularForceToleranceKJPerMol
                && maximumDisplacement <= configuration.convergenceScaledDisplacementTolerance,
            interpretation: interpretation,
            evidenceFingerprint: evidenceID)
    }
}
