import Foundation

/// An explicitly mapped, separately characterized molecular endpoint component.
/// atomIndices gives the saddle atom index for each atom in point, in order.
/// Components are never inferred from a distance cutoff or permuted to improve a fit.
public struct VivoMappedEndpointComponent: Codable, Sendable, Equatable {
    public let atomIndices: [Int]
    public let point: VivoNuclearQualifiedPoint
    public init(atomIndices: [Int], point: VivoNuclearQualifiedPoint) {
        self.atomIndices = atomIndices; self.point = point
    }
}
public struct VivoMappedReactionEndpoint: Codable, Sendable, Equatable {
    public let identifier: String
    public let components: [VivoMappedEndpointComponent]
    public init(identifier: String, components: [VivoMappedEndpointComponent]) {
        self.identifier = identifier; self.components = components
    }
}

/// Mass-weighted proper-rotation RMSD. The largest eigenvalue of the quaternion
/// alignment matrix maximizes a proper rotation, not an arbitrary reflection.
/// Input atom correspondence is retained, including isotope correspondence.
public enum VivoMappedGeometry {
    public static func properRotationRMSD(_ first: [SIMD3<Double>],
                                          _ second: [SIMD3<Double>], masses: [Double]) throws -> Double {
        guard !first.isEmpty, first.count == second.count, masses.count == first.count,
              first.allSatisfy(vivoQMFinite), second.allSatisfy(vivoQMFinite),
              masses.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoChemistryError.invalid("mapped geometry coordinates or isotopic masses")
        }
        let total = masses.reduce(0, +)
        guard total.isFinite, total > 0 else { throw VivoChemistryError.invalid("mapped geometry total mass") }
        var ca = SIMD3<Double>.zero, cb = SIMD3<Double>.zero
        for i in masses.indices { ca += masses[i] * first[i]; cb += masses[i] * second[i] }
        ca /= total; cb /= total
        var covariance = VivoQMMatrix(3, 3), squaredNorm = 0.0
        for i in masses.indices {
            let a = first[i] - ca, b = second[i] - cb, weight = masses[i] / total
            squaredNorm += weight * (vivoQMDot(a, a) + vivoQMDot(b, b))
            for p in 0..<3 { for q in 0..<3 { covariance[p, q] += weight * a[p] * b[q] } }
        }
        guard squaredNorm.isFinite, covariance.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("mapped alignment overflow")
        }
        let s = covariance
        var k = VivoQMMatrix(4, 4)
        k[0,0] = s[0,0] + s[1,1] + s[2,2]
        k[0,1] = s[1,2] - s[2,1]; k[0,2] = s[2,0] - s[0,2]; k[0,3] = s[0,1] - s[1,0]
        k[1,1] = s[0,0] - s[1,1] - s[2,2]; k[1,2] = s[0,1] + s[1,0]; k[1,3] = s[0,2] + s[2,0]
        k[2,2] = -s[0,0] + s[1,1] - s[2,2]; k[2,3] = s[1,2] + s[2,1]
        k[3,3] = -s[0,0] - s[1,1] + s[2,2]
        for i in 0..<4 { for j in 0..<i { k[i,j] = k[j,i] } }
        let scale = max(k.values.map(abs).max() ?? 0, 1e-300)
        let eigen = try VivoQMDenseAlgebra.symmetricEigen(k.scaled(1 / scale), tolerance: 1e-14)
        let squared = squaredNorm - 2 * eigen.values.last! * scale
        guard squared >= -1e-10 * max(1, squaredNorm) else {
            throw VivoChemistryError.convergence("mapped proper-rotation alignment residual")
        }
        return sqrt(max(0, squared))
    }
}

public struct VivoReactionConnectivityConfiguration: Codable, Sendable, Equatable {
    public var initialDisplacementMassWeighted: Double
    public var stepMassWeighted: Double
    public var maximumArcMassWeighted: Double
    public var endpointMaximumGradient: Double
    public var endpointRMSDBohr: Double
    public var endpointEnergyToleranceHartree: Double
    public var minimumIntercomponentSeparationBohr: Double
    public var comparisonSamples: Int
    public var minimumComparisonArcMassWeighted: Double
    public var profileDistanceToleranceBohr: Double
    public var profileEnergyToleranceHartree: Double
    /// Bounds the four descent runs conservatively before any run starts.
    /// Independent endpoint/saddle reconstruction has its own request budgets.
    public var maximumDescentElectronicEvaluations: Int
    public init(initialDisplacementMassWeighted: Double = 0.03, stepMassWeighted: Double = 0.08,
                maximumArcMassWeighted: Double = 20, endpointMaximumGradient: Double = 1e-4,
                endpointRMSDBohr: Double = 0.02, endpointEnergyToleranceHartree: Double = 1e-5,
                minimumIntercomponentSeparationBohr: Double = 8, comparisonSamples: Int = 64,
                minimumComparisonArcMassWeighted: Double = 1,
                profileDistanceToleranceBohr: Double = 0.005, profileEnergyToleranceHartree: Double = 2e-5,
                maximumDescentElectronicEvaluations: Int = 400000) {
        self.initialDisplacementMassWeighted = initialDisplacementMassWeighted; self.stepMassWeighted = stepMassWeighted
        self.maximumArcMassWeighted = maximumArcMassWeighted; self.endpointMaximumGradient = endpointMaximumGradient
        self.endpointRMSDBohr = endpointRMSDBohr; self.endpointEnergyToleranceHartree = endpointEnergyToleranceHartree
        self.minimumIntercomponentSeparationBohr = minimumIntercomponentSeparationBohr; self.comparisonSamples = comparisonSamples
        self.minimumComparisonArcMassWeighted = minimumComparisonArcMassWeighted
        self.profileDistanceToleranceBohr = profileDistanceToleranceBohr; self.profileEnergyToleranceHartree = profileEnergyToleranceHartree
        self.maximumDescentElectronicEvaluations = maximumDescentElectronicEvaluations
    }
    public func validate() throws {
        guard [initialDisplacementMassWeighted, stepMassWeighted, maximumArcMassWeighted, endpointMaximumGradient,
               endpointRMSDBohr, endpointEnergyToleranceHartree, minimumIntercomponentSeparationBohr,
               minimumComparisonArcMassWeighted, profileDistanceToleranceBohr, profileEnergyToleranceHartree]
                .allSatisfy({ $0.isFinite && $0 > 0 }), stepMassWeighted <= 1,
              maximumArcMassWeighted > initialDisplacementMassWeighted + minimumComparisonArcMassWeighted,
              (3...4096).contains(comparisonSamples), (4...40000000).contains(maximumDescentElectronicEvaluations),
              maximumArcMassWeighted / (stepMassWeighted / 2) < 9998 else {
            throw VivoChemistryError.invalid("reaction connectivity refinement configuration")
        }
    }
}
public struct VivoReactionConnectivityRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mapped-reaction-connectivity/v1"
    public let schema: String
    public let atomIdentifiers: [String]
    public let saddle: VivoNuclearQualifiedPoint
    public let endpoints: [VivoMappedReactionEndpoint]
    public let configuration: VivoReactionConnectivityConfiguration
    public init(atomIdentifiers: [String], saddle: VivoNuclearQualifiedPoint,
                endpoints: [VivoMappedReactionEndpoint], configuration: VivoReactionConnectivityConfiguration = .init()) {
        schema = Self.schema; self.atomIdentifiers = atomIdentifiers; self.saddle = saddle
        self.endpoints = endpoints; self.configuration = configuration
    }
}
public struct VivoEndpointAssignment: Codable, Sendable, Equatable {
    public let endpointIdentifier: String
    public let componentRMSDBohr: [Double]
    public let electronicEnergyDefectHartree: Double
    public let minimumIntercomponentSeparationBohr: Double?
    public let maximumGradient: Double
}
public struct VivoConnectivityTrial: Codable, Sendable, Equatable {
    public let displacementScale: Double
    public let stepScale: Double
    /// Endpoint-completed branches retained for assignment and replay.
    public let descent: VivoNuclearDescentResult
    /// The independently integrated branches used for displacement/step
    /// refinement comparisons, before endpoint completion is appended.
    public let refinementDescent: VivoNuclearDescentResult
    public let reverseAssignment: VivoEndpointAssignment
    public let forwardAssignment: VivoEndpointAssignment
}
public struct VivoConnectivityComparison: Codable, Sendable, Equatable {
    public let firstTrial: Int
    public let secondTrial: Int
    public let endpointIdentifier: String
    public let startingArcMassWeighted: Double
    public let endingArcMassWeighted: Double
    public let samples: Int
    public let maximumPairDistanceDifferenceBohr: Double
    public let maximumEnergyDifferenceHartree: Double
    public let passed: Bool
}
public struct VivoReactionConnectivityResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mapped-reaction-connectivity-result/v1"
    public let schema: String
    public let request: VivoReactionConnectivityRequest
    public let trials: [VivoConnectivityTrial]
    public let comparisons: [VivoConnectivityComparison]
    public let descentElectronicEvaluations: Int
    public let converged: Bool
    public let interpretation: String
}

/// Endpoint identity and sampled numerical path convergence are separate from
/// the local Hessian-index test. Four independent runs vary initial displacement
/// and integration step independently; both directions must map unambiguously to
/// different supplied, independently characterized molecular endpoints.
public enum VivoReactionConnectivity {
    public static let interpretation = "mapped endpoints and sampled mass-weighted descent convergence at the declared displacement, step, geometry and energy tolerances; not a global reaction search, dynamical transmission coefficient or kinetic-rate validation"
    private static func compatible(_ reference: VivoNuclearElectronicModel, _ other: VivoNuclearElectronicModel) -> Bool {
        reference.solver == other.solver && reference.scf == other.scf && reference.solvent == other.solvent
            && reference.correlatedSolventConfiguration == other.correlatedSolventConfiguration
            && reference.eccFrame == nil && other.eccFrame == nil
            && reference.basis.identifier == other.basis.identifier && reference.basis.representation == other.basis.representation
    }
    public static func validateRequest(_ request: VivoReactionConnectivityRequest) throws {
        let cfg = request.configuration, saddle = request.saddle, model = saddle.request.model
        try cfg.validate(); try saddle.request.validate()
        let n = model.system.nuclei.count
        guard request.schema == VivoReactionConnectivityRequest.schema, saddle.request.kind == .firstOrderSaddle,
              request.atomIdentifiers.count == n, Set(request.atomIdentifiers).count == n,
              request.atomIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 1024 }),
              request.endpoints.count == 2, Set(request.endpoints.map(\.identifier)).count == 2,
              model.eccFrame == nil else {
            throw VivoChemistryError.invalid("mapped reaction saddle, atom identities or endpoint count; anchored fragment surfaces need a shared reaction model")
        }
        let work = saddle.request.differences.maximumEnergyEvaluations.multipliedReportingOverflow(by: 4)
        guard !work.overflow, work.partialValue <= cfg.maximumDescentElectronicEvaluations else {
            throw VivoChemistryError.resourceLimit("four descent budgets exceed the declared aggregate cap")
        }
        // Bound retained path coordinates and pair-distance comparison workspace.
        let maximumPoints = Int(ceil(cfg.maximumArcMassWeighted / (cfg.stepMassWeighted / 2))) + 2
        _ = try model.budget.elements([8, maximumPoints, n, 3], simultaneousArrays: 4)
        _ = try model.budget.elements([cfg.comparisonSamples, n, n], simultaneousArrays: 8)
        func shells(_ m: VivoNuclearElectronicModel, _ index: Int) -> [VivoGaussianShell] {
            m.basis.shells.filter { $0.nucleusIndex == index }.map {
                .init(nucleusIndex: 0, angularMomentum: $0.angularMomentum, primitives: $0.primitives)
            }
        }
        for endpoint in request.endpoints {
            guard !endpoint.identifier.isEmpty, endpoint.identifier.utf8.count <= 1024,
                  !endpoint.components.isEmpty, endpoint.components.count <= n else {
                throw VivoChemistryError.invalid("mapped endpoint identity or components")
            }
            let mapping = endpoint.components.flatMap(\.atomIndices)
            guard mapping.count == n, Set(mapping) == Set(0..<n) else {
                throw VivoChemistryError.invalid("endpoint atom maps must cover each saddle atom exactly once")
            }
            var alpha = 0, beta = 0
            for component in endpoint.components {
                let point = component.point, other = point.request.model
                try point.request.validate()
                guard point.request.kind == .minimum, compatible(model, other),
                      component.atomIndices.count == other.system.nuclei.count else {
                    throw VivoChemistryError.invalid("endpoint component is not a compatible mapped minimum")
                }
                alpha += other.system.alphaElectrons; beta += other.system.betaElectrons
                for (local, global) in component.atomIndices.enumerated() {
                    guard other.system.nuclei[local].atomicNumber == model.system.nuclei[global].atomicNumber,
                          point.request.massesDa[local] == saddle.request.massesDa[global],
                          !shells(model, global).isEmpty, shells(other, local) == shells(model, global) else {
                        throw VivoChemistryError.invalid("mapped endpoint changes isotope, element or basis primitives")
                    }
                }
            }
            guard alpha == model.system.alphaElectrons, beta == model.system.betaElectrons else {
                throw VivoChemistryError.invalid("mapped endpoint changes the total fixed-spin electron sector")
            }
        }
    }
    private struct EndpointMetrics {
        let rmsds: [Double]
        let energyDefect: Double
        let minimumSeparation: Double?
    }
    private static func metrics(_ point: VivoNuclearDescentPoint, endpoint: VivoMappedReactionEndpoint) throws -> EndpointMetrics {
        var rmsds: [Double] = [], endpointEnergy = 0.0
        for component in endpoint.components {
            let coordinates = component.atomIndices.map { point.positionsBohr[$0] }
            rmsds.append(try VivoMappedGeometry.properRotationRMSD(coordinates, component.point.finalPositionsBohr,
                                                                   masses: component.point.request.massesDa))
            endpointEnergy += component.point.thermochemistry.electronicEnergyHartree
        }
        var minimumSeparation: Double?
        if endpoint.components.count > 1 {
            var distance = Double.greatestFiniteMagnitude
            for i in endpoint.components.indices { for j in (i + 1)..<endpoint.components.count {
                for a in endpoint.components[i].atomIndices { for b in endpoint.components[j].atomIndices {
                    distance = min(distance, vivoQMNorm(point.positionsBohr[a] - point.positionsBohr[b]))
                } }
            } }
            minimumSeparation = distance
        } else { minimumSeparation = nil }
        return .init(rmsds: rmsds, energyDefect: point.energyHartree - endpointEnergy,
                     minimumSeparation: minimumSeparation)
    }
    private static func assignment(_ point: VivoNuclearDescentPoint, endpoint: VivoMappedReactionEndpoint,
                                   cfg: VivoReactionConnectivityConfiguration) throws -> VivoEndpointAssignment? {
        guard point.energyHartree.isFinite, point.maximumGradient.isFinite, point.maximumGradient >= 0,
              point.maximumGradient <= cfg.endpointMaximumGradient, point.positionsBohr.allSatisfy(vivoQMFinite) else { return nil }
        let metrics = try metrics(point, endpoint: endpoint)
        guard metrics.rmsds.allSatisfy({ $0 <= cfg.endpointRMSDBohr }) else { return nil }
        let defect = metrics.energyDefect
        guard defect.isFinite, abs(defect) <= cfg.endpointEnergyToleranceHartree else { return nil }
        if let separation = metrics.minimumSeparation {
            guard separation >= cfg.minimumIntercomponentSeparationBohr else { return nil }
        }
        return .init(endpointIdentifier: endpoint.identifier, componentRMSDBohr: metrics.rmsds,
                     electronicEnergyDefectHartree: defect, minimumIntercomponentSeparationBohr: metrics.minimumSeparation,
                     maximumGradient: point.maximumGradient)
    }
    private struct EndpointCandidate { let endpoint: VivoMappedReactionEndpoint; let score: Double }
    private static func closestEndpoint(_ point: VivoNuclearDescentPoint,
                                        endpoints: [VivoMappedReactionEndpoint]) throws -> VivoMappedReactionEndpoint {
        let candidates = try endpoints.map { endpoint in
            let metrics = try metrics(point, endpoint: endpoint)
            return EndpointCandidate(endpoint: endpoint, score: metrics.rmsds.reduce(0) { $0 + $1 * $1 })
        }.sorted { $0.score < $1.score }
        guard let first = candidates.first, first.score.isFinite,
              candidates.dropFirst().allSatisfy({ $0.score - first.score > 1e-12 }) else {
            throw VivoChemistryError.convergence("mapped endpoint geometry is not uniquely identified before separation")
        }
        return first.endpoint
    }
    private struct CompletedBranch {
        let points: [VivoNuclearDescentPoint]
        let stationary: Bool
        let evaluations: Int
    }
    private static func completeBranch(_ points: [VivoNuclearDescentPoint], endpoint: VivoMappedReactionEndpoint,
                                      request: VivoReactionConnectivityRequest, maximumEnergyEvaluations: Int) throws -> CompletedBranch {
        guard let last = points.last, maximumEnergyEvaluations > 0 else {
            return .init(points: points, stationary: false, evaluations: 0)
        }
        var differences = request.saddle.request.differences
        differences.maximumEnergyEvaluations = min(maximumEnergyEvaluations, 10_000_000)
        let surface = try VivoNuclearElectronicSurface(model: request.saddle.request.model, differences: differences)
        guard maximumEnergyEvaluations >= 13 else {
            return .init(points: points, stationary: false, evaluations: surface.energyEvaluations)
        }
        let cfg = request.configuration, masses = request.saddle.request.massesDa
        let relaxed = try VivoCartesianGeometry.minimize(positionsBohr: last.positionsBohr,
            configuration: .init(maximumIterations: 256, maximumEvaluations: max(1, min(256, maximumEnergyEvaluations / 13)),
                gradientRMSTolerance: cfg.endpointMaximumGradient, maximumGradientTolerance: cfg.endpointMaximumGradient,
                maximumStepBohr: 0.2), budget: request.saddle.request.model.budget,
            evaluate: { try surface.gradient($0) })
        guard relaxed.converged else {
            return .init(points: points, stationary: false, evaluations: surface.energyEvaluations)
        }
        var completed = points
        let sign = last.arcMassWeighted < 0 ? -1.0 : 1.0
        var currentPositions = relaxed.positionsBohr
        var current = try surface.gradient(currentPositions)
        var currentArc = abs(last.arcMassWeighted)
        func massDistance(_ first: [SIMD3<Double>], _ second: [SIMD3<Double>]) -> Double {
            sqrt(zip(first, second).enumerated().reduce(0.0) { total, pair in
                let delta = pair.1.1 - pair.1.0
                return total + masses[pair.0] * vivoQMDot(delta, delta)
            })
        }
        func maximumGradient(_ evaluation: VivoGeometryEvaluation) -> Double {
            evaluation.gradientHartreePerBohr.flatMap { [$0.x, $0.y, $0.z] }.map(abs).max() ?? .infinity
        }
        func append(_ positions: [SIMD3<Double>], evaluation: VivoGeometryEvaluation) {
            currentArc += massDistance(currentPositions, positions)
            currentPositions = positions
            current = evaluation
            completed.append(.init(positionsBohr: positions, energyHartree: evaluation.energyHartree,
                maximumGradient: maximumGradient(evaluation), arcMassWeighted: sign * currentArc))
        }
        if current.energyHartree > last.energyHartree + 1e-10 {
            return .init(points: points, stationary: false, evaluations: surface.energyEvaluations)
        }
        if massDistance(last.positionsBohr, currentPositions) > 1e-10 {
            append(currentPositions, evaluation: current)
        } else {
            // Replace the raw endpoint with the independently evaluated
            // relaxed point even when minimization moved below the coordinate
            // comparison threshold.  Assignment and separation checks must
            // use the evaluation that produced `current`.
            completed[completed.count - 1] = .init(positionsBohr: currentPositions,
                energyHartree: current.energyHartree, maximumGradient: maximumGradient(current),
                arcMassWeighted: last.arcMassWeighted)
        }
        let target = cfg.minimumIntercomponentSeparationBohr
        let endpointEnergy = endpoint.components.reduce(0.0) { $0 + $1.point.thermochemistry.electronicEnergyHartree }
        if let initial = (try metrics(completed.last!, endpoint: endpoint).minimumSeparation), initial < target {
            for _ in 0..<4 {
                guard let pair = endpoint.components.indices.dropLast().flatMap({ i in
                    endpoint.components.indices.dropFirst(i + 1).flatMap { j in
                        endpoint.components[i].atomIndices.flatMap { a in
                            endpoint.components[j].atomIndices.map { b in (vivoQMNorm(currentPositions[a] - currentPositions[b]), i, j, a, b) }
                        }
                    }
                }).min(by: { $0.0 < $1.0 }), pair.0 < target else { break }
                let raw = currentPositions[pair.3] - currentPositions[pair.4]
                let length = vivoQMNorm(raw)
                guard length > 1e-12 else { return .init(points: completed, stationary: false, evaluations: surface.energyEvaluations) }
                let direction = raw / length, gap = target - pair.0
                var local = min(0.25, gap), accepted = false
                for _ in 0..<20 {
                    var candidate = currentPositions
                    let shift = 0.5 * local * direction
                    for atom in endpoint.components[pair.1].atomIndices { candidate[atom] += shift }
                    for atom in endpoint.components[pair.2].atomIndices { candidate[atom] -= shift }
                    let evaluation = try surface.gradient(candidate)
                    // The unconstrained minimum can lie just inside the declared
                    // fragment-separation boundary because of a tiny residual
                    // intercomponent interaction.  Permit only the declared
                    // endpoint energy-defect budget while moving to that boundary;
                    // assignment still enforces the same separation and gradient
                    // criteria at the completed point.
                    let candidateMetrics = try metrics(.init(positionsBohr: candidate,
                        energyHartree: evaluation.energyHartree, maximumGradient: maximumGradient(evaluation),
                        arcMassWeighted: sign * (currentArc + massDistance(currentPositions, candidate))), endpoint: endpoint)
                    if candidateMetrics.minimumSeparation ?? .greatestFiniteMagnitude > pair.0 + 1e-10,
                       evaluation.energyHartree <= endpointEnergy + cfg.endpointEnergyToleranceHartree {
                        append(candidate, evaluation: evaluation); accepted = true; break
                    }
                    local *= 0.5
                }
                guard accepted else { break }
            }
        }
        let final = completed.last!, finalMetrics = try metrics(final, endpoint: endpoint)
        let stationary = final.maximumGradient <= cfg.endpointMaximumGradient
            && finalMetrics.rmsds.allSatisfy({ $0 <= cfg.endpointRMSDBohr })
            && abs(finalMetrics.energyDefect) <= cfg.endpointEnergyToleranceHartree
            && (finalMetrics.minimumSeparation ?? .greatestFiniteMagnitude) >= target
        return .init(points: completed, stationary: stationary, evaluations: surface.energyEvaluations)
    }
    private static func assign(_ points: [VivoNuclearDescentPoint], request: VivoReactionConnectivityRequest) throws -> VivoEndpointAssignment {
        guard let last = points.last else { throw VivoChemistryError.convergence("empty reaction branch") }
        let matches = try request.endpoints.compactMap { try assignment(last, endpoint: $0, cfg: request.configuration) }
        guard matches.count == 1 else { throw VivoChemistryError.convergence("reaction branch has no unique assigned endpoint at the declared geometry, separation, energy and gradient tolerances") }
        return matches[0]
    }
    private struct ProfileSample { let distances: [Double]; let energy: Double }
    private static func sample(_ points: [VivoNuclearDescentPoint], arc: Double) throws -> ProfileSample {
        guard points.count >= 2, arc >= abs(points[0].arcMassWeighted), arc <= abs(points.last!.arcMassWeighted) else {
            throw VivoChemistryError.invalid("reaction profile interpolation would extrapolate")
        }
        var high = 1
        while high + 1 < points.count && abs(points[high].arcMassWeighted) < arc { high += 1 }
        let a = points[high - 1], b = points[high]
        let lower = abs(a.arcMassWeighted), upper = abs(b.arcMassWeighted)
        guard upper > lower, a.positionsBohr.count == b.positionsBohr.count else {
            throw VivoChemistryError.invalid("reaction profile arc ordering or coordinate shape")
        }
        let fraction = (arc - lower) / (upper - lower), n = a.positionsBohr.count
        var distances: [Double] = []
        for i in 0..<n { for j in (i + 1)..<n {
            let da = vivoQMNorm(a.positionsBohr[i] - a.positionsBohr[j])
            let db = vivoQMNorm(b.positionsBohr[i] - b.positionsBohr[j])
            distances.append(da + fraction * (db - da))
        } }
        return .init(distances: distances, energy: a.energyHartree + fraction * (b.energyHartree - a.energyHartree))
    }
    private static func branch(_ trial: VivoConnectivityTrial, endpoint: String) -> [VivoNuclearDescentPoint] {
        trial.reverseAssignment.endpointIdentifier == endpoint ? trial.refinementDescent.reverse : trial.refinementDescent.forward
    }
    public static func run(_ request: VivoReactionConnectivityRequest) throws -> VivoReactionConnectivityResult {
        try validateRequest(request)
        try VivoNuclearQualification.validate(request.saddle, request: request.saddle.request)
        for endpoint in request.endpoints { for component in endpoint.components {
            try VivoNuclearQualification.validate(component.point, request: component.point.request)
        } }
        let cfg = request.configuration
        let allScales: [(Double, Double)] = [(1, 1), (0.5, 1), (1, 0.5), (0.5, 0.5)]
        let scales = allScales
        var trials: [VivoConnectivityTrial] = [], evaluations = 0
        for (displacementScale, stepScale) in scales {
            let displacement = cfg.initialDisplacementMassWeighted * displacementScale, step = cfg.stepMassWeighted * stepScale
            let steps = Int(ceil((cfg.maximumArcMassWeighted - displacement) / step)) + 1
            let remaining = cfg.maximumDescentElectronicEvaluations - evaluations
            guard remaining > 0 else { throw VivoChemistryError.resourceLimit("four descent aggregate electronic solve budget") }
            let descent = try VivoNuclearDescent.traceValidated(request.saddle, configuration: .init(
                initialDisplacementMassWeighted: displacement, stepMassWeighted: step,
                maximumStepsPerDirection: steps, endpointMaximumGradient: cfg.endpointMaximumGradient),
                maximumEnergyEvaluations: remaining)
            evaluations += descent.energyEvaluations
            guard evaluations < cfg.maximumDescentElectronicEvaluations else {
                throw VivoChemistryError.resourceLimit("four descent aggregate electronic solve budget")
            }
            let reverseEndpoint = try closestEndpoint(descent.reverse.last!, endpoints: request.endpoints)
            let reverseCompleted = try completeBranch(descent.reverse, endpoint: reverseEndpoint, request: request,
                                                     maximumEnergyEvaluations: cfg.maximumDescentElectronicEvaluations - evaluations)
            evaluations += reverseCompleted.evaluations
            guard evaluations < cfg.maximumDescentElectronicEvaluations else {
                throw VivoChemistryError.resourceLimit("four descent aggregate electronic solve budget")
            }
            let forwardEndpoint = try closestEndpoint(descent.forward.last!, endpoints: request.endpoints)
            let forwardCompleted = try completeBranch(descent.forward, endpoint: forwardEndpoint, request: request,
                                                      maximumEnergyEvaluations: cfg.maximumDescentElectronicEvaluations - evaluations)
            evaluations += forwardCompleted.evaluations
            let completedDescent = VivoNuclearDescentResult(reverse: reverseCompleted.points, forward: forwardCompleted.points,
                reverseStationary: reverseCompleted.stationary, forwardStationary: forwardCompleted.stationary,
                energyEvaluations: descent.energyEvaluations + reverseCompleted.evaluations + forwardCompleted.evaluations,
                interpretation: "mass-weighted downhill branches with mapped endpoint relaxation and separation tail; no dynamical rate certification")
            guard completedDescent.reverseStationary, completedDescent.forwardStationary else {
                let reverseMetrics = try metrics(completedDescent.reverse.last!, endpoint: reverseEndpoint)
                let forwardMetrics = try metrics(completedDescent.forward.last!, endpoint: forwardEndpoint)
                throw VivoChemistryError.convergence("mapped completion failed: reverse stationary=\(reverseCompleted.stationary), g=\(completedDescent.reverse.last!.maximumGradient), rmsd=\(reverseMetrics.rmsds), defect=\(reverseMetrics.energyDefect), separation=\(String(describing: reverseMetrics.minimumSeparation)), evals=\(reverseCompleted.evaluations); forward stationary=\(forwardCompleted.stationary), g=\(completedDescent.forward.last!.maximumGradient), rmsd=\(forwardMetrics.rmsds), defect=\(forwardMetrics.energyDefect), separation=\(String(describing: forwardMetrics.minimumSeparation)), evals=\(forwardCompleted.evaluations)")
            }
            let reverse = try assign(completedDescent.reverse, request: request), forward = try assign(completedDescent.forward, request: request)
            guard reverse.endpointIdentifier != forward.endpointIdentifier else {
                throw VivoChemistryError.convergence("both saddle branches reach the same mapped endpoint")
            }
            trials.append(.init(displacementScale: displacementScale, stepScale: stepScale, descent: completedDescent,
                                refinementDescent: descent,
                                reverseAssignment: reverse, forwardAssignment: forward))
        }
        // Every pair is compared, not merely a diagonal refinement where errors
        // in displacement and integration step could cancel each other.
        var comparisons: [VivoConnectivityComparison] = []
        for i in trials.indices { for j in (i + 1)..<trials.count { for endpoint in request.endpoints {
            let a = branch(trials[i], endpoint: endpoint.identifier), b = branch(trials[j], endpoint: endpoint.identifier)
            guard a.count >= 2, b.count >= 2 else { throw VivoChemistryError.convergence("insufficient path points for refinement comparison") }
            // Do not score the singular saddle-launch interpolation.  The
            // declared minimum comparison arc is the lower bound for both
            // displacement and step refinement, while all later samples still
            // share the same physical mass-weighted arc.
            let lower = max(cfg.minimumComparisonArcMassWeighted,
                            max(abs(a[0].arcMassWeighted), abs(b[0].arcMassWeighted)))
            let upper = min(abs(a.last!.arcMassWeighted), abs(b.last!.arcMassWeighted))
            guard upper - lower >= cfg.minimumComparisonArcMassWeighted else {
                throw VivoChemistryError.convergence("reaction branches have insufficient common physical arc for convergence testing")
            }
            var geometryError = 0.0, energyError = 0.0
            for k in 0..<cfg.comparisonSamples {
                let arc = k == cfg.comparisonSamples - 1 ? upper : lower + (upper - lower) * Double(k) / Double(cfg.comparisonSamples - 1)
                let x = try sample(a, arc: arc), y = try sample(b, arc: arc)
                guard x.distances.count == y.distances.count else { throw VivoChemistryError.invalid("reaction comparison atom count") }
                for (first, second) in zip(x.distances, y.distances) { geometryError = max(geometryError, abs(first - second)) }
                energyError = max(energyError, abs(x.energy - y.energy))
            }
            guard geometryError.isFinite, energyError.isFinite else { throw VivoChemistryError.convergence("nonfinite reaction refinement error") }
            comparisons.append(.init(firstTrial: i, secondTrial: j, endpointIdentifier: endpoint.identifier,
                startingArcMassWeighted: lower, endingArcMassWeighted: upper, samples: cfg.comparisonSamples,
                maximumPairDistanceDifferenceBohr: geometryError, maximumEnergyDifferenceHartree: energyError,
                passed: geometryError <= cfg.profileDistanceToleranceBohr && energyError <= cfg.profileEnergyToleranceHartree))
        } } }
        return .init(schema: VivoReactionConnectivityResult.schema, request: request, trials: trials, comparisons: comparisons,
                     descentElectronicEvaluations: evaluations, converged: comparisons.allSatisfy(\.passed), interpretation: interpretation)
    }
    public static func validate(_ result: VivoReactionConnectivityResult, request: VivoReactionConnectivityRequest) throws {
        guard result.schema == VivoReactionConnectivityResult.schema, result.request == request, result.converged,
              result.interpretation == interpretation, result.trials.count == 4, result.comparisons.count == 12 else {
            throw VivoChemistryError.invalid("reaction connectivity result identity or qualification contract")
        }
        // Deterministic reconstruction checks endpoint assignment and all four
        // numerical paths. User-edited assignment/acceptance flags are not authority.
        let rebuilt = try run(request)
        guard rebuilt.converged, rebuilt == result else {
            throw VivoChemistryError.invalid("reaction connectivity reconstruction differs or is not converged")
        }
    }
}
