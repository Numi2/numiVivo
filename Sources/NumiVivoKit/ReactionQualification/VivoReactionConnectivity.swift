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
    public let descent: VivoNuclearDescentResult
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
    private static func assignment(_ point: VivoNuclearDescentPoint, endpoint: VivoMappedReactionEndpoint,
                                   cfg: VivoReactionConnectivityConfiguration) throws -> VivoEndpointAssignment? {
        guard point.energyHartree.isFinite, point.maximumGradient.isFinite, point.maximumGradient >= 0,
              point.maximumGradient <= cfg.endpointMaximumGradient, point.positionsBohr.allSatisfy(vivoQMFinite) else { return nil }
        var rmsds: [Double] = [], energy = 0.0
        for component in endpoint.components {
            let coordinates = component.atomIndices.map { point.positionsBohr[$0] }
            let rmsd = try VivoMappedGeometry.properRotationRMSD(coordinates, component.point.finalPositionsBohr,
                                                               masses: component.point.request.massesDa)
            if rmsd > cfg.endpointRMSDBohr { return nil }
            rmsds.append(rmsd); energy += component.point.thermochemistry.electronicEnergyHartree
        }
        let defect = point.energyHartree - energy
        guard defect.isFinite, abs(defect) <= cfg.endpointEnergyToleranceHartree else { return nil }
        var minimum: Double?
        if endpoint.components.count > 1 {
            var distance = Double.greatestFiniteMagnitude
            for i in endpoint.components.indices { for j in (i + 1)..<endpoint.components.count {
                for a in endpoint.components[i].atomIndices { for b in endpoint.components[j].atomIndices {
                    distance = min(distance, vivoQMNorm(point.positionsBohr[a] - point.positionsBohr[b]))
                } }
            } }
            guard distance >= cfg.minimumIntercomponentSeparationBohr else { return nil }
            minimum = distance
        }
        return .init(endpointIdentifier: endpoint.identifier, componentRMSDBohr: rmsds,
                     electronicEnergyDefectHartree: defect, minimumIntercomponentSeparationBohr: minimum,
                     maximumGradient: point.maximumGradient)
    }
    private static func assign(_ points: [VivoNuclearDescentPoint], request: VivoReactionConnectivityRequest) throws -> VivoEndpointAssignment {
        guard let last = points.last else { throw VivoChemistryError.convergence("empty reaction branch") }
        let matches = try request.endpoints.compactMap { try assignment(last, endpoint: $0, cfg: request.configuration) }
        guard matches.count == 1 else {
            throw VivoChemistryError.convergence("reaction branch has no unique assigned endpoint at the declared geometry, separation, energy and gradient tolerances")
        }
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
        trial.reverseAssignment.endpointIdentifier == endpoint ? trial.descent.reverse : trial.descent.forward
    }
    public static func run(_ request: VivoReactionConnectivityRequest) throws -> VivoReactionConnectivityResult {
        try validateRequest(request)
        try VivoNuclearQualification.validate(request.saddle, request: request.saddle.request)
        for endpoint in request.endpoints { for component in endpoint.components {
            try VivoNuclearQualification.validate(component.point, request: component.point.request)
        } }
        let cfg = request.configuration
        let scales: [(Double, Double)] = [(1, 1), (0.5, 1), (1, 0.5), (0.5, 0.5)]
        var trials: [VivoConnectivityTrial] = [], evaluations = 0
        for (displacementScale, stepScale) in scales {
            let displacement = cfg.initialDisplacementMassWeighted * displacementScale, step = cfg.stepMassWeighted * stepScale
            let steps = Int(ceil((cfg.maximumArcMassWeighted - displacement) / step)) + 1
            let descent = try VivoNuclearDescent.trace(request.saddle, configuration: .init(
                initialDisplacementMassWeighted: displacement, stepMassWeighted: step,
                maximumStepsPerDirection: steps, endpointMaximumGradient: cfg.endpointMaximumGradient))
            evaluations += descent.energyEvaluations
            guard evaluations <= cfg.maximumDescentElectronicEvaluations,
                  descent.reverseStationary, descent.forwardStationary else {
                throw VivoChemistryError.convergence("reaction endpoint descent exhausted its budget or did not reach both gradient thresholds")
            }
            let reverse = try assign(descent.reverse, request: request), forward = try assign(descent.forward, request: request)
            guard reverse.endpointIdentifier != forward.endpointIdentifier else {
                throw VivoChemistryError.convergence("both saddle branches reach the same mapped endpoint")
            }
            trials.append(.init(displacementScale: displacementScale, stepScale: stepScale, descent: descent,
                                reverseAssignment: reverse, forwardAssignment: forward))
        }
        // Every pair is compared, not merely a diagonal refinement where errors
        // in displacement and integration step could cancel each other.
        var comparisons: [VivoConnectivityComparison] = []
        for i in trials.indices { for j in (i + 1)..<trials.count { for endpoint in request.endpoints {
            let a = branch(trials[i], endpoint: endpoint.identifier), b = branch(trials[j], endpoint: endpoint.identifier)
            guard a.count >= 2, b.count >= 2 else { throw VivoChemistryError.convergence("insufficient path points for refinement comparison") }
            let lower = max(abs(a[0].arcMassWeighted), abs(b[0].arcMassWeighted))
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
