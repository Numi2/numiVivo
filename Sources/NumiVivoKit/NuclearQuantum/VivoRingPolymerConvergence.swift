import Foundation

public enum VivoRingPolymerConvergenceObservable: Codable, Sendable, Equatable {
    case primitiveTotalEnergyKJPerMol
    case meanPotentialEnergyKJPerMol
    case radiusOfGyrationSquaredNM2
    /// Independent-atom slots in the bound nuclear-potential definition.
    case centroidDistanceNM(atomA: Int, atomB: Int)
}

public struct VivoRingPolymerConvergenceChain: Codable, Sendable, Equatable {
    public let identifier: String
    public let run: VivoRingPolymerRun
    public let recordEvery: Int
    public let discardedObservations: Int
    public init(identifier: String, run: VivoRingPolymerRun, recordEvery: Int = 1,
                discardedObservations: Int) {
        self.identifier=identifier;self.run=run;self.recordEvery=recordEvery
        self.discardedObservations=discardedObservations
    }
}

public struct VivoRingPolymerConvergenceConfiguration: Codable, Sendable, Equatable {
    public let observable: VivoRingPolymerConvergenceObservable
    public let absoluteTolerance: Double
    public let relativeTolerance: Double
    public let confidenceMultiplier: Double
    public let minimumEffectiveSamples: Double
    public let maximumSplitRHat: Double
    public let requiredConsecutiveRefinements: Int
    public let maximumPrimitiveWork: Int
    public init(observable: VivoRingPolymerConvergenceObservable, absoluteTolerance: Double,
                relativeTolerance: Double, confidenceMultiplier: Double = 2,
                minimumEffectiveSamples: Double = 100, maximumSplitRHat: Double = 1.05,
                requiredConsecutiveRefinements: Int = 2, maximumPrimitiveWork: Int = 100_000_000) {
        self.observable=observable;self.absoluteTolerance=absoluteTolerance
        self.relativeTolerance=relativeTolerance;self.confidenceMultiplier=confidenceMultiplier
        self.minimumEffectiveSamples=minimumEffectiveSamples;self.maximumSplitRHat=maximumSplitRHat
        self.requiredConsecutiveRefinements=requiredConsecutiveRefinements
        self.maximumPrimitiveWork=maximumPrimitiveWork
    }
    public func validate() throws {
        guard absoluteTolerance.isFinite, absoluteTolerance >= 0,
              relativeTolerance.isFinite, relativeTolerance >= 0,
              absoluteTolerance+relativeTolerance > 0,
              confidenceMultiplier.isFinite, confidenceMultiplier >= 1, confidenceMultiplier <= 10,
              minimumEffectiveSamples.isFinite, minimumEffectiveSamples >= 4,
              maximumSplitRHat.isFinite, maximumSplitRHat >= 1, maximumSplitRHat <= 2,
              (1...8).contains(requiredConsecutiveRefinements),maximumPrimitiveWork > 0 else {
            throw VivoChemistryError.invalid("ring-polymer convergence criteria")
        }
    }
}

public struct VivoRingPolymerConvergenceRequest: Codable, Sendable, Equatable {
    public let campaignIdentifier: String
    public let selectionProtocol: String
    public let chains: [VivoRingPolymerConvergenceChain]
    public let configuration: VivoRingPolymerConvergenceConfiguration
    public init(campaignIdentifier: String, selectionProtocol: String,
                chains: [VivoRingPolymerConvergenceChain],
                configuration: VivoRingPolymerConvergenceConfiguration) {
        self.campaignIdentifier=campaignIdentifier;self.selectionProtocol=selectionProtocol
        self.chains=chains;self.configuration=configuration
    }
}

public struct VivoRingPolymerBeadDiagnostics: Codable, Sendable, Equatable {
    public let beadCount: Int
    public let chainIdentifiers: [String]
    public let startCheckpointFingerprints: [VivoFingerprint]
    public let endCheckpointFingerprints: [VivoFingerprint]
    public let recordEveryPerChain: [Int]
    public let discardedObservationsPerChain: [Int]
    public let acceptanceFraction: Double
    public let samplingPassed: Bool
    public let statistics: VivoCorrelatedSamplingDiagnostics
}

public struct VivoRingPolymerBeadComparison: Codable, Sendable, Equatable {
    public let lowerBeadCount: Int
    public let upperBeadCount: Int
    public let absoluteMeanDifference: Double
    public let combinedMonteCarloStandardError: Double
    public let upperDifferenceBound: Double
    public let allowedDifference: Double
    public let passed: Bool
}

public struct VivoRingPolymerConvergenceResult: Codable, Sendable, Equatable {
    public let campaignIdentifier: String
    public let potentialDefinitionFingerprint: VivoFingerprint
    public let temperatureK: Double
    public let configuration: VivoRingPolymerConvergenceConfiguration
    public let beadDiagnostics: [VivoRingPolymerBeadDiagnostics]
    public let comparisons: [VivoRingPolymerBeadComparison]
    public let passed: Bool
    public let requestFingerprint: VivoFingerprint
    public let evidenceFingerprint: VivoFingerprint
    public let selectionProtocol: String
    public let interpretation: String
}

/// Post-processes retained Markov states; it does not run more dynamics or
/// silently relabel a finite-P result as a converged quantum observable.
public enum VivoRingPolymerConvergence {
    public static let interpretation = "Independent-chain diagnostics and consecutive finite-bead comparisons for one declared scalar observable. The difference bound includes supplied-chain Monte Carlo error. Distinct seeds and chain identifiers are checked but do not prove overdispersed initialization, equilibration, rank-normalized tail convergence, Hamiltonian accuracy, isotope overlap or a real-time quantum rate."

    public static func assess(campaignIdentifier: String, selectionProtocol: String,
                              chains: [VivoRingPolymerConvergenceChain],
                              configuration: VivoRingPolymerConvergenceConfiguration) throws -> VivoRingPolymerConvergenceResult {
        try configuration.validate()
        func validText(_ value: String) -> Bool {
            !value.trimmingCharacters(in:.whitespacesAndNewlines).isEmpty && value.utf8.count <= 4096
        }
        guard validText(campaignIdentifier), validText(selectionProtocol),
              (4...1024).contains(chains.count), Set(chains.map(\.identifier)).count == chains.count else {
            throw VivoChemistryError.invalid("ring-polymer convergence campaign identity or chains")
        }
        for chain in chains {
            guard validText(chain.identifier), chain.recordEvery > 0,
                  chain.discardedObservations >= 0,
                  chain.discardedObservations < chain.run.observations.count else {
                throw VivoChemistryError.invalid("ring-polymer discard or chain identity")
            }
            try chain.run.validate(recordEvery: chain.recordEvery)
        }
        guard let first = chains.first else { throw VivoChemistryError.invalid("ring-polymer convergence chains") }
        let definition = first.run.start.definition, temperature = first.run.start.configuration.temperatureK
        let definitionID = try definition.fingerprint()
        guard chains.allSatisfy({ $0.run.start.definition == definition
            && $0.run.start.configuration.temperatureK == temperature }),
              Set(chains.map { $0.run.start.fingerprint }).count == chains.count,
              Set(chains.map { $0.run.start.randomState.state }).count == chains.count else {
            throw VivoChemistryError.invalid("ring-polymer convergence Hamiltonian, temperature or independent starts")
        }
        let beadCounts = Array(Set(chains.map { $0.run.start.configuration.beadCount })).sorted()
        guard beadCounts.count >= configuration.requiredConsecutiveRefinements+1 else {
            throw VivoChemistryError.invalid("ring-polymer convergence requires enough bead refinements")
        }
        let retainedCount=chains.reduce(0) { $0+$1.run.observations.count-$1.discardedObservations }
        let work=Double(retainedCount)*20+Double(beadCounts.count*VivoCorrelatedSamplingAnalysis.maximumAutocorrelationProducts)
        guard work <= Double(configuration.maximumPrimitiveWork) else {
            throw VivoChemistryError.resourceLimit("ring-polymer convergence analysis budget")
        }
        func value(_ observation: VivoRingPolymerObservation) throws -> Double {
            switch configuration.observable {
            case .primitiveTotalEnergyKJPerMol: return observation.primitiveTotalEnergyKJPerMol
            case .meanPotentialEnergyKJPerMol: return observation.meanPotentialEnergyKJPerMol
            case .radiusOfGyrationSquaredNM2:
                let masses=definition.massesDa,totalMass=masses.reduce(0,+),p=Double(observation.beadPositionsNM.count)
                let sum=observation.beadPositionsNM.reduce(0.0) { partial,bead in
                    partial+zip(bead,zip(observation.centroidPositionsNM,masses)).reduce(0.0) { total,item in
                        total+item.1.1*(item.0-item.1.0).squaredNorm
                    }
                }
                return sum/(p*totalMass)
            case .centroidDistanceNM(let a,let b):
                guard observation.centroidPositionsNM.indices.contains(a),
                      observation.centroidPositionsNM.indices.contains(b), a != b else {
                    throw VivoChemistryError.invalid("ring-polymer centroid-distance atom slots")
                }
                return (observation.centroidPositionsNM[a]-observation.centroidPositionsNM[b]).norm
            }
        }
        var levels: [VivoRingPolymerBeadDiagnostics] = []
        for beadCount in beadCounts {
            let selected=chains.filter { $0.run.start.configuration.beadCount == beadCount }.sorted { $0.identifier < $1.identifier }
            guard selected.count >= 2 else { throw VivoChemistryError.invalid("ring-polymer bead level requires independent chains") }
            let retained=try selected.map { chain in
                try Array(chain.run.observations.dropFirst(chain.discardedObservations)).map(value)
            }
            guard let retainedCount=retained.first?.count, retained.allSatisfy({ $0.count == retainedCount }) else {
                throw VivoChemistryError.invalid("ring-polymer bead level requires equal retained chain lengths")
            }
            let statistics=try VivoCorrelatedSamplingAnalysis.calculate(chains:retained)
            let accepted=selected.reduce(0) { $0+$1.run.observations.dropFirst($1.discardedObservations).filter(\.accepted).count }
            let total=selected.reduce(0) { $0+$1.run.observations.count-$1.discardedObservations }
            let samplingPassed=statistics.autocorrelationEffectiveSampleSize >= configuration.minimumEffectiveSamples
                && statistics.splitRHat <= configuration.maximumSplitRHat
                && !statistics.autocorrelationSequenceTruncated
            levels.append(.init(beadCount:beadCount,chainIdentifiers:selected.map(\.identifier),
                startCheckpointFingerprints:selected.map { $0.run.start.fingerprint },
                endCheckpointFingerprints:selected.map { $0.run.end.fingerprint },
                recordEveryPerChain:selected.map(\.recordEvery),
                discardedObservationsPerChain:selected.map(\.discardedObservations),
                acceptanceFraction:Double(accepted)/Double(total),samplingPassed:samplingPassed,statistics:statistics))
        }
        var comparisons: [VivoRingPolymerBeadComparison] = []
        for index in 1..<levels.count {
            let lower=levels[index-1],upper=levels[index]
            let difference=abs(upper.statistics.mean-lower.statistics.mean)
            let combined=sqrt(pow(lower.statistics.monteCarloStandardError,2)+pow(upper.statistics.monteCarloStandardError,2))
            let bound=difference+configuration.confidenceMultiplier*combined
            let tolerance=configuration.absoluteTolerance+configuration.relativeTolerance
                * max(abs(lower.statistics.mean),abs(upper.statistics.mean))
            comparisons.append(.init(lowerBeadCount:lower.beadCount,upperBeadCount:upper.beadCount,
                absoluteMeanDifference:difference,combinedMonteCarloStandardError:combined,
                upperDifferenceBound:bound,allowedDifference:tolerance,
                passed:lower.samplingPassed && upper.samplingPassed && bound <= tolerance))
        }
        let passed=comparisons.suffix(configuration.requiredConsecutiveRefinements).allSatisfy(\.passed)
        let request=VivoRingPolymerConvergenceRequest(
            campaignIdentifier:campaignIdentifier,selectionProtocol:selectionProtocol,
            chains:chains,configuration:configuration)
        let requestID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let request:VivoFingerprint;let levels:[VivoRingPolymerBeadDiagnostics]
            let comparisons:[VivoRingPolymerBeadComparison];let passed:Bool
        }
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            request:requestID,levels:levels,comparisons:comparisons,passed:passed)))
        return .init(campaignIdentifier:campaignIdentifier,potentialDefinitionFingerprint:definitionID,
            temperatureK:temperature,configuration:configuration,beadDiagnostics:levels,
            comparisons:comparisons,passed:passed,requestFingerprint:requestID,evidenceFingerprint:evidenceID,
            selectionProtocol:selectionProtocol,interpretation:interpretation)
    }

    public static func assess(_ request: VivoRingPolymerConvergenceRequest) throws -> VivoRingPolymerConvergenceResult {
        try assess(campaignIdentifier:request.campaignIdentifier,selectionProtocol:request.selectionProtocol,
            chains:request.chains,configuration:request.configuration)
    }
    public static func validate(_ result: VivoRingPolymerConvergenceResult,
                                request: VivoRingPolymerConvergenceRequest) throws {
        guard result == (try assess(request)) else {
            throw VivoChemistryError.invalid("ring-polymer convergence evidence does not reconstruct")
        }
    }
}
