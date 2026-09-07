import Foundation

/// Fixed-cell, distinguishable-nucleus equilibrium sampling. Bead velocities are
/// auxiliary variables at P*T; they are NEVER a physical trajectory or a kappa.
public struct VivoRingPolymerConfiguration: Codable, Sendable, Equatable {
    public var temperatureK: Double
    public var beadCount: Int
    public var timeStepPS: Double
    public var integrationSteps: Int
    public var maximumForceEvaluations: Int
    public var maximumPrimitiveWork: Int
    public init(temperatureK: Double, beadCount: Int = 16, timeStepPS: Double = 0.0001,
                integrationSteps: Int = 8, maximumForceEvaluations: Int = 1_000_000,
                maximumPrimitiveWork: Int = 100_000_000) {
        self.temperatureK = temperatureK; self.beadCount = beadCount; self.timeStepPS = timeStepPS
        self.integrationSteps = integrationSteps; self.maximumForceEvaluations = maximumForceEvaluations
        self.maximumPrimitiveWork = maximumPrimitiveWork
    }
    public func validate() throws {
        guard temperatureK.isFinite, temperatureK > 0, temperatureK <= 100_000,
              (1...256).contains(beadCount), timeStepPS.isFinite, timeStepPS > 0, timeStepPS <= 0.1,
              (1...100_000).contains(integrationSteps), (1...100_000_000).contains(maximumForceEvaluations),
              maximumPrimitiveWork > 0 else { throw VivoChemistryError.invalid("ring-polymer numerical configuration") }
    }
}
public struct VivoRingPolymerCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ring-polymer-equilibrium-checkpoint/v1"
    public let schema: String
    public let definition: VivoNuclearPotentialDefinition
    public let configuration: VivoRingPolymerConfiguration
    public let sweep: UInt64
    public let randomState: VivoSplitMix64
    public let beadPositionsNM: [[VivoVector3D]]
    public let fingerprint: VivoFingerprint

    public init(definition: VivoNuclearPotentialDefinition, configuration: VivoRingPolymerConfiguration,
                seed: UInt64, beadPositionsNM: [[VivoVector3D]]) throws {
        self = try Self.make(definition: definition, configuration: configuration, sweep: 0,
            randomState: .init(state: seed), beadPositionsNM: definition.canonicalBeads(beadPositionsNM))
    }
    static func make(definition: VivoNuclearPotentialDefinition, configuration: VivoRingPolymerConfiguration,
                     sweep: UInt64, randomState: VivoSplitMix64, beadPositionsNM: [[VivoVector3D]]) throws -> Self {
        try definition.validate(); try configuration.validate()
        guard beadPositionsNM.count == configuration.beadCount,
              beadPositionsNM.allSatisfy({ $0.count == definition.atomIndices.count && $0.allSatisfy(\.isFinite) }) else {
            throw VivoChemistryError.invalid("ring-polymer bead layout")
        }
        struct Identity: Encodable {
            let schema: String; let definition: VivoNuclearPotentialDefinition; let configuration: VivoRingPolymerConfiguration
            let sweep: UInt64; let randomState: VivoSplitMix64; let beadPositionsNM: [[VivoVector3D]]
        }
        let identity = Identity(schema: schema, definition: definition, configuration: configuration,
            sweep: sweep, randomState: randomState, beadPositionsNM: beadPositionsNM)
        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity))
        return .init(schema: schema, definition: definition, configuration: configuration,
            sweep: sweep, randomState: randomState, beadPositionsNM: beadPositionsNM, fingerprint: fingerprint)
    }
    private init(schema: String, definition: VivoNuclearPotentialDefinition, configuration: VivoRingPolymerConfiguration,
                 sweep: UInt64, randomState: VivoSplitMix64, beadPositionsNM: [[VivoVector3D]], fingerprint: VivoFingerprint) {
        self.schema = schema; self.definition = definition; self.configuration = configuration
        self.sweep = sweep; self.randomState = randomState; self.beadPositionsNM = beadPositionsNM; self.fingerprint = fingerprint
    }
    public func validate() throws {
        guard self == (try Self.make(definition: definition, configuration: configuration, sweep: sweep,
            randomState: randomState, beadPositionsNM: beadPositionsNM)) else {
            throw VivoChemistryError.invalid("ring-polymer checkpoint payload identity")
        }
    }
}
public struct VivoRingPolymerObservation: Codable, Sendable, Equatable {
    public let sweep: UInt64
    public let accepted: Bool
    public let logAcceptanceProbability: Double
    public let dimensionlessIntegrationWork: Double
    public let beadPositionsNM: [[VivoVector3D]]
    public let centroidPositionsNM: [VivoVector3D]
    public let meanPotentialEnergyKJPerMol: Double
    public let springEnergyKJPerMol: Double
    /// Thermodynamic primitive estimator, not fictitious-bead kinetic energy.
    public let primitiveTotalEnergyKJPerMol: Double
}
public struct VivoRingPolymerRun: Codable, Sendable, Equatable {
    public let start: VivoRingPolymerCheckpoint
    public let end: VivoRingPolymerCheckpoint
    public let observations: [VivoRingPolymerObservation]
    public let forceEvaluations: Int
    public let interpretation: String
}

public enum VivoRingPolymerSampling {
    public static let interpretation = "Primitive imaginary-time Boltzmann path integral for distinguishable nuclei, sampled by Metropolized normal-mode split HMC at bead temperature P*T. Fixed cell and fixed atom/mass mapping; no holonomic constraints, nuclear exchange, real-time RPMD or instanton rate. Bead-number, sampling and numerical-potential convergence remain separate. Retained observations include rejected states and are correlated, not independent samples."
    static let rtPerK = VivoAtomicUnits.gasConstantJPerMolK/1000
    // hbar N_A in kJ mol^-1 ps; Da*(nm/ps)^2 = kJ/mol.
    static let hbar = VivoAtomicUnits.planckJS*6.02214076e23*1e9/(2*Double.pi)

    /// Full H_P spring energy. The configurational target is exp[-H_P/(P RT)].
    public static func springEnergy(positions: [[VivoVector3D]], massesDa: [Double], temperatureK: Double) throws -> Double {
        let p = positions.count
        guard p > 0, !massesDa.isEmpty, massesDa.allSatisfy({ $0.isFinite && $0 > 0 }),
              positions.allSatisfy({ $0.count == massesDa.count && $0.allSatisfy(\.isFinite) }),
              temperatureK.isFinite, temperatureK > 0 else { throw VivoChemistryError.invalid("ring spring geometry/masses") }
        let omega = Double(p)*rtPerK*temperatureK/hbar
        var energy = 0.0
        for b in 0..<p { for i in massesDa.indices {
            energy += 0.5*massesDa[i]*omega*omega*(positions[b][i]-positions[(b+1)%p][i]).squaredNorm
        } }
        guard energy.isFinite else { throw VivoChemistryError.convergence("ring spring overflow") }
        return energy
    }
    /// The mass-dependent path-integral normalization must be retained for
    /// isotope reweighting. Common temperature/bead constants are omitted only
    /// because this function is restricted to comparisons at the SAME P and T.
    public static func reducedConfigurationalPotential(positions: [[VivoVector3D]], massesDa: [Double],
        temperatureK: Double, beadPotentialEnergiesKJPerMol: [Double]) throws -> Double {
        guard beadPotentialEnergiesKJPerMol.count == positions.count,
              beadPotentialEnergiesKJPerMol.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("ring potential energy count") }
        let p = Double(positions.count), rt = rtPerK*temperatureK
        let spring = try springEnergy(positions: positions,massesDa: massesDa,temperatureK: temperatureK)
        let value = (spring+beadPotentialEnergiesKJPerMol.reduce(0,+))/(p*rt)
            - 1.5*p*massesDa.reduce(0) { $0+log($1) }
        guard value.isFinite else { throw VivoChemistryError.convergence("ring reduced action overflow") }
        return value
    }
    public static func run(potential: VivoNuclearPotential, checkpoint: VivoRingPolymerCheckpoint,
                           sweeps: Int, recordEvery: Int = 1) async throws -> VivoRingPolymerRun {
        try checkpoint.validate()
        let cfg = checkpoint.configuration, p = cfg.beadCount, n = checkpoint.definition.massesDa.count
        guard checkpoint.definition == potential.definition, (1...1_000_000).contains(sweeps),
              (1...sweeps).contains(recordEvery), checkpoint.sweep <= UInt64.max-UInt64(sweeps) else {
            throw VivoChemistryError.invalid("ring-polymer continuation identity or run size")
        }
        let calls = Double(p)*(1+Double(sweeps)*Double(cfg.integrationSteps))
        let work = Double(sweeps)*Double(cfg.integrationSteps)*Double(p*p)*Double(n)*12
        guard calls <= Double(cfg.maximumForceEvaluations), work <= Double(cfg.maximumPrimitiveWork),
              Double(sweeps/recordEvery+1)*Double(p)*Double(n)*3 <= Double(cfg.maximumPrimitiveWork) else {
            throw VivoChemistryError.resourceLimit("ring-polymer aggregate force, transform or observation budget")
        }
        let rt = rtPerK*cfg.temperatureK, beadRT = Double(p)*rt, omegaP = beadRT/hbar
        let masses = potential.definition.massesDa
        var transform = [[Double]](repeating: [Double](repeating: 0,count: p),count: p)
        var frequencies = [Double](repeating: 0,count: p)
        for b in 0..<p { transform[0][b] = 1/sqrt(Double(p)) }
        var row = 1
        if p >= 3 {
            for k in 1...((p-1)/2) {
                for b in 0..<p {
                    let angle = 2*Double.pi*Double(k*b)/Double(p)
                    transform[row][b] = sqrt(2/Double(p))*cos(angle)
                    transform[row+1][b] = sqrt(2/Double(p))*sin(angle)
                }
                frequencies[row] = 2*omegaP*sin(Double.pi*Double(k)/Double(p))
                frequencies[row+1] = frequencies[row]; row += 2
            }
        }
        if p % 2 == 0 {
            for b in 0..<p { transform[p-1][b] = (b%2 == 0 ? 1 : -1)/sqrt(Double(p)) }
            frequencies[p-1] = 2*omegaP
        }
        func freeStep(_ q: inout [[VivoVector3D]], _ v: inout [[VivoVector3D]]) {
            var modesQ = q, modesV = v
            for k in 0..<p { for i in 0..<n {
                var a = VivoVector3D.zero, b = VivoVector3D.zero
                for bead in 0..<p { a = a+q[bead][i]*transform[k][bead]; b = b+v[bead][i]*transform[k][bead] }
                let w = frequencies[k], angle = w*cfg.timeStepPS
                if k == 0 { modesQ[k][i] = a+b*cfg.timeStepPS; modesV[k][i] = b }
                else { modesQ[k][i] = a*cos(angle)+b*(sin(angle)/w); modesV[k][i] = b*cos(angle)-a*(w*sin(angle)) }
            } }
            for bead in 0..<p { for i in 0..<n {
                var a = VivoVector3D.zero, b = VivoVector3D.zero
                for k in 0..<p { a = a+modesQ[k][i]*transform[k][bead]; b = b+modesV[k][i]*transform[k][bead] }
                q[bead][i] = a; v[bead][i] = b
            } }
        }
        func kinetic(_ v: [[VivoVector3D]]) -> Double {
            v.reduce(0.0) { sum, bead in sum+zip(bead,masses).reduce(0) { $0+0.5*$1.1*$1.0.squaredNorm } }
        }
        var count = 0
        func evaluate(_ q: [[VivoVector3D]]) async throws -> [VivoNuclearPotentialEvaluation] {
            var values: [VivoNuclearPotentialEvaluation] = []; values.reserveCapacity(p)
            for bead in q {
                guard count < cfg.maximumForceEvaluations else { throw VivoChemistryError.resourceLimit("ring force budget") }
                count += 1; values.append(try await potential.checked(bead))
            }
            return values
        }
        var q = checkpoint.beadPositionsNM, rng = checkpoint.randomState
        var current = try await evaluate(q), records: [VivoRingPolymerObservation] = []
        for sweep in 1...sweeps {
            try Task.checkCancellation()
            var candidate = q, force = current
            var velocity = q.map { _ in masses.map { mass in
                let s = sqrt(beadRT/mass); return VivoVector3D(rng.normal()*s,rng.normal()*s,rng.normal()*s)
            } }
            let initialH = current.reduce(0) { $0+$1.energyKJPerMol }+kinetic(velocity)
                + (try springEnergy(positions: q,massesDa: masses,temperatureK: cfg.temperatureK))
            for _ in 0..<cfg.integrationSteps {
                for b in 0..<p { for i in 0..<n { velocity[b][i] = velocity[b][i]+force[b].forcesKJPerMolNM[i]*(0.5*cfg.timeStepPS/masses[i]) } }
                freeStep(&candidate,&velocity)
                candidate = try potential.definition.canonicalBeads(candidate)
                force = try await evaluate(candidate)
                for b in 0..<p { for i in 0..<n { velocity[b][i] = velocity[b][i]+force[b].forcesKJPerMolNM[i]*(0.5*cfg.timeStepPS/masses[i]) } }
            }
            let finalH = force.reduce(0) { $0+$1.energyKJPerMol }+kinetic(velocity)
                + (try springEnergy(positions: candidate,massesDa: masses,temperatureK: cfg.temperatureK))
            let work = (finalH-initialH)/beadRT
            guard work.isFinite else { throw VivoChemistryError.convergence("ring-polymer integration energy overflow") }
            let logA = min(0,-work), accepted = log(max(rng.unitInterval(),Double.leastNonzeroMagnitude)) < logA
            if accepted { q = candidate; current = force }
            let index = checkpoint.sweep+UInt64(sweep)
            if index % UInt64(recordEvery) == 0 {
                let spring = try springEnergy(positions: q,massesDa: masses,temperatureK: cfg.temperatureK)
                let mean = current.reduce(0) { $0+$1.energyKJPerMol }/Double(p)
                let centroid = (0..<n).map { i in q.reduce(VivoVector3D.zero) { $0+$1[i] }/Double(p) }
                records.append(.init(sweep: index,accepted: accepted,logAcceptanceProbability: logA,
                    dimensionlessIntegrationWork: work,beadPositionsNM: q,centroidPositionsNM: centroid,
                    meanPotentialEnergyKJPerMol: mean,springEnergyKJPerMol: spring,
                    primitiveTotalEnergyKJPerMol: 1.5*Double(n*p)*rt-spring/Double(p)+mean))
            }
        }
        let end = try VivoRingPolymerCheckpoint.make(definition: potential.definition,configuration: cfg,
            sweep: checkpoint.sweep+UInt64(sweeps),randomState: rng,beadPositionsNM: q)
        return .init(start: checkpoint,end: end,observations: records,forceEvaluations: count,interpretation: interpretation)
    }
}
