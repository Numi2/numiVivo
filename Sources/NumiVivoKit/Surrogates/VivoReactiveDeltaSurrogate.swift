import Foundation

public struct VivoReactiveDistanceFeature: Codable, Sendable, Equatable {
    /// Slots in the explicitly bound independent-atom layout, NOT particle IDs.
    public let atomA: Int
    public let atomB: Int
    public let scaleNM: Double
    public init(atomA: Int, atomB: Int, scaleNM: Double) { self.atomA = atomA; self.atomB = atomB; self.scaleNM = scaleNM }
}
public struct VivoReactiveTrainingLabel: Codable, Sendable, Equatable {
    public let identifier: String
    public let sourceGroup: String
    public let authority: VivoNuclearPotentialEvaluation
    public let baseline: VivoNuclearPotentialEvaluation
    public init(identifier: String, sourceGroup: String, authority: VivoNuclearPotentialEvaluation, baseline: VivoNuclearPotentialEvaluation) {
        self.identifier = identifier; self.sourceGroup = sourceGroup; self.authority = authority; self.baseline = baseline
    }
}
public struct VivoReactiveSurrogateConfiguration: Codable, Sendable, Equatable {
    public var features: [VivoReactiveDistanceFeature]
    public var maximumCenters: Int
    public var committeeSize: Int
    public var ridge: Double
    public var energyFitScaleKJPerMol: Double
    public var forceFitScaleKJPerMolNM: Double
    public var maximumHeldOutEnergyErrorKJPerMol: Double
    public var maximumHeldOutForceErrorKJPerMolNM: Double
    public var maximumEnergyDisagreementKJPerMol: Double
    public var maximumForceDisagreementKJPerMolNM: Double
    public var maximumNormalizedDistance: Double
    public var seed: UInt64
    public var maximumPrimitiveWork: Int
    public init(features: [VivoReactiveDistanceFeature], maximumCenters: Int = 32, committeeSize: Int = 4,
                ridge: Double = 1e-8, energyFitScaleKJPerMol: Double = 1, forceFitScaleKJPerMolNM: Double = 10,
                maximumHeldOutEnergyErrorKJPerMol: Double = 1, maximumHeldOutForceErrorKJPerMolNM: Double = 10,
                maximumEnergyDisagreementKJPerMol: Double = 1, maximumForceDisagreementKJPerMolNM: Double = 10,
                maximumNormalizedDistance: Double = 2, seed: UInt64 = 0, maximumPrimitiveWork: Int = 100_000_000) {
        self.features = features; self.maximumCenters = maximumCenters; self.committeeSize = committeeSize; self.ridge = ridge
        self.energyFitScaleKJPerMol = energyFitScaleKJPerMol; self.forceFitScaleKJPerMolNM = forceFitScaleKJPerMolNM
        self.maximumHeldOutEnergyErrorKJPerMol = maximumHeldOutEnergyErrorKJPerMol
        self.maximumHeldOutForceErrorKJPerMolNM = maximumHeldOutForceErrorKJPerMolNM
        self.maximumEnergyDisagreementKJPerMol = maximumEnergyDisagreementKJPerMol
        self.maximumForceDisagreementKJPerMolNM = maximumForceDisagreementKJPerMolNM
        self.maximumNormalizedDistance = maximumNormalizedDistance; self.seed = seed; self.maximumPrimitiveWork = maximumPrimitiveWork
    }
    public func validate(atomCount: Int) throws {
        let keys = features.map { "\(min($0.atomA,$0.atomB)):\(max($0.atomA,$0.atomB))" }
        guard (1...256).contains(features.count), Set(keys).count == keys.count,
              features.allSatisfy({ $0.atomA >= 0 && $0.atomA < atomCount && $0.atomB >= 0 && $0.atomB < atomCount
                && $0.atomA != $0.atomB && $0.scaleNM.isFinite && $0.scaleNM >= 1e-6 && $0.scaleNM <= 100 }),
              (1...128).contains(maximumCenters), (2...8).contains(committeeSize), ridge.isFinite, ridge >= 1e-12, ridge <= 1,
              [energyFitScaleKJPerMol,forceFitScaleKJPerMolNM,maximumHeldOutEnergyErrorKJPerMol,
               maximumHeldOutForceErrorKJPerMolNM,maximumEnergyDisagreementKJPerMol,
               maximumForceDisagreementKJPerMolNM,maximumNormalizedDistance].allSatisfy({ $0.isFinite && $0 > 0 }),
              maximumPrimitiveWork > 0 else { throw VivoChemistryError.invalid("reactive surrogate features, fit scales or limits") }
    }
}
public struct VivoReactiveSurrogateQualification: Codable, Sendable, Equatable {
    public let heldOutGroups: [String]
    public let trainingGroupCount: Int
    public let heldOutCount: Int
    public let maximumEnergyErrorKJPerMol: Double
    public let maximumForceErrorKJPerMolNM: Double
    public let passed: Bool
}
public struct VivoReactiveSurrogatePayload: Codable, Sendable, Equatable {
    public let schema: String
    public let authorityDefinition: VivoNuclearPotentialDefinition
    public let baselineDefinition: VivoNuclearPotentialDefinition
    public let configuration: VivoReactiveSurrogateConfiguration
    public let trainingDataFingerprint: VivoFingerprint
    public let centersNM: [[Double]]
    /// Constant followed by radial-basis coefficients, in kJ/mol.
    public let coefficientsKJPerMol: [[Double]]
    public let qualification: VivoReactiveSurrogateQualification
}
public struct VivoReactiveSurrogateModel: Codable, Sendable, Equatable {
    public let payload: VivoReactiveSurrogatePayload
    public let fingerprint: VivoFingerprint
    init(payload: VivoReactiveSurrogatePayload) throws {
        self.payload = payload; fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(payload))
        try validate()
    }
    public func validate() throws {
        let p = payload, n = p.authorityDefinition.atomIndices.count
        try p.authorityDefinition.validate(); try p.baselineDefinition.validate(); try p.configuration.validate(atomCount: n)
        guard p.schema == "numivivo.org/reactive-delta-rbf/v1", fingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(p))),
              p.authorityDefinition.atomIndices == p.baselineDefinition.atomIndices,
              p.authorityDefinition.particleIndices == p.baselineDefinition.particleIndices,
              p.authorityDefinition.massesDa == p.baselineDefinition.massesDa,
              p.authorityDefinition.periodicCell == p.baselineDefinition.periodicCell,
              p.authorityDefinition.periodicMoleculeGroups == p.baselineDefinition.periodicMoleculeGroups,
              !p.centersNM.isEmpty, p.centersNM.count <= p.configuration.maximumCenters,
              p.centersNM.allSatisfy({ $0.count == p.configuration.features.count && $0.allSatisfy({ $0.isFinite && $0 > 0 }) }),
              p.coefficientsKJPerMol.count == p.configuration.committeeSize,
              p.coefficientsKJPerMol.allSatisfy({ $0.count == p.centersNM.count+1 && $0.allSatisfy(\.isFinite) }),
              p.qualification.trainingGroupCount >= 2, p.qualification.heldOutCount > 0,
              !p.qualification.heldOutGroups.isEmpty,
              [p.qualification.maximumEnergyErrorKJPerMol,p.qualification.maximumForceErrorKJPerMolNM].allSatisfy({ $0.isFinite && $0 >= 0 }),
              p.qualification.passed == (p.qualification.maximumEnergyErrorKJPerMol <= p.configuration.maximumHeldOutEnergyErrorKJPerMol
                && p.qualification.maximumForceErrorKJPerMolNM <= p.configuration.maximumHeldOutForceErrorKJPerMolNM) else {
            throw VivoChemistryError.invalid("reactive surrogate payload, mapping or qualification integrity")
        }
    }
}
public struct VivoReactiveSurrogatePrediction: Codable, Sendable, Equatable {
    public let modelFingerprint: VivoFingerprint
    public let positionsNM: [VivoVector3D]
    public let deltaEnergyKJPerMol: Double
    public let deltaForcesKJPerMolNM: [VivoVector3D]
    public let energyDisagreementKJPerMol: Double
    public let maximumForceDisagreementKJPerMolNM: Double
    public let nearestNormalizedDistance: Double
    public let eligible: Bool
    public let reasons: [String]
}

/// Energy-conserving local delta model: forces are exact derivatives of the
/// learned scalar energy. Fixed atom identities; no distance-based bond changes,
/// transferability claim or substitution for periodic electrostatics.
public enum VivoReactiveDeltaSurrogate {
    struct Geometry {
        let distances: [Double]
        let directions: [VivoVector3D]
    }
    static func geometry(_ positions: [VivoVector3D], features: [VivoReactiveDistanceFeature]) throws -> Geometry {
        guard positions.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("reactive descriptor geometry") }
        var distances: [Double] = [], directions: [VivoVector3D] = []
        for feature in features {
            guard positions.indices.contains(feature.atomA), positions.indices.contains(feature.atomB) else { throw VivoChemistryError.invalid("reactive descriptor atom slot") }
            // Deliberately no atomwise periodic reimaging. Inputs use one mapped
            // unwrapped chart; image transfer outside training requests new labels.
            let delta = positions[feature.atomA]-positions[feature.atomB], r = delta.norm
            guard r.isFinite, r > 1e-10 else { throw VivoChemistryError.invalid("degenerate reactive distance feature") }
            distances.append(r); directions.append(delta/r)
        }
        return .init(distances: distances,directions: directions)
    }
    static func basis(_ g: Geometry, features: [VivoReactiveDistanceFeature], centers: [[Double]], atomCount: Int)
    -> (energy: [Double], forces: [[VivoVector3D]], nearest: Double) {
        var energy = [Double](repeating: 1,count: centers.count+1)
        var forces = [[VivoVector3D]](repeating: [VivoVector3D](repeating: .zero,count: atomCount),count: centers.count+1)
        var nearest = Double.infinity
        for (c,center) in centers.enumerated() {
            var distance2 = 0.0
            for j in features.indices { distance2 += pow((g.distances[j]-center[j])/features[j].scaleNM,2) }
            nearest = min(nearest,sqrt(distance2))
            let phi = exp(-0.5*distance2); energy[c+1] = phi
            for j in features.indices {
                let f = features[j]
                let value = g.directions[j]*(phi*(g.distances[j]-center[j])/(f.scaleNM*f.scaleNM))
                forces[c+1][f.atomA] = forces[c+1][f.atomA]+value
                forces[c+1][f.atomB] = forces[c+1][f.atomB]-value
            }
        }
        return (energy,forces,nearest)
    }
    static func members(geometry g: Geometry, payload: VivoReactiveSurrogatePayload) -> ([Double],[[VivoVector3D]],Double) {
        let n = payload.authorityDefinition.atomIndices.count
        let b = basis(g,features: payload.configuration.features,centers: payload.centersNM,atomCount: n)
        var energies: [Double] = [], forces: [[VivoVector3D]] = []
        for coefficients in payload.coefficientsKJPerMol {
            energies.append(zip(b.energy,coefficients).reduce(0) { $0+$1.0*$1.1 })
            forces.append((0..<n).map { i in coefficients.indices.reduce(VivoVector3D.zero) { $0+b.forces[$1][i]*coefficients[$1] } })
        }
        return (energies,forces,b.nearest)
    }
    public static func predict(_ model: VivoReactiveSurrogateModel, positionsNM: [VivoVector3D]) throws -> VivoReactiveSurrogatePrediction {
        try model.validate()
        guard positionsNM.count == model.payload.authorityDefinition.atomIndices.count else { throw VivoChemistryError.invalid("reactive model atom count") }
        let g = try geometry(positionsNM,features: model.payload.configuration.features)
        let (energies,forces,nearest) = members(geometry: g,payload: model.payload)
        return try summarize(model: model,positions: positionsNM,energies: energies,forces: forces,nearest: nearest)
    }
    static func summarize(model: VivoReactiveSurrogateModel, positions: [VivoVector3D], energies: [Double],
                          forces: [[VivoVector3D]], nearest: Double) throws -> VivoReactiveSurrogatePrediction {
        let k = Double(energies.count), n = positions.count
        let energy = energies.reduce(0,+)/k
        let force = (0..<n).map { i in forces.reduce(VivoVector3D.zero) { $0+$1[i] }/k }
        let sd = sqrt(energies.reduce(0) { $0+pow($1-energy,2) }/(k-1))
        let fsd = (0..<n).map { i in sqrt(forces.reduce(0) { $0+($1[i]-force[i]).squaredNorm }/(k-1)) }.max() ?? 0
        guard energy.isFinite, force.allSatisfy(\.isFinite), sd.isFinite, fsd.isFinite, nearest.isFinite else {
            throw VivoChemistryError.convergence("reactive prediction overflow")
        }
        let cfg = model.payload.configuration
        var reasons: [String] = []
        if !model.payload.qualification.passed { reasons.append("held-out qualification failed") }
        if nearest > cfg.maximumNormalizedDistance { reasons.append("outside training-descriptor support") }
        if sd > cfg.maximumEnergyDisagreementKJPerMol || fsd > cfg.maximumForceDisagreementKJPerMolNM {
            reasons.append("committee disagreement requires authority labels")
        }
        return .init(modelFingerprint: model.fingerprint,positionsNM: positions,deltaEnergyKJPerMol: energy,
            deltaForcesKJPerMolNM: force,energyDisagreementKJPerMol: sd,maximumForceDisagreementKJPerMolNM: fsd,
            nearestNormalizedDistance: nearest,eligible: reasons.isEmpty,reasons: reasons)
    }

    public static func train(authority: VivoNuclearPotentialDefinition, baseline: VivoNuclearPotentialDefinition,
                             labels: [VivoReactiveTrainingLabel], heldOutGroups: [String],
                             configuration cfg: VivoReactiveSurrogateConfiguration) throws -> VivoReactiveSurrogateModel {
        try authority.validate(); try baseline.validate(); try cfg.validate(atomCount: authority.atomIndices.count)
        let n = authority.atomIndices.count
        guard authority.atomIndices == baseline.atomIndices, authority.particleIndices == baseline.particleIndices,
              authority.massesDa == baseline.massesDa, authority.periodicCell == baseline.periodicCell,
              authority.periodicMoleculeGroups == baseline.periodicMoleculeGroups,
              (3...4096).contains(labels.count), Set(labels.map(\.identifier)).count == labels.count,
              labels.allSatisfy({ !$0.identifier.isEmpty && !$0.sourceGroup.isEmpty }), !heldOutGroups.isEmpty,
              Set(heldOutGroups).count == heldOutGroups.count,
              Set(heldOutGroups).isSubset(of: Set(labels.map(\.sourceGroup))) else {
            throw VivoChemistryError.invalid("reactive training identity or disjoint group split")
        }
        var geometries: [Geometry] = [], hashes = Set<VivoFingerprint>()
        for label in labels {
            let positions = label.authority.requestedPositionsNM
            try label.authority.validate(definition: authority,positionsNM: positions)
            try label.baseline.validate(definition: baseline,positionsNM: positions)
            guard hashes.insert(try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(positions))).inserted else {
                throw VivoChemistryError.invalid("duplicate geometry in reactive training/held-out data")
            }
            geometries.append(try geometry(positions,features: cfg.features))
        }
        let training = labels.indices.filter { !heldOutGroups.contains(labels[$0].sourceGroup) }.sorted { labels[$0].identifier < labels[$1].identifier }
        let held = labels.indices.filter { heldOutGroups.contains(labels[$0].sourceGroup) }
        let groups = Set(training.map { labels[$0].sourceGroup }).sorted()
        guard training.count >= 2, groups.count >= 2, !held.isEmpty else { throw VivoChemistryError.invalid("reactive model requires at least two training groups and held-out data") }
        let centerCount = min(cfg.maximumCenters,training.count), width = centerCount+1
        let work = Double(cfg.committeeSize)*Double(training.count)*Double(1+3*n)*Double(width*width)
        guard work <= Double(cfg.maximumPrimitiveWork) else { throw VivoChemistryError.resourceLimit("reactive energy/force fitting budget") }
        var selected = [training[0]], nearest2 = [Double](repeating: .infinity,count: labels.count)
        while selected.count < centerCount {
            let latest = selected.last!
            for i in training {
                let d = cfg.features.indices.reduce(0.0) { $0+pow((geometries[i].distances[$1]-geometries[latest].distances[$1])/cfg.features[$1].scaleNM,2) }
                nearest2[i] = min(nearest2[i],d)
            }
            let candidates = training.filter { !selected.contains($0) }
            selected.append(candidates.max { nearest2[$0] < nearest2[$1] }!)
        }
        let centers = selected.map { geometries[$0].distances }
        let bases = training.map { basis(geometries[$0],features: cfg.features,centers: centers,atomCount: n) }
        let referenceEnergy = training.reduce(0) { $0+labels[$1].authority.energyKJPerMol-labels[$1].baseline.energyKJPerMol }/Double(training.count)
        var coefficients: [[Double]] = [], rng = VivoSplitMix64(state: cfg.seed)
        for member in 0..<cfg.committeeSize {
            var multiplicity = Dictionary(uniqueKeysWithValues: groups.map { ($0,member == 0 ? 1 : 0) })
            if member > 0 {
                for _ in groups { multiplicity[groups[Int(rng.next()%UInt64(groups.count))],default: 0] += 1 }
            }
            var normal = VivoQMMatrix(width,width), rhs = [Double](repeating: 0,count: width)
            func add(_ row: [Double], _ target: Double, _ weight: Double) {
                for i in 0..<width {
                    rhs[i] += weight*row[i]*target
                    for j in 0...i { normal[i,j] += weight*row[i]*row[j] }
                }
            }
            for (slot,index) in training.enumerated() {
                let copies = Double(multiplicity[labels[index].sourceGroup]!)
                if copies == 0 { continue }
                let label = labels[index], b = bases[slot]
                add(b.energy,label.authority.energyKJPerMol-label.baseline.energyKJPerMol-referenceEnergy,
                    copies/(cfg.energyFitScaleKJPerMol*cfg.energyFitScaleKJPerMol))
                for i in 0..<n {
                    let target = label.authority.forcesKJPerMolNM[i]-label.baseline.forcesKJPerMolNM[i]
                    let w = copies/(cfg.forceFitScaleKJPerMolNM*cfg.forceFitScaleKJPerMolNM)
                    add(b.forces.map { $0[i].x },target.x,w)
                    add(b.forces.map { $0[i].y },target.y,w)
                    add(b.forces.map { $0[i].z },target.z,w)
                }
            }
            for i in 0..<width {
                if i > 0 { normal[i,i] += cfg.ridge }
                for j in 0..<i { normal[j,i] = normal[i,j] }
            }
            var beta = try VivoQMDenseAlgebra.solve(normal,rhs: rhs)
            beta[0] += referenceEnergy; coefficients.append(beta)
        }
        let dataID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(labels))
        func payload(_ qualification: VivoReactiveSurrogateQualification) -> VivoReactiveSurrogatePayload {
            .init(schema: "numivivo.org/reactive-delta-rbf/v1",authorityDefinition: authority,baselineDefinition: baseline,
                configuration: cfg,trainingDataFingerprint: dataID,centersNM: centers,coefficientsKJPerMol: coefficients,qualification: qualification)
        }
        let preliminary = payload(.init(heldOutGroups: heldOutGroups.sorted(),trainingGroupCount: groups.count,
            heldOutCount: held.count,maximumEnergyErrorKJPerMol: 0,maximumForceErrorKJPerMolNM: 0,passed: true))
        var maxE = 0.0, maxF = 0.0
        for i in held {
            let (e,f,_) = members(geometry: geometries[i],payload: preliminary), k = Double(cfg.committeeSize)
            maxE = max(maxE,abs(e.reduce(0,+)/k-(labels[i].authority.energyKJPerMol-labels[i].baseline.energyKJPerMol)))
            for a in 0..<n {
                let prediction = f.reduce(VivoVector3D.zero) { $0+$1[a] }/k
                maxF = max(maxF,(prediction-(labels[i].authority.forcesKJPerMolNM[a]-labels[i].baseline.forcesKJPerMolNM[a])).norm)
            }
        }
        return try .init(payload: payload(.init(heldOutGroups: heldOutGroups.sorted(),trainingGroupCount: groups.count,
            heldOutCount: held.count,maximumEnergyErrorKJPerMol: maxE,maximumForceErrorKJPerMolNM: maxF,
            passed: maxE <= cfg.maximumHeldOutEnergyErrorKJPerMol && maxF <= cfg.maximumHeldOutForceErrorKJPerMolNM)))
    }
}
