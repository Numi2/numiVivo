import Foundation

/// Forces are minus energy derivatives, NOT gradients. An embedding backend must
/// supply reactions on ALL point-charge centers, including virtual-site centers.
public struct VivoQMMMElectronicForceResult: Codable, Sendable, Equatable {
    public static let convention = "E-electronic+nuclear+nuclear-MM; no-MM-self; no-QM-MM-LJ"
    public let electronicSystemFingerprint: VivoFingerprint
    public let energyConvention: String
    public let derivativeMethod: String
    public let energyHartree: Double
    public let nucleusForcesHartreePerBohr: [VivoVector3D]
    public let pointChargeForcesHartreePerBohr: [VivoVector3D]

    public init(system: VivoElectronicSystem, energyHartree: Double,
                nucleusForcesHartreePerBohr: [VivoVector3D], pointChargeForcesHartreePerBohr: [VivoVector3D],
                derivativeMethod: String) throws {
        electronicSystemFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(system))
        energyConvention = Self.convention; self.derivativeMethod = derivativeMethod
        self.energyHartree = energyHartree; self.nucleusForcesHartreePerBohr = nucleusForcesHartreePerBohr
        self.pointChargeForcesHartreePerBohr = pointChargeForcesHartreePerBohr
        try validate(system: system)
    }
    public func validate(system: VivoElectronicSystem) throws {
        try system.validate()
        guard energyConvention == Self.convention, !derivativeMethod.isEmpty, energyHartree.isFinite,
              electronicSystemFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(system))),
              nucleusForcesHartreePerBohr.count == system.nuclei.count,
              pointChargeForcesHartreePerBohr.count == system.pointCharges.count,
              nucleusForcesHartreePerBohr.allSatisfy(\.isFinite), pointChargeForcesHartreePerBohr.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("QM/MM electronic force identity, convention, missing centers or nonfinite forces")
        }
    }
}

public struct VivoQMMMLinkForceContribution: Codable, Sendable, Equatable {
    public let qmForceHartreePerBohr: VivoVector3D
    public let mmForceHartreePerBohr: VivoVector3D
}

public enum VivoQMMMLinkForceProjection {
    /// rL = rQ + d (rM-rQ)/|rM-rQ|, with FIXED d, as in VivoQMMMCompiler.
    /// JM = (d/r)(I-u*uT), JQ = I-JM; physical forces are J^T FL.
    /// Fixed-distance projection is not constant-ratio force splitting.
    public static func project(qmPositionBohr q: VivoVector3D, mmPositionBohr m: VivoVector3D,
                               linkDistanceBohr d: Double, linkForceHartreePerBohr f: VivoVector3D) throws -> VivoQMMMLinkForceContribution {
        let delta = m - q, r = delta.norm
        guard q.isFinite, m.isFinite, f.isFinite, d.isFinite, d > 0, r.isFinite, r > d else {
            throw VivoChemistryError.invalid("QM/MM link force projection geometry")
        }
        let u = delta / r
        let mm = (f - u * u.dot(f)) * (d / r)
        let qm = f - mm
        guard qm.isFinite, mm.isFinite else { throw VivoChemistryError.invalid("QM/MM link projection overflow") }
        return .init(qmForceHartreePerBohr: qm, mmForceHartreePerBohr: mm)
    }
}

public struct VivoQMMMForceResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-forces/v1"
    public let schema: String
    public let sourceFrameFingerprint: VivoFingerprint
    public let sourceCheckpointFingerprint: VivoFingerprint?
    public let finiteClusterFingerprint: VivoFingerprint
    public let preparedSnapshotFingerprint: VivoFingerprint
    public let energyConvention: String
    public let electronicDerivativeMethod: String
    public let energyHartree: Double
    public let electronicEnergyHartree: Double
    public let lennardJonesEnergyHartree: Double
    public let atomToParticle: [UInt32]
    /// Physical forces in the ORIGINAL MD particle layout; virtual slots are zero.
    /// All link, embedding-center and LJ virtual-site forces are folded exactly once.
    public let particleForcesHartreePerBohr: [VivoVector3D]
    public let netForceHartreePerBohr: VivoVector3D
    public let netTorqueHartree: VivoVector3D
    public var atomForcesHartreePerBohr: [VivoVector3D] {
        atomToParticle.map { particleForcesHartreePerBohr[Int($0)] }
    }
    public var particleForcesKJPerMolPerNM: [VivoVector3D] {
        let factor = VivoAtomicUnits.hartreeInKJPerMol / VivoAtomicUnits.bohrInNM
        return particleForcesHartreePerBohr.map { $0 * factor }
    }
}

public enum VivoQMMMForceMapper {
    public static func assemble(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                frame: VivoMDQMMMPreparedFrame, electronic: VivoQMMMElectronicForceResult,
                                budget: VivoChemistryBudget = .init()) throws -> VivoQMMMForceResult {
        try frame.validate(document: document, system: system)
        try electronic.validate(system: frame.region.electronicSystem)
        let region = frame.region, positions = frame.cluster.particlePositionsNM, map = frame.cluster.atomToParticle
        let lj = try VivoQMMMCompiler.lennardJones(document: document, system: system,
            particlePositionsNM: positions, request: frame.request, budget: budget)
        guard var forces = lj.particleForcesHartreePerBohr, forces.count == system.particles.count else {
            throw VivoChemistryError.invalid("QM/MM force assembly requires a full-particle LJ result")
        }
        let projected = try projectCenters(electronicSystem: region.electronicSystem,links: region.links,
            atomToParticle: map,particlePositionsNM: positions,nucleusForces: electronic.nucleusForcesHartreePerBohr,
            pointChargeForces: electronic.pointChargeForcesHartreePerBohr)
        for i in forces.indices { forces[i] = forces[i]+projected[i] }
        let topology = try VivoQMMMTopology(document: document, system: system)
        forces = try topology.siteGraph.redistribute(rawForces: forces,state: topology.siteGraph.construct(positionsNM: positions))
        let energy = electronic.energyHartree + lj.energyHartree
        guard energy.isFinite, forces.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite assembled QM/MM result") }
        let origin = map.reduce(VivoVector3D.zero) { $0 + positions[Int($1)] } / Double(map.count)
        var net = VivoVector3D.zero, torque = VivoVector3D.zero
        for particle in map {
            let i = Int(particle)
            net = net + forces[i]
            torque = torque + ((positions[i] - origin) / VivoAtomicUnits.bohrInNM).cross(forces[i])
        }
        guard net.isFinite, torque.isFinite else { throw VivoChemistryError.invalid("QM/MM force diagnostics overflow") }
        return .init(schema: VivoQMMMForceResult.schema, sourceFrameFingerprint: frame.cluster.sourceFrameFingerprint,
            sourceCheckpointFingerprint: frame.sourceCheckpointFingerprint,
            finiteClusterFingerprint: try frame.cluster.fingerprint(), preparedSnapshotFingerprint: region.snapshotFingerprint,
            energyConvention: region.energyConvention, electronicDerivativeMethod: electronic.derivativeMethod, energyHartree: energy, electronicEnergyHartree: electronic.energyHartree,
            lennardJonesEnergyHartree: lj.energyHartree, atomToParticle: map, particleForcesHartreePerBohr: forces,
            netForceHartreePerBohr: net, netTorqueHartree: torque)
    }

    /// Backend-independent entry point. Native solvers and external adapters share
    /// the same geometry, identity, units, force sign and redistribution contract.
    public static func evaluate(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                frame: VivoMDQMMMPreparedFrame, budget: VivoChemistryBudget = .init(),
                                electronicEvaluator: (VivoElectronicSystem) throws -> VivoQMMMElectronicForceResult) throws -> VivoQMMMForceResult {
        try frame.validate(document: document, system: system)
        return try assemble(document: document, system: system, frame: frame,
            electronic: electronicEvaluator(frame.region.electronicSystem), budget: budget)
    }
}

/// A bounded reference adapter for energy-only backends. This differentiates all
/// independent electronic centers, including MM charges; the mapper subsequently
/// applies link and virtual-site Jacobians. It is not an analytic-force claim.
public enum VivoQMMMElectronicFiniteDifferences {
    public static func evaluate(system: VivoElectronicSystem, stepBohr: Double = 1e-4,
                                maximumEnergyEvaluations: Int = 4096,
                                energy: (VivoElectronicSystem) throws -> Double) throws -> VivoQMMMElectronicForceResult {
        try system.validate()
        let centers = system.nuclei.count + system.pointCharges.count
        let (calls, overflow) = centers.multipliedReportingOverflow(by: 6)
        guard stepBohr.isFinite, stepBohr >= 1e-8, stepBohr <= 0.01,
              !overflow, calls < maximumEnergyEvaluations else {
            throw VivoChemistryError.resourceLimit("QM/MM center finite-difference step or energy-call budget")
        }
        func value(_ candidate: VivoElectronicSystem) throws -> Double {
            let e = try energy(candidate)
            guard e.isFinite else { throw VivoChemistryError.convergence("nonfinite electronic finite-difference energy") }
            return e
        }
        let baseline = try value(system)
        var nuclear = [VivoVector3D](repeating: .zero, count: system.nuclei.count)
        var charges = [VivoVector3D](repeating: .zero, count: system.pointCharges.count)
        for center in 0..<centers {
            var components = [Double](repeating: 0, count: 3)
            for axis in 0..<3 {
                var plus = system, minus = system
                if center < system.nuclei.count {
                    plus.nuclei[center].positionBohr[axis] += stepBohr
                    minus.nuclei[center].positionBohr[axis] -= stepBohr
                } else {
                    let q = center - system.nuclei.count
                    plus.pointCharges[q].positionBohr[axis] += stepBohr
                    minus.pointCharges[q].positionBohr[axis] -= stepBohr
                }
                components[axis] = -(try value(plus) - value(minus)) / (2 * stepBohr)
            }
            let f = VivoVector3D(components[0], components[1], components[2])
            if center < system.nuclei.count { nuclear[center] = f }
            else { charges[center - system.nuclei.count] = f }
        }
        return try .init(system: system, energyHartree: baseline, nucleusForcesHartreePerBohr: nuclear,
            pointChargeForcesHartreePerBohr: charges, derivativeMethod: "central-difference-all-electronic-centers; step-Bohr=\(stepBohr)")
    }
}

public extension VivoQMMMForceMapper {
    /// Geometry-only projection shared by finite and periodic Hamiltonians. The
    /// returned array contains RAW site forces; a caller combines other raw
    /// contributions before applying the dependent-site graph exactly once.
    static func projectCenters(electronicSystem: VivoElectronicSystem,links: [VivoQMMMLinkAtom],
                               atomToParticle: [UInt32],particlePositionsNM: [VivoVector3D],
                               nucleusForces: [VivoVector3D],pointChargeForces: [VivoVector3D]) throws -> [VivoVector3D] {
        try electronicSystem.validate()
        guard nucleusForces.count == electronicSystem.nuclei.count,pointChargeForces.count == electronicSystem.pointCharges.count,
              nucleusForces.allSatisfy(\.isFinite),pointChargeForces.allSatisfy(\.isFinite),particlePositionsNM.allSatisfy(\.isFinite),
              Set(atomToParticle).count == atomToParticle.count,atomToParticle.allSatisfy({ Int($0) < particlePositionsNM.count }),
              Set(links.map(\.linkNucleusIndex)).count == links.count else { throw VivoChemistryError.invalid("electronic force projection shape or mapping") }
        let linkMap = Dictionary(uniqueKeysWithValues: links.map { ($0.linkNucleusIndex,$0) })
        var output = [VivoVector3D](repeating: .zero,count: particlePositionsNM.count)
        for (index,nucleus) in electronicSystem.nuclei.enumerated() {
            if let atom = nucleus.structureAtomIndex {
                guard Int(atom) < atomToParticle.count else { throw VivoChemistryError.invalid("electronic nucleus has no physical source mapping") }
                let slot = Int(atomToParticle[Int(atom)]),expected = particlePositionsNM[slot]/VivoAtomicUnits.bohrInNM
                guard (expected-VivoVector3D(nucleus.positionBohr.x,nucleus.positionBohr.y,nucleus.positionBohr.z)).norm < 1e-8 else {
                    throw VivoChemistryError.invalid("electronic physical-nucleus geometry does not match source")
                }
                output[slot] = output[slot]+nucleusForces[index]
            } else {
                guard let link = linkMap[index],Int(link.qmAtomIndex) < atomToParticle.count,Int(link.mmAtomIndex) < atomToParticle.count else {
                    throw VivoChemistryError.invalid("unmapped artificial electronic nucleus")
                }
                let q = Int(atomToParticle[Int(link.qmAtomIndex)]),m = Int(atomToParticle[Int(link.mmAtomIndex)])
                let rq = particlePositionsNM[q]/VivoAtomicUnits.bohrInNM,rm = particlePositionsNM[m]/VivoAtomicUnits.bohrInNM
                let expected = rq+(rm-rq)*(link.distanceBohr/(rm-rq).norm)
                guard (expected-VivoVector3D(nucleus.positionBohr.x,nucleus.positionBohr.y,nucleus.positionBohr.z)).norm < 1e-8 else {
                    throw VivoChemistryError.invalid("link geometry differs from its fixed-distance construction")
                }
                let projected = try VivoQMMMLinkForceProjection.project(qmPositionBohr: rq,mmPositionBohr: rm,
                    linkDistanceBohr: link.distanceBohr,linkForceHartreePerBohr: nucleusForces[index])
                output[q] = output[q]+projected.qmForceHartreePerBohr;output[m] = output[m]+projected.mmForceHartreePerBohr
            }
        }
        for (index,charge) in electronicSystem.pointCharges.enumerated() {
            guard let particle = charge.classicalParticleIndex,Int(particle) < output.count else { throw VivoChemistryError.invalid("unmapped electronic charge center") }
            output[Int(particle)] = output[Int(particle)]+pointChargeForces[index]
        }
        guard output.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("electronic force projection overflow") }
        return output
    }
}
