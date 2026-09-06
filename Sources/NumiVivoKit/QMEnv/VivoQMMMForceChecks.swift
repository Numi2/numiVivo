import Foundation

public struct VivoQMMMForceCheckConfiguration: Codable, Sendable, Equatable {
    public var stepBohr: Double
    public var absoluteForceTolerance: Double
    public var relativeForceTolerance: Double
    public var netForceTolerance: Double
    public var netTorqueTolerance: Double
    public var energyToleranceHartree: Double
    public var maximumEnergyEvaluations: Int
    /// nil checks every physical atom. Explicit subsets are recorded as partial.
    public var atomIndices: [UInt32]?
    public init(stepBohr: Double = 1e-4, absoluteForceTolerance: Double = 1e-5,
                relativeForceTolerance: Double = 1e-4, netForceTolerance: Double = 1e-7,
                netTorqueTolerance: Double = 1e-6, energyToleranceHartree: Double = 1e-9,
                maximumEnergyEvaluations: Int = 4096, atomIndices: [UInt32]? = nil) {
        self.stepBohr = stepBohr; self.absoluteForceTolerance = absoluteForceTolerance
        self.relativeForceTolerance = relativeForceTolerance; self.netForceTolerance = netForceTolerance
        self.netTorqueTolerance = netTorqueTolerance; self.energyToleranceHartree = energyToleranceHartree
        self.maximumEnergyEvaluations = maximumEnergyEvaluations; self.atomIndices = atomIndices
    }
    public func validate() throws {
        guard stepBohr.isFinite, stepBohr >= 1e-8, stepBohr <= 0.01,
              [absoluteForceTolerance, netForceTolerance, netTorqueTolerance, energyToleranceHartree].allSatisfy({ $0.isFinite && $0 > 0 }),
              relativeForceTolerance.isFinite, relativeForceTolerance >= 0, maximumEnergyEvaluations > 0 else {
            throw VivoChemistryError.invalid("QM/MM force-check configuration")
        }
    }
}

public struct VivoQMMMForceCheckComponent: Codable, Sendable, Equatable {
    public let atomIndex: UInt32
    public let particleIndex: UInt32
    public let axis: Int
    public let assembledForceHartreePerBohr: Double
    public let finiteDifferenceForceHartreePerBohr: Double
    public let absoluteError: Double
    public let stepHalvingChange: Double
    public let tolerance: Double
    public let passed: Bool
}
public struct VivoQMMMForceCheckReport: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-force-check/v1"
    public let schema: String
    public let passed: Bool
    public let electronicDerivativeMethod: String
    public let configuration: VivoQMMMForceCheckConfiguration
    public let sourceFrameFingerprint: VivoFingerprint
    public let finiteClusterFingerprint: VivoFingerprint
    public let preparedSnapshotFingerprint: VivoFingerprint
    public let checkedAtomIndices: [UInt32]
    public let completePhysicalAtomCheck: Bool
    public let energyEvaluations: Int
    public let electronicBaselineEnergyDifferenceHartree: Double
    public let maximumAbsoluteError: Double
    public let rmsError: Double
    public let maximumStepHalvingChange: Double
    public let netForceHartreePerBohr: VivoVector3D
    public let netTorqueHartree: VivoVector3D
    public let components: [VivoQMMMForceCheckComponent]
    public let interpretation: String
}

public enum VivoQMMMForceChecks {
    /// Compare assembled physical forces to independent TOTAL-energy differences.
    /// Keep membership and lattice images fixed. At each perturbation rebuild all
    /// virtual sites, link nuclei and embedding charges from the moved physical
    /// coordinates. Moving a link independently or freezing MM charges is invalid.
    public static func check(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                             frame: VivoMDQMMMPreparedFrame, electronic: VivoQMMMElectronicForceResult,
                             configuration cfg: VivoQMMMForceCheckConfiguration = .init(),
                             budget: VivoChemistryBudget = .init(),
                             electronicEnergy: (VivoElectronicSystem) throws -> Double) throws -> VivoQMMMForceCheckReport {
        try cfg.validate(); try budget.validate()
        let result = try VivoQMMMForceMapper.assemble(document: document, system: system, frame: frame,
                                                     electronic: electronic, budget: budget)
        let topology = try VivoQMMMTopology(document: document, system: system)
        let atoms = cfg.atomIndices ?? Array(0..<UInt32(topology.atomToParticle.count))
        guard !atoms.isEmpty, Set(atoms).count == atoms.count, atoms.allSatisfy({ Int($0) < topology.atomToParticle.count }) else {
            throw VivoChemistryError.invalid("QM/MM force-check atom selection")
        }
        let (required, overflow) = atoms.count.multipliedReportingOverflow(by: 12)
        guard !overflow, required < cfg.maximumEnergyEvaluations else {
            throw VivoChemistryError.resourceLimit("QM/MM force check exceeds energy-call budget")
        }
        var calls = 0
        func value(_ candidate: VivoElectronicSystem) throws -> Double {
            guard calls < cfg.maximumEnergyEvaluations else { throw VivoChemistryError.resourceLimit("QM/MM force-check energy-call limit") }
            calls += 1
            let value = try electronicEnergy(candidate)
            guard value.isFinite else { throw VivoChemistryError.convergence("nonfinite QM/MM force-check electronic energy") }
            return value
        }
        let baselineDifference = abs(try value(frame.region.electronicSystem) - electronic.energyHartree)
        func totalEnergy(_ positions: [VivoVector3D]) throws -> Double {
            var rebuilt = positions
            try topology.rebuildVirtualSites(&rebuilt)
            let region = try VivoQMMMCompiler.prepare(document: document, system: system,
                particlePositionsNM: rebuilt, request: frame.request)
            guard region.regionPolicyFingerprint == frame.region.regionPolicyFingerprint,
                  region.links == frame.region.links, region.excludedQMParticles == frame.region.excludedQMParticles,
                  region.z1ZeroedMMParticles == frame.region.z1ZeroedMMParticles else {
                throw VivoChemistryError.invalid("QM/MM membership changed during force differentiation")
            }
            let e = try value(region.electronicSystem)
            let lj = try VivoQMMMCompiler.lennardJones(document: document, system: system,
                particlePositionsNM: rebuilt, request: frame.request, budget: budget)
            guard (e + lj.energyHartree).isFinite else { throw VivoChemistryError.convergence("nonfinite total QM/MM check energy") }
            return e + lj.energyHartree
        }
        func component(_ p: VivoVector3D, _ axis: Int) -> Double {
            axis == 0 ? p.x : (axis == 1 ? p.y : p.z)
        }
        func displaced(_ p: VivoVector3D, _ axis: Int, _ d: Double) -> VivoVector3D {
            switch axis {
            case 0: return .init(p.x + d, p.y, p.z)
            case 1: return .init(p.x, p.y + d, p.z)
            default: return .init(p.x, p.y, p.z + d)
            }
        }
        var checks: [VivoQMMMForceCheckComponent] = [], sumSquares = 0.0, maxError = 0.0, maxChange = 0.0
        for atom in atoms.sorted() {
            let particle = topology.atomToParticle[Int(atom)], slot = Int(particle)
            for axis in 0..<3 {
                func derivative(_ step: Double) throws -> Double {
                    var plus = frame.cluster.particlePositionsNM, minus = plus
                    plus[slot] = displaced(plus[slot], axis, step * VivoAtomicUnits.bohrInNM)
                    minus[slot] = displaced(minus[slot], axis, -step * VivoAtomicUnits.bohrInNM)
                    return -(try totalEnergy(plus) - totalEnergy(minus)) / (2 * step)
                }
                let coarse = try derivative(cfg.stepBohr), fine = try derivative(cfg.stepBohr / 2)
                let assembled = component(result.particleForcesHartreePerBohr[slot], axis)
                let error = abs(assembled - fine), change = abs(fine - coarse)
                let tolerance = cfg.absoluteForceTolerance + cfg.relativeForceTolerance * max(abs(assembled), abs(fine))
                guard [coarse, fine, error, change, tolerance].allSatisfy(\.isFinite) else {
                    throw VivoChemistryError.convergence("nonfinite force-check difference")
                }
                maxError = max(maxError, error); maxChange = max(maxChange, change); sumSquares += error * error
                checks.append(.init(atomIndex: atom, particleIndex: particle, axis: axis,
                    assembledForceHartreePerBohr: assembled, finiteDifferenceForceHartreePerBohr: fine,
                    absoluteError: error, stepHalvingChange: change, tolerance: tolerance,
                    passed: error <= tolerance && change <= tolerance))
            }
        }
        guard sumSquares.isFinite else { throw VivoChemistryError.convergence("force-check error norm overflow") }
        let passed = checks.allSatisfy(\.passed) && baselineDifference <= cfg.energyToleranceHartree
            && result.netForceHartreePerBohr.norm <= cfg.netForceTolerance && result.netTorqueHartree.norm <= cfg.netTorqueTolerance
        return .init(schema: VivoQMMMForceCheckReport.schema, passed: passed,
            electronicDerivativeMethod: electronic.derivativeMethod, configuration: cfg,
            sourceFrameFingerprint: result.sourceFrameFingerprint, finiteClusterFingerprint: result.finiteClusterFingerprint,
            preparedSnapshotFingerprint: result.preparedSnapshotFingerprint, checkedAtomIndices: atoms.sorted(),
            completePhysicalAtomCheck: atoms.count == topology.atomToParticle.count, energyEvaluations: calls,
            electronicBaselineEnergyDifferenceHartree: baselineDifference, maximumAbsoluteError: maxError,
            rmsError: sqrt(sumSquares / Double(checks.count)), maximumStepHalvingChange: maxChange,
            netForceHartreePerBohr: result.netForceHartreePerBohr, netTorqueHartree: result.netTorqueHartree,
            components: checks,
            interpretation: "fixed-membership, fixed-image finite-cluster energy/force check with link and virtual-site reconstruction; two step sizes; all-center translation/rotation diagnostics; not periodic Ewald QM/MM, adaptive-region dynamics or a complete classical-force replacement")
    }
}
