import Foundation

public enum VivoMDElectrostatics: String, Codable, Sendable, CaseIterable { case cutoff, reactionField, pme }
public enum VivoMDEnsemble: String, Codable, Sendable, CaseIterable { case nve, nvt, npt }
public enum VivoMDThermostat: String, Codable, Sendable, CaseIterable { case none, langevinMiddle, velocityRescale }
public enum VivoMDBarostat: String, Codable, Sendable, CaseIterable { case none, monteCarloIsotropic }

public struct VivoMDConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-configuration/v1"
    public var schema: String
    public var timeStepPS: Double
    public var cutoffNM: Double
    public var neighborSkinNM: Double
    public var electrostatics: VivoMDElectrostatics
    public var relativeDielectric: Double
    public var reactionFieldDielectric: Double
    public var ensemble: VivoMDEnsemble
    public var thermostat: VivoMDThermostat
    public var targetTemperatureK: Double?
    public var frictionPerPS: Double?
    public var barostat: VivoMDBarostat
    public var targetPressureBar: Double?
    public var barostatInterval: UInt32
    public var constraintTolerance: Double
    public var maximumConstraintIterations: UInt32
    public var neighborRebuildInterval: UInt32
    /// Optional for v1 backward decoding. nil means enabled.
    public var neighborListEnabled: Bool?
    /// Per-particle list capacity. Overflow rejects the candidate; neighbors are never dropped.
    /// Optional for v1 backward decoding. nil means 512.
    public var maximumNeighborsPerParticle: UInt32?
    public var randomSeed: UInt64

    public init(timeStepPS: Double = 0.002, cutoffNM: Double = 1.0, neighborSkinNM: Double = 0.15,
                electrostatics: VivoMDElectrostatics = .pme, relativeDielectric: Double = 1,
                reactionFieldDielectric: Double = 78.3, ensemble: VivoMDEnsemble = .nvt,
                thermostat: VivoMDThermostat = .langevinMiddle, targetTemperatureK: Double? = 300,
                frictionPerPS: Double? = 1, barostat: VivoMDBarostat = .none,
                targetPressureBar: Double? = nil, barostatInterval: UInt32 = 25,
                constraintTolerance: Double = 1e-6, maximumConstraintIterations: UInt32 = 32,
                neighborRebuildInterval: UInt32 = 10, neighborListEnabled: Bool? = true,
                maximumNeighborsPerParticle: UInt32? = 512,
                randomSeed: UInt64 = 0x4e554d495649564f) {
        schema=Self.schema;self.timeStepPS=timeStepPS;self.cutoffNM=cutoffNM;self.neighborSkinNM=neighborSkinNM
        self.electrostatics=electrostatics;self.relativeDielectric=relativeDielectric;self.reactionFieldDielectric=reactionFieldDielectric
        self.ensemble=ensemble;self.thermostat=thermostat;self.targetTemperatureK=targetTemperatureK;self.frictionPerPS=frictionPerPS
        self.barostat=barostat;self.targetPressureBar=targetPressureBar;self.barostatInterval=barostatInterval
        self.constraintTolerance=constraintTolerance;self.maximumConstraintIterations=maximumConstraintIterations
        self.neighborRebuildInterval=neighborRebuildInterval;self.neighborListEnabled=neighborListEnabled
        self.maximumNeighborsPerParticle=maximumNeighborsPerParticle;self.randomSeed=randomSeed
    }

    public var resolvedNeighborListEnabled: Bool { neighborListEnabled ?? true }
    public var resolvedMaximumNeighborsPerParticle: UInt32 { maximumNeighborsPerParticle ?? 512 }

    public func validate() throws {
        guard schema==Self.schema,timeStepPS.isFinite,timeStepPS>0,cutoffNM.isFinite,cutoffNM>0,
              neighborSkinNM.isFinite,neighborSkinNM>=0,relativeDielectric.isFinite,relativeDielectric>0,
              reactionFieldDielectric.isFinite,reactionFieldDielectric>0,constraintTolerance.isFinite,constraintTolerance>0,
              maximumConstraintIterations>0,neighborRebuildInterval>0,barostatInterval>0,
              resolvedMaximumNeighborsPerParticle>0 else {
            throw VivoArtifactValidationError.invalid("MD configuration contains invalid numerical settings")
        }
        if resolvedNeighborListEnabled, neighborSkinNM <= 0 {
            throw VivoArtifactValidationError.invalid("Verlet neighbor-list execution requires a positive skin")
        }
        if thermostat == .none {
            guard ensemble == .nve, targetTemperatureK == nil else {
                throw VivoArtifactValidationError.invalid("thermostat-free execution is NVE and must not declare a target temperature")
            }
        } else {
            guard ensemble != .nve,let t=targetTemperatureK,t.isFinite,t>0 else {
                throw VivoArtifactValidationError.invalid("thermostatted execution requires NVT/NPT and positive temperature")
            }
        }
        if thermostat == .langevinMiddle {
            guard let f=frictionPerPS,f.isFinite,f>=0 else {
                throw VivoArtifactValidationError.invalid("Langevin-middle requires finite nonnegative friction")
            }
        }
        if ensemble == .npt || barostat != .none {
            guard ensemble == .npt,barostat != .none,let p=targetPressureBar,p.isFinite,p>0 else {
                throw VivoArtifactValidationError.invalid("NPT requires a barostat and positive pressure")
            }
        }
        if electrostatics == .reactionField,reactionFieldDielectric<relativeDielectric {
            throw VivoArtifactValidationError.invalid("reaction-field dielectric cannot be smaller than internal dielectric")
        }
    }
    public func fingerprint() throws -> VivoFingerprint { try validate();return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self)) }
}

public struct VivoMDStateSnapshot: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/md-state/v1"
    public var schema:String;public var systemFingerprint:VivoFingerprint;public var configurationFingerprint:VivoFingerprint
    public var stepIndex:UInt64;public var timePS:Double;public var positionsNM:[VivoVector3D];public var velocitiesNMPerPS:[VivoVector3D]
    public var periodicCell:VivoPeriodicCell?;public var potentialEnergyKJPerMol:Double?;public var kineticEnergyKJPerMol:Double?;public var temperatureK:Double?
    public init(systemFingerprint:VivoFingerprint,configurationFingerprint:VivoFingerprint,stepIndex:UInt64,timePS:Double,
                positionsNM:[VivoVector3D],velocitiesNMPerPS:[VivoVector3D],periodicCell:VivoPeriodicCell?,potentialEnergyKJPerMol:Double?=nil,
                kineticEnergyKJPerMol:Double?=nil,temperatureK:Double?=nil){schema=Self.schema;self.systemFingerprint=systemFingerprint;self.configurationFingerprint=configurationFingerprint;self.stepIndex=stepIndex;self.timePS=timePS;self.positionsNM=positionsNM;self.velocitiesNMPerPS=velocitiesNMPerPS;self.periodicCell=periodicCell;self.potentialEnergyKJPerMol=potentialEnergyKJPerMol;self.kineticEnergyKJPerMol=kineticEnergyKJPerMol;self.temperatureK=temperatureK}
    public func validate(particleCount:Int)throws{guard schema==Self.schema,positionsNM.count==particleCount,velocitiesNMPerPS.count==particleCount,positionsNM.allSatisfy(\.isFinite),velocitiesNMPerPS.allSatisfy(\.isFinite),timePS.isFinite,timePS>=0,periodicCell?.isValid != false,potentialEnergyKJPerMol?.isFinite != false,kineticEnergyKJPerMol?.isFinite != false,temperatureK?.isFinite != false else{throw VivoArtifactValidationError.invalid("MD state snapshot is malformed or non-finite")}}
}

public struct VivoMDCapabilityReport: Codable, Sendable, Equatable {
    public var executable:Bool;public var particleCount:UInt32;public var massiveParticleCount:UInt32;public var virtualSiteCount:UInt32
    public var bondCount:UInt32;public var angleCount:UInt32;public var torsionCount:UInt32;public var constraintCount:UInt32
    public var explicitPairTypeCount:UInt32;public var blockers:[String];public var notes:[String]
}

public enum VivoMDCapabilityAnalyzer {
    public static func analyze(system:VivoClassicalSystem,initialState:VivoClassicalInitialState,configuration:VivoMDConfiguration)throws->VivoMDCapabilityReport{
        try configuration.validate();try VivoClassicalSystemValidator.validate(system);try initialState.validate(particleCount:system.particles.count)
        guard initialState.systemFingerprint == system.fingerprint() else{throw VivoArtifactValidationError.incompatible("MD initial state does not identify the supplied classical system")}
        let virtualSites=system.particles.filter{$0.role == .virtualSite};var blockers:[String]=[],notes:[String]=[]
        if !virtualSites.isEmpty{blockers.append("virtual-site construction rules are not yet represented")}
        if system.particles.contains(where:{$0.role == .drude}){blockers.append("Drude polarization is not implemented")}
        if configuration.electrostatics == .pme,initialState.periodicCell==nil{blockers.append("PME requires a periodic cell")}
        if configuration.ensemble == .npt,initialState.periodicCell==nil{blockers.append("NPT requires a periodic cell")}
        let pairs=system.nonbondedTypePairs ?? [];if system.mixingRule == .explicitPairTable,pairs.isEmpty{blockers.append("explicit pair table is absent")}
        if configuration.resolvedNeighborListEnabled,UInt64(configuration.resolvedMaximumNeighborsPerParticle)>UInt64(max(system.particles.count-1,0)){notes.append("neighbor capacity exceeds maximum possible neighbors; allocator may clamp it")}
        if configuration.electrostatics == .pme{notes.append("PME requires reciprocal-space kernels in the execution backend")}
        if configuration.barostat != .none{notes.append("barostat requires transactional cell scaling")}
        return .init(executable:blockers.isEmpty,particleCount:UInt32(system.particles.count),massiveParticleCount:UInt32(system.particles.filter{$0.massDa>0}.count),virtualSiteCount:UInt32(virtualSites.count),bondCount:UInt32(system.bonds.count),angleCount:UInt32(system.angles.count),torsionCount:UInt32(system.torsions.count),constraintCount:UInt32(system.constraints.count),explicitPairTypeCount:UInt32(pairs.count),blockers:blockers,notes:notes)
    }
}
