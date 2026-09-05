import Foundation

public enum VivoMDElectrostatics:String,Codable,Sendable,CaseIterable{case cutoff,reactionField,pme}
public enum VivoMDEnsemble:String,Codable,Sendable,CaseIterable{case nve,nvt,npt}
public enum VivoMDThermostat:String,Codable,Sendable,CaseIterable{case none,langevinMiddle,velocityRescale}
public enum VivoMDBarostat:String,Codable,Sendable,CaseIterable{case none,monteCarloIsotropic}

public struct VivoMDConfiguration:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/md-configuration/v1"
    public var schema:String
    public var timeStepPS:Double
    public var cutoffNM:Double
    public var neighborSkinNM:Double
    public var electrostatics:VivoMDElectrostatics
    public var relativeDielectric:Double
    public var reactionFieldDielectric:Double
    public var pmeTolerance:Double?
    public var pmeGridSpacingNM:Double?
    public var ensemble:VivoMDEnsemble
    public var thermostat:VivoMDThermostat
    public var targetTemperatureK:Double?
    public var frictionPerPS:Double?
    public var barostat:VivoMDBarostat
    public var targetPressureBar:Double?
    public var barostatInterval:UInt32
    public var barostatMaximumLogVolumeStep:Double?
    public var constraintTolerance:Double
    public var maximumConstraintIterations:UInt32
    public var neighborRebuildInterval:UInt32
    public var neighborListEnabled:Bool?
    public var maximumNeighborsPerParticle:UInt32?
    public var randomSeed:UInt64

    public init(timeStepPS:Double=0.002,cutoffNM:Double=1,neighborSkinNM:Double=0.15,
                electrostatics:VivoMDElectrostatics = .pme,relativeDielectric:Double=1,
                reactionFieldDielectric:Double=78.3,pmeTolerance:Double?=1e-5,
                pmeGridSpacingNM:Double?=0.12,ensemble:VivoMDEnsemble = .nvt,
                thermostat:VivoMDThermostat = .langevinMiddle,targetTemperatureK:Double?=300,
                frictionPerPS:Double?=1,barostat:VivoMDBarostat = .none,
                targetPressureBar:Double?=nil,barostatInterval:UInt32=25,
                barostatMaximumLogVolumeStep:Double?=0.01,
                constraintTolerance:Double=1e-6,maximumConstraintIterations:UInt32=32,
                neighborRebuildInterval:UInt32=10,neighborListEnabled:Bool?=true,
                maximumNeighborsPerParticle:UInt32?=512,randomSeed:UInt64=0x4e554d495649564f) {
        schema=Self.schema;self.timeStepPS=timeStepPS;self.cutoffNM=cutoffNM;self.neighborSkinNM=neighborSkinNM
        self.electrostatics=electrostatics;self.relativeDielectric=relativeDielectric;self.reactionFieldDielectric=reactionFieldDielectric
        self.pmeTolerance=pmeTolerance;self.pmeGridSpacingNM=pmeGridSpacingNM;self.ensemble=ensemble;self.thermostat=thermostat
        self.targetTemperatureK=targetTemperatureK;self.frictionPerPS=frictionPerPS;self.barostat=barostat;self.targetPressureBar=targetPressureBar
        self.barostatInterval=barostatInterval;self.barostatMaximumLogVolumeStep=barostatMaximumLogVolumeStep
        self.constraintTolerance=constraintTolerance;self.maximumConstraintIterations=maximumConstraintIterations
        self.neighborRebuildInterval=neighborRebuildInterval;self.neighborListEnabled=neighborListEnabled
        self.maximumNeighborsPerParticle=maximumNeighborsPerParticle;self.randomSeed=randomSeed
    }
    public var resolvedNeighborListEnabled:Bool{neighborListEnabled ?? true}
    public var resolvedMaximumNeighborsPerParticle:UInt32{maximumNeighborsPerParticle ?? 512}
    public var resolvedPMETolerance:Double{pmeTolerance ?? 1e-5}
    public var resolvedPMEGridSpacingNM:Double{pmeGridSpacingNM ?? 0.12}
    public var resolvedBarostatMaximumLogVolumeStep:Double{barostatMaximumLogVolumeStep ?? 0.01}

    public func validate()throws {
        func positiveFP32(_ x:Double)->Bool{x.isFinite && x>0 && Float(x).isFinite && Float(x)>0}
        guard schema==Self.schema,positiveFP32(timeStepPS),positiveFP32(cutoffNM),
              neighborSkinNM.isFinite,neighborSkinNM>=0,Float(neighborSkinNM).isFinite,
              positiveFP32(relativeDielectric),positiveFP32(reactionFieldDielectric),positiveFP32(constraintTolerance),
              maximumConstraintIterations>0,maximumConstraintIterations<=16_384,
              neighborRebuildInterval>0,barostatInterval>0,resolvedMaximumNeighborsPerParticle>0 else {
            throw VivoArtifactValidationError.invalid("MD configuration contains invalid/unsupported numerical settings")
        }
        if electrostatics == .pme {
            guard resolvedPMETolerance.isFinite,resolvedPMETolerance>0,resolvedPMETolerance<0.1,
                  positiveFP32(resolvedPMEGridSpacingNM) else {throw VivoArtifactValidationError.invalid("PME planning parameters are invalid")}
        }
        if resolvedNeighborListEnabled, !positiveFP32(neighborSkinNM) {throw VivoArtifactValidationError.invalid("neighbor lists require a positive representable skin")}
        if thermostat == .none {
            guard ensemble == .nve,targetTemperatureK==nil else {throw VivoArtifactValidationError.invalid("thermostat-free execution must be NVE without a target temperature")}
        } else {
            guard ensemble != .nve,let t=targetTemperatureK,positiveFP32(t) else {throw VivoArtifactValidationError.invalid("thermostatted execution requires a positive temperature")}
        }
        if thermostat == .langevinMiddle {
            guard let f=frictionPerPS,f.isFinite,f>=0,(f*timeStepPS).isFinite else {throw VivoArtifactValidationError.invalid("invalid Langevin friction")}
        }
        if ensemble == .npt || barostat != .none {
            guard ensemble == .npt,barostat == .monteCarloIsotropic,let p=targetPressureBar,p.isFinite,p>0,
                  resolvedBarostatMaximumLogVolumeStep.isFinite,resolvedBarostatMaximumLogVolumeStep>0,
                  resolvedBarostatMaximumLogVolumeStep<=0.25 else {throw VivoArtifactValidationError.invalid("invalid NPT pressure/proposal settings")}
        }
        if electrostatics == .reactionField,reactionFieldDielectric<relativeDielectric {
            throw VivoArtifactValidationError.invalid("reaction-field dielectric is below internal dielectric")
        }
    }
    public func fingerprint()throws->VivoFingerprint {
        try validate();return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoMDStateSnapshot:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/md-state/v1"
    public var schema:String
    public var systemFingerprint:VivoFingerprint
    public var configurationFingerprint:VivoFingerprint
    public var stepIndex:UInt64
    public var timePS:Double
    public var positionsNM:[VivoVector3D]
    public var velocitiesNMPerPS:[VivoVector3D]
    public var periodicCell:VivoPeriodicCell?
    public var potentialEnergyKJPerMol:Double?
    public var kineticEnergyKJPerMol:Double?
    public var temperatureK:Double?
    public init(systemFingerprint:VivoFingerprint,configurationFingerprint:VivoFingerprint,stepIndex:UInt64,timePS:Double,
                positionsNM:[VivoVector3D],velocitiesNMPerPS:[VivoVector3D],periodicCell:VivoPeriodicCell?,
                potentialEnergyKJPerMol:Double?=nil,kineticEnergyKJPerMol:Double?=nil,temperatureK:Double?=nil) {
        schema=Self.schema;self.systemFingerprint=systemFingerprint;self.configurationFingerprint=configurationFingerprint
        self.stepIndex=stepIndex;self.timePS=timePS;self.positionsNM=positionsNM;self.velocitiesNMPerPS=velocitiesNMPerPS
        self.periodicCell=periodicCell;self.potentialEnergyKJPerMol=potentialEnergyKJPerMol
        self.kineticEnergyKJPerMol=kineticEnergyKJPerMol;self.temperatureK=temperatureK
    }
    public func validate(particleCount:Int)throws {
        guard schema==Self.schema,particleCount>0,positionsNM.count==particleCount,velocitiesNMPerPS.count==particleCount,
              positionsNM.allSatisfy(\.isFinite),velocitiesNMPerPS.allSatisfy(\.isFinite),timePS.isFinite,timePS>=0,
              periodicCell?.isValid != false,potentialEnergyKJPerMol?.isFinite != false,
              kineticEnergyKJPerMol?.isFinite != false,temperatureK?.isFinite != false else {
            throw VivoArtifactValidationError.invalid("MD snapshot shape, clock or numerical values are invalid")
        }
    }
}

public struct VivoMDCapabilityReport:Codable,Sendable,Equatable {
    public var executable:Bool
    public var particleCount:UInt32
    public var massiveParticleCount:UInt32
    public var virtualSiteCount:UInt32
    public var bondCount:UInt32
    public var angleCount:UInt32
    public var torsionCount:UInt32
    public var constraintCount:UInt32
    public var explicitPairTypeCount:UInt32
    public var blockers:[String]
    public var notes:[String]
    public var numericalContract:String = VivoMDExecutionIdentity.current
    public var qualification:String = "source-implemented; Apple compilation and numerical qualification pending"
}
public enum VivoMDCapabilityAnalyzer {
    public static func analyze(system:VivoClassicalSystem,initialState:VivoClassicalInitialState,
                               configuration:VivoMDConfiguration)throws->VivoMDCapabilityReport {
        try configuration.validate();try VivoClassicalSystemValidator.validate(system)
        try initialState.validate(particleCount:system.particles.count)
        guard initialState.systemFingerprint == (try system.fingerprint()) else {
            throw VivoArtifactValidationError.incompatible("MD initial state does not identify the supplied system")
        }
        let virtual=system.particles.filter{$0.role == .virtualSite}
        let resolved=Set((system.linearVirtualSites ?? []).map(\.siteParticle))
        let unresolved=virtual.map(\.index).filter{!resolved.contains($0)}
        var blockers=VivoMDExecutionPreflight.blockers(system:system,initial:initialState,configuration:configuration)
        var notes:[String]=[]
        if !unresolved.isEmpty {blockers.append("unresolved virtual-site rules: \(unresolved)")}
        if system.particles.contains(where:{$0.role == .drude}) {blockers.append("Drude polarization is not implemented")}
        if !virtual.isEmpty {notes.append("virtual forces redistribute to physical parents; virtual slots are not minimization degrees of freedom")}
        if configuration.electrostatics == .pme {notes.append("cubic mesh/FFT PME; the input tolerance is a planning heuristic, not a demonstrated force-error bound")}
        if configuration.ensemble == .npt {notes.append("molecular-center log-volume proposals; operational failures abort rather than count as Metropolis rejections")}
        if configuration.resolvedNeighborListEnabled {notes.append("neighbor lists are rebuilt at both force positions; rebuildInterval is retained for compatibility, not yet an amortization guarantee")}
        notes.append("reported thermal degrees of freedom assume independent distance constraints")
        notes.append("executable means no known source-contract blocker; it does not establish Metal availability or numerical correctness")
        return .init(executable:blockers.isEmpty,particleCount:UInt32(system.particles.count),
            massiveParticleCount:UInt32(system.particles.filter{$0.massDa>0}.count),virtualSiteCount:UInt32(virtual.count),
            bondCount:UInt32(system.bonds.count),angleCount:UInt32(system.angles.count),torsionCount:UInt32(system.torsions.count),
            constraintCount:UInt32(system.constraints.count),explicitPairTypeCount:UInt32((system.nonbondedTypePairs ?? []).count),
            blockers:Array(Set(blockers)).sorted(),notes:notes)
    }
}
