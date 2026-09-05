import Foundation

public struct VivoMDObservables: Codable, Sendable, Equatable {
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let stepIndex: UInt64
    public let timePS: Double
    public let potentialEnergyKJPerMol: Double
    public let kineticEnergyKJPerMol: Double
    public let totalEnergyKJPerMol: Double
    public let temperatureK: Double?
    public let degreesOfFreedom: UInt64

    public init(systemFingerprint: VivoFingerprint, configurationFingerprint: VivoFingerprint,
                stepIndex: UInt64, timePS: Double,
                potentialEnergyKJPerMol: Double, kineticEnergyKJPerMol: Double,
                temperatureK: Double?, degreesOfFreedom: UInt64) {
        self.systemFingerprint=systemFingerprint; self.configurationFingerprint=configurationFingerprint
        self.stepIndex=stepIndex; self.timePS=timePS
        self.potentialEnergyKJPerMol=potentialEnergyKJPerMol; self.kineticEnergyKJPerMol=kineticEnergyKJPerMol
        totalEnergyKJPerMol=potentialEnergyKJPerMol+kineticEnergyKJPerMol
        self.temperatureK=temperatureK; self.degreesOfFreedom=degreesOfFreedom
    }
}

/// Accepted-boundary restart. Random evolution is reproduced from acceptedStep +
/// configuration seed, so no mutable RNG state exists outside this payload.
public struct VivoMDCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-checkpoint/v1"
    public var schema: String
    public var systemFingerprint: VivoFingerprint
    public var configurationFingerprint: VivoFingerprint
    public var acceptedStep: UInt64
    public var timePS: Double
    public var positionsNM: [VivoVector3D]
    public var velocitiesNMPerPS: [VivoVector3D]
    public var periodicCell: VivoPeriodicCell?

    public init(systemFingerprint: VivoFingerprint, configurationFingerprint: VivoFingerprint,
                acceptedStep: UInt64, timePS: Double,
                positionsNM: [VivoVector3D], velocitiesNMPerPS: [VivoVector3D],
                periodicCell: VivoPeriodicCell?) {
        schema=Self.schema; self.systemFingerprint=systemFingerprint
        self.configurationFingerprint=configurationFingerprint; self.acceptedStep=acceptedStep
        self.timePS=timePS; self.positionsNM=positionsNM; self.velocitiesNMPerPS=velocitiesNMPerPS
        self.periodicCell=periodicCell
    }

    public func validate(particleCount: Int) throws {
        guard schema==Self.schema, positionsNM.count==particleCount,
              velocitiesNMPerPS.count==particleCount, positionsNM.allSatisfy(\.isFinite),
              velocitiesNMPerPS.allSatisfy(\.isFinite), timePS.isFinite, timePS>=0,
              periodicCell?.isValid != false else {
            throw VivoArtifactValidationError.invalid("MD checkpoint shape, clock, cell, or numerical values are invalid")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate(particleCount: positionsNM.count)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoMDTrajectorySample: Codable, Sendable, Equatable {
    public let stepIndex: UInt64
    public let timePS: Double
    public let positionsNM: [VivoVector3D]
    public let periodicCell: VivoPeriodicCell?
    public let observables: VivoMDObservables?
}

public struct VivoMDRunReport: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-run-report/v1"
    public var schema: String
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let requestedSteps: UInt64
    public let committedSteps: UInt64
    public let startStep: UInt64
    public let endStep: UInt64
    public let startTimePS: Double
    public let endTimePS: Double
    public let rejected: VivoMDStepCertificate?
    public let samples: [VivoMDTrajectorySample]
    public let finalCheckpoint: VivoMDCheckpoint

    public init(systemFingerprint: VivoFingerprint, configurationFingerprint: VivoFingerprint,
                deviceName: String, deviceRegistryID: UInt64,
                requestedSteps: UInt64, committedSteps: UInt64,
                startStep: UInt64, endStep: UInt64,
                startTimePS: Double, endTimePS: Double,
                rejected: VivoMDStepCertificate?, samples: [VivoMDTrajectorySample],
                finalCheckpoint: VivoMDCheckpoint) {
        schema=Self.schema; self.systemFingerprint=systemFingerprint
        self.configurationFingerprint=configurationFingerprint; self.deviceName=deviceName
        self.deviceRegistryID=deviceRegistryID; self.requestedSteps=requestedSteps
        self.committedSteps=committedSteps; self.startStep=startStep; self.endStep=endStep
        self.startTimePS=startTimePS; self.endTimePS=endTimePS
        self.rejected=rejected; self.samples=samples; self.finalCheckpoint=finalCheckpoint
    }
}
