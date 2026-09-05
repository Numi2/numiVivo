import Foundation

public struct VivoMeasurementSummary: Codable, Sendable, Equatable {
    public let measurementIdentifier: String
    public let replicateIndex: UInt32
    public let sampleCount: UInt64
    public let firstTimeSeconds: Double
    public let lastTimeSeconds: Double
    public let mean: [Double]
    public let minimum: [Double]
    public let maximum: [Double]
    public let last: [Double]
    public let unit: String
    public init(measurementIdentifier: String, replicateIndex: UInt32, sampleCount: UInt64,
                firstTimeSeconds: Double, lastTimeSeconds: Double, mean: [Double],
                minimum: [Double], maximum: [Double], last: [Double], unit: String) {
        self.measurementIdentifier=measurementIdentifier;self.replicateIndex=replicateIndex
        self.sampleCount=sampleCount;self.firstTimeSeconds=firstTimeSeconds;self.lastTimeSeconds=lastTimeSeconds
        self.mean=mean;self.minimum=minimum;self.maximum=maximum;self.last=last;self.unit=unit
    }
    public func validate() throws {
        guard !measurementIdentifier.isEmpty,!unit.isEmpty,sampleCount>0,
              firstTimeSeconds.isFinite,lastTimeSeconds.isFinite,firstTimeSeconds>=0,lastTimeSeconds>=firstTimeSeconds,
              !mean.isEmpty,minimum.count==mean.count,maximum.count==mean.count,last.count==mean.count,
              (mean+minimum+maximum+last).allSatisfy(\.isFinite) else {
            throw VivoArtifactValidationError.invalid("measurement summary shape, values or time")
        }
        for i in mean.indices {
            let tolerance=1e-12*max(1,abs(minimum[i]),abs(maximum[i]))
            guard minimum[i]<=maximum[i],mean[i]>=minimum[i]-tolerance,mean[i]<=maximum[i]+tolerance,
                  last[i]>=minimum[i]-tolerance,last[i]<=maximum[i]+tolerance else {
                throw VivoArtifactValidationError.invalid("measurement summary bounds")
            }
        }
    }
}
public struct VivoRecordedEvent: Codable, Sendable, Equatable {
    public let replicateIndex: UInt32
    public let event: VivoEvent
    public init(replicateIndex: UInt32,event: VivoEvent) {self.replicateIndex=replicateIndex;self.event=event}
}
/// The legacy checkpoint pack omits spatial layout, velocity and volume fields.
/// Experiment exports now retain the authoritative configuration-bound payload
/// rather than down-converting it to an incomplete state snapshot.
public struct VivoExperimentCheckpoint: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/experiment-resume-checkpoint/v1"
    public let schema: String
    public let state: VivoMolecularResumeCheckpoint
    public let replicateIndex: UInt32
    public let experimentFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint
    public let couplingFingerprint: VivoFingerprint?
    public let ledgerHead: VivoFingerprint
    public init(state: VivoMolecularResumeCheckpoint, replicateIndex: UInt32,experimentFingerprint: VivoFingerprint,
                hostContextFingerprint: VivoFingerprint,couplingFingerprint: VivoFingerprint?,ledgerHead: VivoFingerprint) throws {
        self.schema=Self.schema;self.state=state;self.replicateIndex=replicateIndex
        self.experimentFingerprint=experimentFingerprint;self.hostContextFingerprint=hostContextFingerprint
        self.couplingFingerprint=couplingFingerprint;self.ledgerHead=ledgerHead
        try validate()
    }
    public func validate() throws {
        guard schema==Self.schema else {throw VivoArtifactValidationError.invalid("experiment checkpoint schema")}
        try state.validate()
    }
    public func fingerprint() throws -> VivoFingerprint {
        try validate();return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}
