import Foundation

public struct VivoContextApplicationReceipt: Codable, Sendable, Equatable {
    public let contextFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let environmentCount: UInt32
    public let laneCount: UInt32
    public let parameterValueCount: Int
    public let transportRecordCount: Int
    public let pendingCouplingCount: Int
    public let diagnostics: [VivoContextDiagnostic]

    public init(
        contextFingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint,
        environmentCount: UInt32,
        laneCount: UInt32,
        parameterValueCount: Int,
        transportRecordCount: Int,
        pendingCouplingCount: Int,
        diagnostics: [VivoContextDiagnostic]
    ) {
        self.contextFingerprint = contextFingerprint
        self.programFingerprint = programFingerprint
        self.environmentCount = environmentCount
        self.laneCount = laneCount
        self.parameterValueCount = parameterValueCount
        self.transportRecordCount = transportRecordCount
        self.pendingCouplingCount = pendingCouplingCount
        self.diagnostics = diagnostics
    }
}
