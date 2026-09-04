import Foundation

public enum VivoNumiLabChannelID {
    // NumanX -> NumiVivo
    public static let tissuePressure = "numanx.tissue.pressure"
    public static let tissueTemperature = "numanx.tissue.temperature"
    public static let tissueOxygen = "numanx.tissue.oxygen"
    public static let tissuePerfusion = "numanx.tissue.perfusion"
    public static let tissueShearRate = "numanx.tissue.shear-rate"
    public static let tissueVolumetricStrain = "numanx.tissue.volumetric-strain"

    // NumiVivo -> NumanX
    public static let osmoticPressure = "numivivo.tissue.osmotic-pressure"
    public static let activeStress = "numivivo.tissue.active-stress"
    public static let fluidSource = "numivivo.tissue.fluid-source"
    public static let permeabilityModifier = "numivivo.tissue.permeability-modifier"
    public static let viscosityModifier = "numivivo.tissue.viscosity-modifier"

    // NumiTissue -> NumiVivo
    public static let cellVolumeFraction = "numitissue.cell.volume-fraction"
    public static let membraneAreaDensity = "numitissue.cell.membrane-area-density"
    public static let extracellularVolumeFraction = "numitissue.ecm.volume-fraction"
    public static let phenotypeDistribution = "numitissue.cell.phenotype-distribution"

    // NumiVivo -> NumiTissue
    public static let proliferationDrive = "numivivo.cell.proliferation-drive"
    public static let differentiationDrive = "numivivo.cell.differentiation-drive"
    public static let apoptosisDrive = "numivivo.cell.apoptosis-drive"
    public static let migrationDrive = "numivivo.cell.migration-drive"
    public static let matrixSynthesisRate = "numivivo.ecm.synthesis-rate"
    public static let matrixDegradationRate = "numivivo.ecm.degradation-rate"

    // NumiBrain -> NumiVivo
    public static let sympatheticTone = "numibrain.autonomic.sympathetic-tone"
    public static let parasympatheticTone = "numibrain.autonomic.parasympathetic-tone"
    public static let endocrineCommand = "numibrain.endocrine.command"
    public static let neuralActivity = "numibrain.local.activity"

    // NumiVivo -> NumiBrain
    public static let nociceptiveDrive = "numivivo.neural.nociceptive-drive"
    public static let inflammatoryDrive = "numivivo.neural.inflammatory-drive"
    public static let metabolicDrive = "numivivo.neural.metabolic-drive"
    public static let endocrineFeedback = "numivivo.endocrine.feedback"
}

public struct VivoNumiLabStandardChannelShape: Codable, Equatable, Sendable {
    public let spatialSampleCount: UInt32
    public let phenotypeCount: UInt32
    public let endocrineChannelCount: UInt32
    public let neuralChannelCount: UInt32

    public init(
        spatialSampleCount: UInt32,
        phenotypeCount: UInt32,
        endocrineChannelCount: UInt32,
        neuralChannelCount: UInt32
    ) {
        self.spatialSampleCount = spatialSampleCount
        self.phenotypeCount = phenotypeCount
        self.endocrineChannelCount = endocrineChannelCount
        self.neuralChannelCount = neuralChannelCount
    }

    public func validate() throws {
        guard spatialSampleCount > 0,
              phenotypeCount > 0,
              endocrineChannelCount > 0,
              neuralChannelCount > 0 else {
            throw VivoNumiLabCouplingError.invalidChannel("standard channel shapes must be non-zero")
        }
    }
}

/// Canonical channel declarations. All fields are numeric simulation contracts;
/// biological identity and ontology references remain in the CouplingPack.
public enum VivoNumiLabStandardChannels {
    public static func make(shape: VivoNumiLabStandardChannelShape) throws -> [VivoCouplingChannelDescriptor] {
        try shape.validate()
        let spatial = [shape.spatialSampleCount]
        let phenotypeSpatial = [shape.phenotypeCount, shape.spatialSampleCount]
        let endocrine = [shape.endocrineChannelCount, shape.spatialSampleCount]
        let neural = [shape.neuralChannelCount, shape.spatialSampleCount]

        return [
            .init(id: VivoNumiLabChannelID.tissuePressure, producer: .numanX, consumer: .numiVivo, unit: "Pa", semantics: .replace, shape: spatial),
            .init(id: VivoNumiLabChannelID.tissueTemperature, producer: .numanX, consumer: .numiVivo, unit: "K", semantics: .replace, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.tissueOxygen, producer: .numanX, consumer: .numiVivo, unit: "mol/m3", semantics: .replace, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.tissuePerfusion, producer: .numanX, consumer: .numiVivo, unit: "1/s", semantics: .replace, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.tissueShearRate, producer: .numanX, consumer: .numiVivo, unit: "1/s", semantics: .replace, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.tissueVolumetricStrain, producer: .numanX, consumer: .numiVivo, unit: "1", semantics: .replace, shape: spatial),

            .init(id: VivoNumiLabChannelID.osmoticPressure, producer: .numiVivo, consumer: .numanX, unit: "Pa", semantics: .add, shape: spatial),
            .init(id: VivoNumiLabChannelID.activeStress, producer: .numiVivo, consumer: .numanX, unit: "Pa", semantics: .add, shape: [6, shape.spatialSampleCount]),
            .init(id: VivoNumiLabChannelID.fluidSource, producer: .numiVivo, consumer: .numanX, unit: "1/s", semantics: .rate, shape: spatial),
            .init(id: VivoNumiLabChannelID.permeabilityModifier, producer: .numiVivo, consumer: .numanX, unit: "1", semantics: .replace, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.viscosityModifier, producer: .numiVivo, consumer: .numanX, unit: "1", semantics: .replace, shape: spatial, lowerBound: 0),

            .init(id: VivoNumiLabChannelID.cellVolumeFraction, producer: .numiTissue, consumer: .numiVivo, unit: "1", semantics: .replace, shape: phenotypeSpatial, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.membraneAreaDensity, producer: .numiTissue, consumer: .numiVivo, unit: "m2/m3", semantics: .replace, shape: phenotypeSpatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.extracellularVolumeFraction, producer: .numiTissue, consumer: .numiVivo, unit: "1", semantics: .replace, shape: spatial, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.phenotypeDistribution, producer: .numiTissue, consumer: .numiVivo, unit: "1", semantics: .replace, shape: phenotypeSpatial, lowerBound: 0, upperBound: 1),

            .init(id: VivoNumiLabChannelID.proliferationDrive, producer: .numiVivo, consumer: .numiTissue, unit: "1/s", semantics: .rate, shape: phenotypeSpatial),
            .init(id: VivoNumiLabChannelID.differentiationDrive, producer: .numiVivo, consumer: .numiTissue, unit: "1/s", semantics: .rate, shape: phenotypeSpatial),
            .init(id: VivoNumiLabChannelID.apoptosisDrive, producer: .numiVivo, consumer: .numiTissue, unit: "1/s", semantics: .rate, shape: phenotypeSpatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.migrationDrive, producer: .numiVivo, consumer: .numiTissue, unit: "m/s", semantics: .replace, shape: [3, shape.phenotypeCount, shape.spatialSampleCount]),
            .init(id: VivoNumiLabChannelID.matrixSynthesisRate, producer: .numiVivo, consumer: .numiTissue, unit: "kg/(m3*s)", semantics: .rate, shape: spatial, lowerBound: 0),
            .init(id: VivoNumiLabChannelID.matrixDegradationRate, producer: .numiVivo, consumer: .numiTissue, unit: "1/s", semantics: .rate, shape: spatial, lowerBound: 0),

            .init(id: VivoNumiLabChannelID.sympatheticTone, producer: .numiBrain, consumer: .numiVivo, unit: "1", semantics: .replace, shape: spatial, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.parasympatheticTone, producer: .numiBrain, consumer: .numiVivo, unit: "1", semantics: .replace, shape: spatial, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.endocrineCommand, producer: .numiBrain, consumer: .numiVivo, unit: "1", semantics: .replace, shape: endocrine),
            .init(id: VivoNumiLabChannelID.neuralActivity, producer: .numiBrain, consumer: .numiVivo, unit: "Hz", semantics: .replace, shape: neural, lowerBound: 0),

            .init(id: VivoNumiLabChannelID.nociceptiveDrive, producer: .numiVivo, consumer: .numiBrain, unit: "1", semantics: .replace, shape: neural, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.inflammatoryDrive, producer: .numiVivo, consumer: .numiBrain, unit: "1", semantics: .replace, shape: neural, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.metabolicDrive, producer: .numiVivo, consumer: .numiBrain, unit: "1", semantics: .replace, shape: neural, lowerBound: 0, upperBound: 1),
            .init(id: VivoNumiLabChannelID.endocrineFeedback, producer: .numiVivo, consumer: .numiBrain, unit: "1", semantics: .replace, shape: endocrine)
        ]
    }
}

/// Adapter for integrating a runtime without introducing a compile-time dependency
/// on the NumanX, NumiTissue or NumiBrain repositories. The closure owner retains
/// candidate state and must obey the two-phase publication contract.
public actor VivoClosureCoupledParticipant: VivoCoupledParticipant {
    public nonisolated let participantID: VivoNumiLabParticipant

    public typealias Prepare = @Sendable (VivoCouplingTransaction) async throws -> VivoPreparedCouplingState
    public typealias Stage = @Sendable (String) async throws -> Void
    public typealias Release = @Sendable (String) async -> Void
    public typealias Rollback = @Sendable (String) async -> Void

    private let prepareBody: Prepare
    private let stageBody: Stage
    private let releaseBody: Release
    private let rollbackBody: Rollback

    public init(
        participantID: VivoNumiLabParticipant,
        prepare: @escaping Prepare,
        stageCommit: @escaping Stage,
        releaseCommit: @escaping Release,
        rollback: @escaping Rollback
    ) {
        self.participantID = participantID
        self.prepareBody = prepare
        self.stageBody = stageCommit
        self.releaseBody = releaseCommit
        self.rollbackBody = rollback
    }

    public func prepare(_ transaction: VivoCouplingTransaction) async throws -> VivoPreparedCouplingState {
        let result = try await prepareBody(transaction)
        guard result.participant == participantID, result.transactionID == transaction.id else {
            throw VivoNumiLabCouplingError.participantMismatch
        }
        return result
    }

    public func stageCommit(transactionID: String) async throws {
        try await stageBody(transactionID)
    }

    public func releaseCommit(transactionID: String) async {
        await releaseBody(transactionID)
    }

    public func rollback(transactionID: String) async {
        await rollbackBody(transactionID)
    }
}
